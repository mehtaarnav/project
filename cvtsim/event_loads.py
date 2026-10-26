"""Loads and durability over a whole event, and full-model replays of its worst moments (17).

    python -m cvtsim.event_loads run default --event endurance --laps 106 --replays 3
    python -m cvtsim.event_loads sweep default cvtech_stock cvtech_0651_1100 --event endurance --final-drive 8.1 8.96

The fast model gives the mean belt loads over the whole event. That covers shaft bending,
bearing life and the PTO side load. It has no firing pulses or PTO twist, so peak torque comes
from the full model, replayed over the windows where the fast model's loads peak.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from cvtsim.cvt.system import CH
from cvtsim.events import run_event
from cvtsim.fast.model import FCH
from cvtsim.fe.shaft_loads import shaft_histories
from cvtsim.io.params import load
from cvtsim.post import durability as du

WARMUP = 0.2      # s discarded after a warm start (17 section 2)
WINDOW = 1.0      # s replayed around each peak
LAUNCH_T = 3.0    # s: a peak this early is replayed from rest


def hub_forces(Tm, F_b, r_p, r_s, C):
    """Belt force on the primary hub (+x towards the secondary) from mean span tension Tm
    (centrifugal part removed) and belt pull F_b (17 section 1)."""
    gam = np.arcsin(np.clip((r_s - r_p) / C, -1.0, 1.0))
    return 2.0 * Tm * np.cos(gam), F_b * np.sin(gam)


def fast_to_full_channels(res, p):
    """A full-model channel array (CH) filled from fast-model output, for the shaft model."""
    c = lambda k: res[:, FCH.index(k)]
    out = np.zeros((res.shape[0], len(CH)))
    put = lambda k, v: out.__setitem__((slice(None), CH.index(k)), v)
    Fx, Fy = hub_forces(c("T_mean"), c("F_b"), c("r_p"), c("r_s"), p["derived"]["geometry"]["C"])
    put("t", c("t")); put("w_e", c("w_e")); put("w_c", c("w_e")); put("w_s", c("w_s"))
    put("Fx_p", Fx); put("Fy_p", Fy); put("Fx_s", -Fx); put("Fy_s", -Fy)
    put("T_pto", c("F_b") * c("r_p")); put("tq_p", -c("F_b") * c("r_p"))
    put("tq_s_fix", 0.5 * c("F_b") * c("r_s")); put("tq_s_mov", 0.5 * c("F_b") * c("r_s"))
    put("v", c("v")); put("x_p", c("x_p")); put("x_s", c("x_s"))
    return out


def event_durability(res, p, hours=None):
    """Durability of one event from the fast model (mean belt loads) (17 section 1).

    hours: duty to scale per-event damage to (e.g. 4 for the endurance race)."""
    full = fast_to_full_channels(res, p)
    sh = shaft_histories(full, p, CH)
    d = du.launch_durability(full, CH, sh, p)
    T = float(res[-1, 0])
    m = du.summary_metrics(d)
    m.pop("durability_sampled_ok")   # the firing-resolution rule applies to the full model only
    m.update(event_time_s=T, pto_torque_peak_mean_model_Nm=float(np.abs(full[:, CH.index("T_pto")]).max()))
    for name, tag in (("primary", "prim"), ("secondary", "sec")):
        for j, b in enumerate(d[name]["bearings"]):
            m[f"bearing_{tag}{j+1}_rpm_mean"] = b["rpm_mean"]
            if "L10_h" in b:
                m[f"L10_h_{tag}{j+1}"] = b["L10_h"]
    if hours:
        k = hours * 3600.0 / T
        m["duty_hours"] = hours
        for tag in ("prim", "sec"):
            m[f"damage_{tag}_per_duty"] = m[f"damage_{tag}"] * k
    return m, sh


def pick_windows(t, load_, k, sep=3.0):
    """Times of the k largest local peaks of load_, at least sep apart (17 section 2)."""
    order = np.argsort(load_)[::-1]
    picked = []
    for i in order:
        if all(abs(t[i] - t[j]) >= sep for j in picked):
            picked.append(i)
        if len(picked) == k:
            break
    return [float(t[i]) for i in picked]


def _replay(args):
    tune, car, t_peak, state, knots, t_end_launch = args
    from cvtsim.cvt.system import build, relax, simulate, warm_start
    p = load(car, tune=tune) if car else load(tune=tune)
    m = build(p, scenario=dict(knots=knots))
    if state is None:      # peak near the start: replay the launch from rest
        relax(m)
        res, _ = simulate(m, t_end_launch, nout=int(t_end_launch / 2.5e-4) + 1, n_frames=2)
        keep = res[:, CH.index("t")] >= 0.0
    else:
        warm_start(m, **state)
        res, _ = simulate(m, WINDOW, nout=int(WINDOW / 2.5e-4) + 1, n_frames=2)
        keep = res[:, CH.index("t")] >= state["t0"] + WARMUP
    res = res[keep]
    sh = shaft_histories(res, p, CH)
    d = du.launch_durability(res, CH, sh, p)
    c = lambda k: res[:, CH.index(k)]
    return dict(t_peak=t_peak, from_rest=state is None,
                pto_torque_peak_Nm=float(np.abs(c("T_pto")).max()),
                pto_side_peak_N=float(np.hypot(c("Fx_p"), c("Fy_p")).max()),
                belt_tension_peak_N=float(c("T_tight").max()),
                bearing_peak_N=float(max(sh[k].max() for k in sh["bearing_cols"])),
                n_goodman_prim=d["primary"]["fatigue"]["n_goodman"],
                n_goodman_sec=d["secondary"]["fatigue"]["n_goodman"],
                n_yield_prim=d["primary"]["fatigue"]["n_yield"],
                n_yield_sec=d["secondary"]["fatigue"]["n_yield"])


def replays(res, p, tune, car=None, k=3, workers=3):
    """Full-model replays of the k windows where the fast model's PTO side load peaks (17 s2)."""
    c = lambda key: res[:, FCH.index(key)]
    t = c("t")
    Fx, Fy = hub_forces(c("T_mean"), c("F_b"), c("r_p"), c("r_s"), p["derived"]["geometry"]["C"])
    peaks = pick_windows(t, np.hypot(Fx, Fy), k)
    jobs = []
    for tp in peaks:
        lo = max(tp - 0.5 * (WINDOW - WARMUP) - WARMUP, 0.0)
        hi = (LAUNCH_T + 0.5) if tp < LAUNCH_T else lo + WINDOW + 0.05
        sel = (t >= 0.0 if tp < LAUNCH_T else t >= lo - 0.05) & (t <= hi)
        knots = [[float(a), float(b), float(math.degrees(g)), 1.0] for a, b, g in
                 zip(t[sel][::2], c("throttle")[sel][::2], c("grade")[sel][::2])]
        if len(knots) < 2:
            knots = [[0.0, 1.0, 0.0, 1.0], [hi, 1.0, 0.0, 1.0]]
        if tp < LAUNCH_T:
            state = None
        else:
            i = int(np.argmin(abs(t - lo)))
            state = dict(w_e=float(c("w_e")[i]), x_p=float(c("x_p")[i]), x_s=float(c("x_s")[i]),
                         w_s=float(c("w_s")[i]), v=float(c("v")[i]), T_belt=float(c("T_belt")[i]),
                         t0=float(t[i]))
        jobs.append((tune, car, tp, state, knots, hi))
    with ProcessPoolExecutor(min(workers, len(jobs) or 1)) as ex:
        return list(ex.map(_replay, jobs))


def run(tune, event, laps=None, car=None, n_replays=3, hours=None, workers=3):
    """Event metrics, event durability and replays for one set-up."""
    p = load(car, tune=tune) if car else load(tune=tune)
    m, res, _ = run_event(tune, event, car, laps)
    dur, _ = event_durability(res, p, hours)
    rep = replays(res, p, tune, car, n_replays, workers) if n_replays else []
    return dict(tune=tune, event=event, metrics=m, durability=dur, replays=rep)


def sweep(tunes, event, final_drives=(None,), laps=None, car=None, hours=None):
    """Event metrics and durability (no replays) for every tune x final drive (17 section 3)."""
    rows = []
    for tune in tunes:
        for fd in final_drives:
            ov = {"vehicle.final_drive": fd} if fd else None
            p = load(car, tune=tune, overrides=ov) if car else load(tune=tune, overrides=ov)
            m, res, _ = run_event(tune, event, car, laps, params=p)
            dur, _ = event_durability(res, p, hours)
            rows.append(dict(tune=tune, final_drive=p["vehicle"]["final_drive"], **m,
                             pto_side_peak_N=dur["pto_side_peak_N"], P_eq_prim2_N=dur["P_eq_prim2_N"],
                             damage_sec=dur["damage_sec"], n_goodman_sec=dur["n_goodman_sec"]))
    return rows


def _table(rows, keys):
    lines = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    for r in rows:
        lines.append("| " + " | ".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + " |")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.event_loads", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="loads, durability and full-model replays for one set-up")
    r.add_argument("tune")
    s = sub.add_parser("sweep", help="compare set-ups (and final drives) on one event")
    s.add_argument("tunes", nargs="+")
    s.add_argument("--final-drive", nargs="+", type=float, default=[None])
    s.add_argument("--csv", default=None)
    for q in (r, s):
        q.add_argument("--event", default="endurance")
        q.add_argument("--laps", type=int, default=None)
        q.add_argument("--car", default=None)
        q.add_argument("--hours", type=float, default=None, help="duty to scale damage to, e.g. 4")
    r.add_argument("--replays", type=int, default=3)
    r.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    from cvtsim import runs
    if a.cmd == "run":
        out = run(a.tune, a.event, a.laps, a.car, a.replays, a.hours)
        rec = runs.Run(f"event-{a.tune}-{a.event}", dict(tune=a.tune, event=a.event, laps=a.laps, car=a.car,
                                                          replays=a.replays, hours=a.hours), group="event")
        rec.log(**{k: v for k, v in out["metrics"].items() if isinstance(v, (int, float))})
        (rec.dir / "event.json").write_text(json.dumps(out, indent=1, default=float))
        rec.finish()
        runs.write_index()
        if a.json:
            Path(a.json).write_text(json.dumps(out, indent=1, default=float))
        print(rec.dir)
        print(json.dumps(out, indent=1, default=float))
    else:
        rows = sweep(a.tunes, a.event, a.final_drive, a.laps, a.car, a.hours)
        keys = ["tune", "final_drive", "time_s", "v_avg_kmh", "rpm_mean", "belt_slip_time_frac",
                "corner_overspeed_kmh", "pto_side_peak_N", "P_eq_prim2_N", "n_goodman_sec"]
        md = _table(rows, keys)
        print(md)
        rec = runs.Run(f"sweep-{a.event}", dict(tunes=a.tunes, final_drives=a.final_drive, event=a.event,
                                                 laps=a.laps, car=a.car), group="event")
        (rec.dir / "sweep.md").write_text(md + "\n")
        with open(rec.dir / "sweep.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        rec.finish()
        runs.write_index()
        print(rec.dir)
        if a.csv:
            with open(a.csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    main()
