"""Baja events on the fast model: acceleration, hill climb and endurance (docs/theory/16).

    python -m cvtsim.events run default --event acceleration
    python -m cvtsim.events compare default cvtech_stock cvtech_0651_1100 --event endurance

A course is a list of segments [length, grade start, grade end, corner radius, traction]. The
driver follows a target speed set by corner grip and by how fast the car can slow without brakes.
"""
from __future__ import annotations

import argparse
import math
import pathlib

import numpy as np
import yaml

from cvtsim.drivetrain.vehicle import G_ACC
from cvtsim.fast.model import FCH, build_fast, simulate_fast
from cvtsim.io.params import load

TRACKS = pathlib.Path(__file__).resolve().parent / "io" / "tracks"
DS = 1.0             # course sampling [m]
V_FREE = 100.0       # target on straights: no limit [m/s]


def track_names():
    return sorted(q.stem for q in TRACKS.glob("*.yaml"))


def load_track(name):
    q = pathlib.Path(name)
    q = q if q.suffix in (".yaml", ".yml") and q.exists() else TRACKS / f"{name}.yaml"
    d = yaml.safe_load(q.read_text())
    d["name"] = q.stem
    return d


def course(track, laps=None):
    """Course sampled every DS: arrays s, grade [rad], traction, corner radius (0 = straight) (16 s1)."""
    seg = np.array(track["segments"], float)
    if seg.ndim != 2 or seg.shape[1] != 5 or np.any(seg[:, 0] <= 0) or np.any(seg[:, 4] < 0):
        raise ValueError("segments must be rows [length > 0, grade0, grade1, radius, traction >= 0]")
    s, g, tr, R = [0.0], [math.radians(seg[0, 1])], [seg[0, 4]], [seg[0, 3]]
    for L, g0, g1, rad, trc in seg:
        n = max(int(math.ceil(L / DS)), 1)
        for i in range(1, n + 1):
            s.append(s[-1] + L / n)
            g.append(math.radians(g0 + (g1 - g0) * i / n))
            tr.append(trc)
            R.append(rad)
    s, g, tr, R = map(np.array, (s, g, tr, R))
    n_laps = int(laps or track.get("laps", 1))
    if n_laps > 1:
        L = s[-1]
        s = np.concatenate([s] + [s[1:] + k * L for k in range(1, n_laps)])
        g, tr, R = (np.concatenate([a] + [a[1:]] * (n_laps - 1)) for a in (g, tr, R))
    return s, g, tr, R


def target_speed(s, g, tr, R, mu_lat, a_coast):
    """Driver target: corner limit sqrt(mu g R), reached by lifting off (16 s2).

    Lift-off deceleration on a grade g is a_coast + g sin(g); where that is not positive (steep
    descent) the driver cannot slow, and the target is left at the corner limit downstream.
    """
    v = np.where(R > 0, np.sqrt(np.maximum(mu_lat * tr, 0.0) * G_ACC * np.maximum(R, 1e-9)), V_FREE)
    for i in range(len(s) - 2, -1, -1):          # backward pass: lift early enough to make the corner
        a = max(a_coast + G_ACC * math.sin(g[i]), 0.0)
        v[i] = min(v[i], math.sqrt(v[i + 1] ** 2 + 2.0 * a * (s[i + 1] - s[i])))
    return v


def run_event(tune="default", event="acceleration", car=None, laps=None, cal=None, out_dt=0.05, params=None):
    """Simulate one event; returns (metrics dict, result array, course dict).

    params: a loaded parameter dict to use instead of loading `tune` (e.g. with overrides)."""
    p = params if params is not None else (load(car, tune=tune) if car else load(tune=tune))
    tk = load_track(event)
    s, g, tr, R = course(tk, laps)
    drv = p["driver"]
    v_t = target_speed(s, g, tr, R, p["vehicle"]["mu_t"], drv["a_coast"])
    trk = np.column_stack([s, g, tr, v_t])
    n_laps = int(laps or tk.get("laps", 1))
    t_lim = float(tk["t_limit"]) * n_laps
    fm = build_fast(p, cal=cal)
    res = simulate_fast(fm, t_lim, nout=int(t_lim / out_dt) + 1,
                        track=dict(trk=trk, gain=drv["gain"], lag=drv["throttle_lag"], s_end=s[-1]))
    return metrics(res, s[-1], n_laps, R, s, v_t, p["vehicle"]["tire_radius"]), res, dict(s=s, grade=g, traction=tr, radius=R, v_target=v_t)


def metrics(res, s_end, n_laps, R, s, v_t, r_t):
    c = lambda k: res[:, FCH.index(k)]
    t, v, rpm = c("t"), c("v"), c("w_e") * 30 / math.pi
    done = c("s")[-1] >= s_end - 1e-6
    rp = c("r_p")
    ok = rp > 0
    ratio = c("r_s")[ok] / rp[ok]
    corner = R > 0
    v_at = np.interp(s[corner], c("s"), v) if corner.any() else np.array([0.0])
    over = float(np.max(v_at - v_t[corner])) if corner.any() else 0.0
    moving = v > 0.5
    eng = c("W_in")[-1]
    return dict(
        finished=bool(done), time_s=float(t[-1]) if done else math.nan,
        lap_time_s=float(t[-1]) / n_laps if done else math.nan,
        distance_m=float(c("s").max()), stalled=bool(not done and v[-1] < 0.1),
        v_avg_kmh=float(3.6 * c("s").max() / t[-1]),
        v_max_kmh=float(3.6 * v.max()), rpm_mean=float(rpm[moving].mean()) if moving.any() else math.nan,
        rpm_p10=float(np.percentile(rpm[moving], 10)) if moving.any() else math.nan,
        rpm_p90=float(np.percentile(rpm[moving], 90)) if moving.any() else math.nan,
        belt_slip_time_frac=float(np.mean(c("stick")[moving] < 0.5)) if moving.any() else math.nan,
        belt_T_max_C=float(c("T_belt").max() - 273.15), engine_work_kJ=float(eng / 1e3),
        corner_overspeed_kmh=float(3.6 * over),
        wheelspin_time_frac=float(np.mean((c("w_wheel") * r_t - v)[moving] > 0.2 * np.maximum(v[moving], 1.0)))
        if moving.any() else math.nan,
        ratio_min=float(ratio.min()), ratio_max=float(ratio.max()))


def _fmt(m):
    return {k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.events", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list events and tunes")
    for cmd in ("run", "compare"):
        sp = sub.add_parser(cmd, help="simulate one tune" if cmd == "run" else "simulate several tunes")
        sp.add_argument("tunes", nargs=1 if cmd == "run" else "+")
        sp.add_argument("--event", default="acceleration", help="acceleration, hill_climb, endurance or a track file")
        sp.add_argument("--laps", type=int, default=None)
        sp.add_argument("--car", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "list":
        from cvtsim.io.params import tune_names
        print("events:", " ".join(track_names()))
        print("tunes:", " ".join(tune_names()))
        return
    for tune in a.tunes:
        m, _, _ = run_event(tune, a.event, a.car, a.laps)
        print(tune, _fmt(m))


if __name__ == "__main__":
    main()
