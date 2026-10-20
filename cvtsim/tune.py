"""Tune runs and comparisons (docs/theory/09_tuning_vehicle_belt.md, section 5).

    python -m cvtsim.tune list
    python -m cvtsim.tune run default [--t-end 8] [--surface dirt]
    python -m cvtsim.tune compare default heavy_flyweights stiff_primary [--t-end 8]

Each simulation is a tracked run (group "tune") with config, series, snapshots and metrics;
`compare` adds a summary run holding a markdown table and an overlay figure.
"""
import argparse
from pathlib import Path

import numpy as np

from cvtsim import runs
from cvtsim.io.params import component, load, load_raw, load_tune_raw, tune_names
from cvtsim.cvt.system import build, relax, simulate, CH

METRIC_KEYS = ("engage_rpm", "t_shift_s", "v_end_kmh", "dist_m", "rpm_end", "x_s_end",
               "T_tight_max", "slip_s_mean", "T_belt_max_C", "tire_slip_max", "energy_resid")


OUT_DT = 2.5e-4  # s; resolves firing (<= 1.67 ms at 3600 rpm) and the PTO mode (<= 0.53 ms), 12 sec 7.6


def out_samples(t_end):
    """Output samples for a launch of length t_end that meets the durability sampling rule."""
    return int(np.ceil(t_end / OUT_DT)) + 1


def durability_metrics(res, p, sh):
    """Shaft, bearing and PTO side-load results (12 sec 7) as flat metrics."""
    from cvtsim.post.durability import launch_durability, summary_metrics
    return summary_metrics(launch_durability(res, CH, sh, p))


def _col(res, k):
    return res[:, CH.index(k)]


def metrics(res, p):
    """Scalar results of a launch (the numbers a tune comparison ranks)."""
    t = _col(res, "t")
    rpm = _col(res, "w_e") * 30 / np.pi
    v = _col(res, "v")
    E, Win, Wd = _col(res, "E_mech"), _col(res, "W_in"), _col(res, "W_diss")
    shifted = _col(res, "x_p") > 0.95 * p["derived"]["geometry"]["x_p_stop"]
    engaged = _col(res, "Fa_p") > 1.0
    v_w = _col(res, "w_wheel") * p["vehicle"]["tire_radius"]
    slip = (v_w - v) / np.maximum(np.abs(v), p["vehicle"]["v_min"])
    return dict(engage_rpm=rpm[np.argmax(engaged)] if engaged.any() else np.nan,
                t_shift_s=t[np.argmax(shifted)] if shifted.any() else np.nan,
                v_end=v[-1], v_end_kmh=3.6 * v[-1], dist_m=float(np.sum(0.5 * (v[1:] + v[:-1]) * np.diff(t))), rpm_end=rpm[-1],
                x_s_end=_col(res, "x_s")[-1] * 1e3, T_tight_max=_col(res, "T_tight").max(),
                slip_s_mean=np.mean(_col(res, "slip_s")[t > 0.5]) if t[-1] > 0.5 else np.nan,
                T_belt_max_C=_col(res, "T_belt").max() - 273.15, tire_slip_max=slip.max(),
                energy_resid=np.max(np.abs((E - E[0]) - (Win - Wd))) / max(Win[-1], 1e-12))


def run_tune(name, t_end, surface=None, dt=None, nout=None, n_frames=400, dashboard=True):
    """Simulate one tune as a tracked run. Returns (p, res, snaps, run, metrics)."""
    nout = nout or out_samples(t_end)
    overrides = {f"vehicle.{k}": v for k, v in component("vehicle_surfaces", surface).items()} if surface else {}
    p = load(tune=name, overrides=overrides)
    model = build(p, dt=dt)
    relax(model)
    res, snaps = simulate(model, t_end, nout=nout, n_frames=n_frames)
    cfg = dict(base=load_raw(), **load_tune_raw(name))
    run = runs.Run(f"tune-{Path(name).stem}" + (f"-{surface}" if surface else ""), cfg,
                   overrides=dict(overrides, dt=model["h"], t_end=t_end), group="tune", tags=(name,))
    run.series(CH, res)
    np.savez_compressed(run.dir / "snaps.npz", **snaps)
    m = metrics(res, p)
    run.log(**m)
    try:
        if dashboard:
            _dashboard(run, res, snaps, p, m)
            run.log(**m)  # adds shaft metrics computed by the dashboard
    finally:
        run.finish()  # the run record survives a plotting failure
    return p, res, snaps, run, m


def _dashboard(run, res, snaps, p, m):
    try:
        import matplotlib.pyplot as plt
        from cvtsim.fe.shaft_loads import shaft_histories
        from cvtsim.viz import dashboard
    except ImportError:
        return
    sh = shaft_histories(res, p, CH)
    m.update(bearing_peak_N=max(sh[k].max() for k in sh["bearing_cols"]),
             vm_peak_MPa=max(sh["vm_prim"].max(), sh["vm_sec"].max()) / 1e6)
    m.update(durability_metrics(res, p, sh))
    fig = dashboard.launch(res, CH, sh, p, m, run.dir.name)
    run.figure(fig, "dashboard")
    plt.close(fig)
    try:
        from cvtsim.viz.export import write_viewer
    except ImportError:
        return
    write_viewer(run.dir, res, CH, snaps, p, sh=sh)


def table(rows, ref):
    """Markdown table of metrics, with relative change against the reference tune."""
    lines = ["| tune | " + " | ".join(METRIC_KEYS) + " |", "|---|" + "---|" * len(METRIC_KEYS)]
    base = dict(rows)[ref]
    for name, m in rows:
        cells = []
        for k in METRIC_KEYS:
            cell = f"{m[k]:.4g}"
            if name != ref and np.isfinite(m[k]) and np.isfinite(base[k]) and base[k] != 0:
                cell += f" ({(m[k] - base[k]) / abs(base[k]):+.1%})"
            cells.append(cell)
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def compare(names, t_end, surface=None, **kw):
    """Run each tune (first is the reference) and write a summary run with a markdown table."""
    cases, rows = [], []
    for name in names:
        _, res, _, run, m = run_tune(name, t_end, surface=surface, **kw)
        cases.append((name, res))
        rows.append((name, m))
        print(name, run.dir.name, {k: float(f"{m[k]:.4g}") for k in METRIC_KEYS}, flush=True)
    summary = runs.Run("tune-compare-" + "-".join(names), dict(tunes=list(names), t_end=t_end, surface=surface),
                       group="tune")
    md = f"# Tune comparison ({t_end:g} s, reference `{names[0]}`)\n\n" + table(rows, names[0])
    (summary.dir / "compare.md").write_text(md)
    try:
        import matplotlib.pyplot as plt
        from cvtsim.viz import dashboard
        fig = getattr(dashboard, "compare", dashboard.sweep)(cases, CH, names[0], rows)
        summary.figure(fig, "compare")
        plt.close(fig)
    except ImportError:
        pass
    summary.finish()
    runs.write_index()
    return md


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.tune")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    for cmd in ("run", "compare"):
        sp = sub.add_parser(cmd)
        sp.add_argument("tunes", nargs="+" if cmd == "compare" else 1)
        sp.add_argument("--t-end", type=float, default=None)
        sp.add_argument("--surface", default=None)
        sp.add_argument("--dt", type=float, default=None)
        sp.add_argument("--nout", type=int, default=None, help="output samples (default: 4 kHz, the 12 sec 7.6 rule)")
        sp.add_argument("--frames", type=int, default=400)
    a = ap.parse_args(argv)
    if a.cmd == "list":
        print("\n".join(tune_names()))
        return
    t_end = a.t_end or load()["sim"]["t_end"]
    kw = dict(dt=a.dt, nout=a.nout, n_frames=a.frames)
    if a.cmd == "run":
        _, _, _, run, m = run_tune(a.tunes[0], t_end, surface=a.surface, **kw)
        runs.write_index()
        print(run.dir, {k: float(f"{m[k]:.4g}") for k in METRIC_KEYS})
    else:
        print(compare(a.tunes, t_end, surface=a.surface, **kw))


if __name__ == "__main__":
    main()
