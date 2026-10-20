"""Full-throttle launch from idle through upshift, tracked as runs.

    python examples/launch.py            # baseline (8 s) + sweep (8 cases, 2 s window)
    python examples/launch.py --quick    # 0.3 s smoke test, no sweep

Each simulation is a run under runs/ (see cvtsim/runs.py). Rebuilds runs/INDEX.md.
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt

from cvtsim import runs
from cvtsim.io.params import load, load_raw
from cvtsim.cvt.system import build, relax, simulate, CH
from cvtsim.fe.shaft_loads import shaft_histories
from cvtsim.viz import dashboard


def col(res, k):
    return res[:, CH.index(k)]


def metrics_of(res, p, sh=None):
    t = col(res, "t")
    rpm = col(res, "w_e") * 30 / np.pi
    v = col(res, "w_s") * p["vehicle"]["tire_radius"] / p["vehicle"]["final_drive"]
    E, Win, Wd = col(res, "E_mech"), col(res, "W_in"), col(res, "W_diss")
    x_stop = p["derived"]["geometry"]["x_p_stop"]
    shifted = col(res, "x_p") > 0.95 * x_stop
    m = dict(engage_rpm=rpm[np.argmax(col(res, "Fa_p") > 1.0)],
             t_shift_s=t[np.argmax(shifted)] if shifted.any() else np.nan,
             v_end=v[-1], v_end_kmh=3.6 * v[-1], rpm_end=rpm[-1],
             x_s_end=col(res, "x_s")[-1] * 1e3,
             T_tight_max=col(res, "T_tight").max(),
             slip_s_mean=np.mean(col(res, "slip_s")[t > 0.5]) if t[-1] > 0.5 else np.nan,
             energy_resid=np.max(np.abs((E - E[0]) - (Win - Wd))) / max(Win[-1], 1e-12))
    if sh is not None:
        m.update(bearing_peak_N=max(sh[k].max() for k in sh["bearing_cols"]),
                 vm_peak_MPa=max(sh["vm_prim"].max(), sh["vm_sec"].max()) / 1e6,
                 M_prim_peak=sh["M_prim"].max(), M_sec_peak=sh["M_sec"].max())
    return m


def simulate_run(name, t_end, overrides=None, dt=None, kn_scale=1.0, group="", nout=3000):
    p = load(overrides=overrides)
    model = build(p, dt=dt, kn_scale=kn_scale)
    relax(model)
    res = simulate(model, t_end, nout=nout)
    run = runs.Run(name, load_raw(), overrides=dict(overrides or {}, dt=model["h"], kn_scale=kn_scale),
                   group=group)
    run.series(CH, res)
    return p, res, run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    p0 = load()
    t_end = 0.3 if a.quick else p0["sim"]["t_end"]

    p, res, run = simulate_run("launch", t_end, group="baseline")
    sh = shaft_histories(res, p, CH)
    m = metrics_of(res, p, sh)
    run.log(**m)
    fig = dashboard.launch(res, CH, sh, p, m, run.dir.name)
    run.figure(fig, "dashboard")
    plt.close(fig)
    print("baseline", run.finish(), {k: round(v, 4) for k, v in m.items()}, flush=True)

    if not a.quick:
        t_var = 2.0
        dt0 = p0["sim"]["dt"]
        variants = [("baseline", {}, {}), ("dt/2", {"dt": dt0 / 2}, {}),
                    ("k_n x0.5", {"kn_scale": 0.5}, {}), ("k_n x2", {"kn_scale": 2.0}, {}),
                    ("k_n x4", {"kn_scale": 4.0}, {}),
                    ("damping x0.5", {}, {"primary.c_axial": 5.0, "secondary.c_axial": 5.0}),
                    ("damping x2", {}, {"primary.c_axial": 20.0, "secondary.c_axial": 20.0}),
                    ("beta 11.64 deg", {}, {"sheave.beta_deg": 11.64})]
        cases, rows = [], []
        for name, kw, ov in variants:
            pv, rv, rr = simulate_run(f"sweep-{name}", t_var, overrides=ov, group="sweep", nout=600, **kw)
            mv = metrics_of(rv, pv)
            rr.log(**mv)
            rr.finish()
            cases.append((name, rv))
            rows.append((name, mv))
            print("sweep", name, {k: round(v, 4) for k, v in mv.items()}, flush=True)
        fig = dashboard.sweep(cases, CH, "baseline", rows)
        summary = runs.Run("sweep-summary", load_raw(), group="sweep")
        summary.figure(fig, "sweep")
        summary.finish()
        plt.close(fig)
    runs.write_index()


if __name__ == "__main__":
    main()
