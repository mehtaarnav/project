"""Full-throttle launch from idle through upshift (BRIEF Phase 2).

Outputs (examples/out/): launch.csv, shafts.csv, convergence.csv, energy.csv, *.png.
Usage: python examples/launch.py [--t-end 6] [--quick]
"""
import argparse
import pathlib
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cvtsim.io.params import load
from cvtsim.cvt.system import build, relax, simulate, CH
from cvtsim.fe.shaft_loads import shaft_histories

OUT = pathlib.Path(__file__).resolve().parent / "out"


def run(p, t_end, dt=None, nout=3000, tag=""):
    t0 = time.time()
    model = build(p, dt=dt)
    relax(model)
    res = simulate(model, t_end, nout=nout)
    print(f"[{tag}] dt={model['h']:.2e} steps={int(t_end/model['h'])} wall={time.time()-t0:.1f}s")
    return res, model


def col(res, name):
    return res[:, CH.index(name)]


def summarize(res, p):
    t = col(res, "t")
    we, ws = col(res, "w_e"), col(res, "w_s")
    v = ws * p["vehicle"]["tire_radius"] / p["vehicle"]["final_drive"]
    E, Win, Wd = col(res, "E_mech"), col(res, "W_in"), col(res, "W_diss")
    resid = (E - E[0]) - (Win - Wd)
    return dict(t=t, rpm=we * 30 / np.pi, v=v,
                ratio_speed=np.divide(we, ws, out=np.full_like(we, np.nan), where=np.abs(ws) > 1.0),
                ratio_geo=col(res, "rstar_s") / col(res, "rstar_p"),
                resid=resid, Win=Win)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-end", type=float, default=None)
    ap.add_argument("--quick", action="store_true", help="short run for smoke testing")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    p = load()
    t_end = a.t_end or (0.3 if a.quick else p["sim"]["t_end"])

    # ---- baseline
    res, model = run(p, t_end, tag="baseline")
    s = summarize(res, p)
    np.savetxt(OUT / "launch.csv", res, delimiter=",", header=",".join(CH), comments="")
    sh = shaft_histories(res, p, CH)
    np.savetxt(OUT / "shafts.csv", sh["table"], delimiter=",", header=",".join(sh["cols"]), comments="")

    # ---- energy balance
    rel = np.max(np.abs(s["resid"])) / max(np.max(np.abs(s["Win"])), 1e-12)
    print(f"energy residual max |dE-(Win-Wd)|/Win_final = {rel:.3e}")

    # ---- convergence: dt and penalty stiffness (shorter window covering engagement + early shift)
    t_conv = min(t_end, 2.0 if not a.quick else 0.2)
    conv_rows = []
    cases = [("dt", dict(dt=model["h"])), ("dt/2", dict(dt=model["h"] / 2)),
             ("kn x0.5", dict(kn_scale=0.5)), ("kn x2", dict(kn_scale=2.0))]
    for name, kw in cases:
        pp = load()
        m2 = build(pp, dt=kw.get("dt"), kn_scale=kw.get("kn_scale", 1.0))
        relax(m2)
        r2 = simulate(m2, t_conv, nout=500)
        s2 = summarize(r2, pp)
        conv_rows.append([name, s2["rpm"][-1], s2["v"][-1], col(r2, "x_p")[-1] * 1e3,
                          col(r2, "x_s")[-1] * 1e3, np.max(col(r2, "T_tight")),
                          np.max(np.hypot(col(r2, "Fx_s"), col(r2, "Fy_s")))])
        print("conv", conv_rows[-1])
    with open(OUT / "convergence.csv", "w") as f:
        f.write(f"# t_window={t_conv}\ncase,rpm_end,v_end,x_p_end_mm,x_s_end_mm,T_tight_max,F_hub_s_max\n")
        for r in conv_rows:
            f.write(",".join(str(x) if isinstance(x, str) else f"{x:.6g}" for x in r) + "\n")

    plots(res, s, sh, p)


def plots(res, s, sh, p):
    t = s["t"]
    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    ax = ax.ravel()
    ax[0].plot(s["v"] * 3.6, s["rpm"]); ax[0].set(xlabel="vehicle speed [km/h]", ylabel="engine rpm", title="Shift curve")
    ax[1].plot(t, s["ratio_speed"], label="speed ratio w_e/w_s"); ax[1].plot(t, s["ratio_geo"], "--", label="geometric r*_s/r*_p")
    ax[1].set(xlabel="t [s]", ylabel="ratio", ylim=(0, 4), title="Ratio"); ax[1].legend()
    ax[2].plot(t, col(res, "x_p") * 1e3, label="primary closing x_p"); ax[2].plot(t, col(res, "x_s") * 1e3, label="secondary opening x_s")
    ax[2].set(xlabel="t [s]", ylabel="mm", title="Sheave positions"); ax[2].legend()
    ax[3].plot(t, col(res, "T_tight"), label="upper span"); ax[3].plot(t, col(res, "T_slack"), label="lower span")
    ax[3].plot(t, col(res, "T_max"), ":", label="max element"); ax[3].set(xlabel="t [s]", ylabel="N", title="Belt tension"); ax[3].legend()
    ax[4].plot(t, col(res, "slip_p"), label="primary"); ax[4].plot(t, col(res, "slip_s"), label="secondary")
    ax[4].set(xlabel="t [s]", ylabel="m/s", title="Mean belt slip velocity (belt - sheave, tangential)"); ax[4].legend()
    for k in sh["bearing_cols"]:
        ax[5].plot(t, sh[k], label=k)
    ax[5].set(xlabel="t [s]", ylabel="N", title="Bearing radial loads"); ax[5].legend(fontsize=7)
    ax[6].plot(t, sh["T_prim"], label="primary torque"); ax[6].plot(t, sh["T_sec"], label="secondary torque")
    ax[6].plot(t, sh["M_prim"], "--", label="primary M @crit"); ax[6].plot(t, sh["M_sec"], "--", label="secondary M @crit")
    ax[6].set(xlabel="t [s]", ylabel="N·m", title="Shaft torque and bending"); ax[6].legend(fontsize=7)
    ax[7].plot(t, sh["vm_prim"] / 1e6, label="primary"); ax[7].plot(t, sh["vm_sec"] / 1e6, label="secondary")
    ax[7].set(xlabel="t [s]", ylabel="MPa", title="Peak von Mises at critical shoulder (Kt incl.)"); ax[7].legend()
    ax[8].plot(t, s["resid"], label="E-E0-(Win-Wdiss)"); ax[8].plot(t, s["Win"] * 1e-3, label="Win x1e-3")
    ax[8].set(xlabel="t [s]", ylabel="J", title="Energy balance"); ax[8].legend()
    for a_ in ax:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "launch.png", dpi=110)
    print("wrote", OUT / "launch.png")


if __name__ == "__main__":
    main()
