"""Convergence of peak loads in time step, friction smoothing and output rate (docs/CONVERGENCE.md).

    python examples/convergence_peaks.py [out_dir]

Runs 2.5 s launches at 40 kHz output with dt and belt.v_eps scaled, then prints a table of
peak loads relative to the default. Uses up to 4 processes.
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from cvtsim.cvt.system import CH, build, relax, simulate
from cvtsim.fe.shaft_loads import shaft_histories
from cvtsim.io.params import load
from cvtsim.post import durability as du

CASES = {"default": (1, 1), "dt/2": (0.5, 1), "v_eps x2": (1, 2), "v_eps x0.5": (1, 0.5),
         "v_eps x0.5, dt/2": (0.5, 0.5)}
T_END, NOUT = 2.5, 100001


def run(case):
    fdt, fve = CASES[case]
    p = load(tune="default")
    p["belt"]["v_eps"] *= fve
    m = build(p, dt=p["sim"]["dt"] * fdt)
    relax(m)
    return simulate(m, T_END, nout=NOUT, n_frames=2)[0]


def peaks(res, p):
    c = lambda k: res[:, CH.index(k)]
    sh = shaft_histories(res, p, CH)
    d = du.launch_durability(res, CH, sh, p)
    return {"PTO torque [N m]": np.abs(c("T_pto")).max(), "belt tension [N]": c("T_tight").max(),
            "PTO side load [N]": np.hypot(c("Fx_p"), c("Fy_p")).max(),
            "bearing load [N]": max(sh[k].max() for k in sh["bearing_cols"]),
            "jackshaft stress [MPa]": d["secondary"]["fatigue"]["vm_nominal_peak_MPa"],
            "crank PTO Goodman": d["primary"]["fatigue"]["n_goodman"],
            "jackshaft Goodman": d["secondary"]["fatigue"]["n_goodman"]}


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("convergence")
    out.mkdir(exist_ok=True)
    with ProcessPoolExecutor(4) as ex:
        results = dict(zip(CASES, ex.map(run, CASES)))
    p = load(tune="default")
    rows = {k: peaks(r, p) for k, r in results.items()}
    rows["default, 4 kHz"] = peaks(results["default"][::10], p)
    ref = rows["default"]
    print("| peak | " + " | ".join(rows) + " |\n|---|" + "---|" * len(rows))
    for k in ref:
        print(f"| {k} | " + " | ".join(f"{r[k]:.4g}" + ("" if r is ref else f" ({(r[k] - ref[k]) / abs(ref[k]):+.1%})")
                                        for r in rows.values()) + " |")
    for k, r in results.items():
        np.save(out / (k.replace("/", "").replace(" ", "_").replace(",", "") + ".npy"), r)


if __name__ == "__main__":
    main()
