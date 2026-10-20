"""Verification dashboard: discretization error against refinement, plotted on log-log axes
the way a training curve plots loss against steps. Each panel shows the observed error and
the theoretical slope. Tracked as a run under runs/ (group "verification").

    python examples/verification.py      # about 2 min
"""
import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests"))
import test_contact as tc  # noqa: E402  (reuses the spec rigs; no duplicate physics)
import test_fe as tf  # noqa: E402
import test_integrators as ti  # noqa: E402

from cvtsim import runs  # noqa: E402
from cvtsim.contact.sheave import capstan_ratio  # noqa: E402
from cvtsim.io.params import load_raw  # noqa: E402
from cvtsim.viz import style as S  # noqa: E402


def plain_x(ax):
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())


def slope(ax, x, y0, order, label):
    x = np.asarray(x, float)
    ax.plot(x, y0 * (x / x[0]) ** order, color=S.MUTED, lw=0.8, ls=(0, (3, 3)), label=label)


def kirsch_scf(n):
    m, stress = tf._kirsch_solve(n, n)
    node = np.argmin(np.hypot(m.nodes[:, 0], m.nodes[:, 1] - tf.A_HOLE))
    return stress[node, 0]


def main():
    S.use()
    run = runs.Run("verification", load_raw(), group="verification")
    fig = plt.figure(figsize=(16, 9.5))
    S.header(fig, "cvtsim / verification", "discretization error vs refinement · dashed: theoretical order")
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3, left=0.05, right=0.985, top=0.86, bottom=0.08)

    # integrators: max trajectory error on x'' + w^2 x = 0
    ax = fig.add_subplot(gs[0, 0]); S.title(ax, "time integration")
    ns = np.array([50, 100, 200, 400, 800])
    errs = {"gen_alpha": [ti._genalpha_sdof(0.8, 1 / n) for n in ns],
            "verlet": [ti._verlet_sdof(1 / n) for n in ns],
            "symplectic_euler": [ti._euler_sdof(1 / n) for n in ns]}
    for name, e in errs.items():
        ax.loglog(ns, e, marker="o", ms=3, label=name.replace("_", " "))
        run.log(**{f"order_{name}": np.log2(e[-2] / e[-1])})
    slope(ax, ns, errs["verlet"][0], -2, "O(h^2)")
    slope(ax, ns, errs["symplectic_euler"][0], -1, "O(h)")
    plain_x(ax)
    ax.set(xlabel="steps per period", ylabel="max |x − cos ωt|"); ax.legend()

    # FE: Kirsch stress concentration
    ax = fig.add_subplot(gs[0, 1]); S.title(ax, "tet10 · kirsch scf")
    nk = np.array([4, 8, 16, 32])
    ek = [abs(kirsch_scf(n) - 3.0) / 3.0 for n in nk]
    ax.loglog(nk, ek, marker="o", ms=3, color=S.WHITE, label="|σ_xx(0,a) − 3σ∞| / 3σ∞")
    slope(ax, nk, ek[0], -2, "O(n^-2)")
    plain_x(ax)
    ax.set(xlabel="elements per 45° (n_r = n_θ)", ylabel="relative error"); ax.legend()
    run.log(kirsch_err_finest=ek[-1])

    # FE: cantilever tip
    ax = fig.add_subplot(gs[0, 2]); S.title(ax, "tet10 · cantilever tip")
    meshes = [(5, 1, 1), (10, 2, 2), (20, 2, 2), (40, 4, 4)]
    ec = [abs(tf._cantilever_tip(*m)[0] - 4.024) / 4.024 for m in meshes]
    nx = [m[0] for m in meshes]
    ax.loglog(nx, ec, marker="o", ms=3, color=S.WHITE, label="|δ − δ_Timoshenko| / δ")
    plain_x(ax)
    ax.set(xlabel="elements along span", ylabel="relative error"); ax.legend()
    run.log(cantilever_err_finest=ec[-1])

    # contact: penalty convergence of the capstan rig
    ax = fig.add_subplot(gs[1, 0]); S.title(ax, "contact · penalty stiffness")
    kn = np.array([1e5, 1e6, 1e7])
    target = capstan_ratio(0.3, np.pi, np.pi / 4)
    rr = [tc._rig_ratio_and_thrust(k) for k in kn]
    er = [abs(r[0] - target) / target for r in rr]
    ax.loglog(kn, er, marker="o", ms=3, color=S.WHITE, label="capstan ratio error")
    ax.set(xlabel="k_n [N/m]", ylabel="relative error"); ax.legend()
    run.log(capstan_err_k1e7=er[-1])

    # integrator stability map: spectral radius of generalized-α
    ax = fig.add_subplot(gs[1, 1]); S.title(ax, "gen-alpha · spectral radius")
    z = np.logspace(-2, 4, 200)
    for rho in (0.0, 0.5, 0.8, 1.0):
        sr = [np.max(np.abs(np.linalg.eigvals(ti.genalpha_amplification(v, rho)))) for v in z]
        ax.semilogx(z, sr, label=f"ρ∞ = {rho}")
    ax.set(xlabel="ω h", ylabel="spectral radius", ylim=(-0.05, 1.08)); ax.legend()

    # summary column
    ax = fig.add_subplot(gs[1, 2]); ax.set_axis_off(); S.title(ax, "observed orders")
    y = 0.85
    for k, v in run.metrics.items():
        ax.text(0.0, y, k.replace("_", " ").upper(), fontsize=8, color=S.MUTED, transform=ax.transAxes)
        ax.text(1.0, y, f"{v:.3g}", fontsize=13, fontfamily="Doto", fontweight=900, ha="right",
                color=S.WHITE, transform=ax.transAxes)
        y -= 0.12
    run.figure(fig, "verification")
    print(run.finish(), run.metrics)
    runs.write_index()


if __name__ == "__main__":
    main()
