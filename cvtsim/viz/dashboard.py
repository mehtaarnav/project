"""Run dashboards in the cvtsim visual system."""
import numpy as np
import matplotlib.pyplot as plt

from cvtsim.viz import style as S


def _col(res, CH, k):
    return res[:, CH.index(k)]


def launch(res, CH, sh, p, metrics, run_name):
    """One-page launch dashboard: metric cards over nine time-history panels."""
    S.use()
    t = _col(res, CH, "t")
    rpm = _col(res, CH, "w_e") * 30 / np.pi
    ws = _col(res, CH, "w_s")
    v = (_col(res, CH, "v") if "v" in CH else ws * p["vehicle"]["tire_radius"] / p["vehicle"]["final_drive"]) * 3.6
    ratio = np.where(np.abs(ws) > 1.0, _col(res, CH, "w_e") / np.where(np.abs(ws) > 1.0, ws, 1.0), np.nan)
    geo = _col(res, CH, "rstar_s") / _col(res, CH, "rstar_p")

    fig = plt.figure(figsize=(16, 11.5))
    S.header(fig, "cvtsim / launch", f"{run_name}   ·   full throttle from idle   ·   dt {p['sim']['dt']:.1e} s   ·   {int(p['sim']['belt_nodes'])} belt nodes")
    gs = fig.add_gridspec(4, 6, height_ratios=[0.42, 1, 1, 1], hspace=0.55, wspace=0.45,
                          left=0.045, right=0.985, top=0.9, bottom=0.05)
    cards = [("engagement", f"{metrics['engage_rpm']:.0f}", "rpm", False),
             ("shift complete", "—" if np.isnan(metrics["t_shift_s"]) else f"{metrics['t_shift_s']:.2f}", "s", False),
             (f"speed @ {t[-1]:g} s", f"{metrics['v_end_kmh']:.1f}", "km/h", False),
             ("peak bearing", f"{metrics['bearing_peak_N']:.0f}", "N", True),
             ("peak von mises", f"{metrics['vm_peak_MPa']:.0f}", "MPa", True),
             ("energy residual / w_in", f"{metrics['energy_resid']:.0e}", "", False)]
    for i, (lab, val, unit, acc) in enumerate(cards):
        S.card(fig.add_subplot(gs[0, i]), lab, val, unit, acc)

    def panel(r, c, name, xlabel="t [s]", ylabel=""):
        ax = fig.add_subplot(gs[r, 2 * c:2 * c + 2])
        S.title(ax, name)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        return ax

    ax = panel(1, 0, "shift curve", "vehicle speed [km/h]", "engine [rpm]")
    ax.plot(v, rpm)
    ax.axhline(metrics["engage_rpm"], color=S.RED, lw=0.8, ls=(0, (2, 3)))
    ax = panel(1, 1, "ratio", ylabel="r_s / r_p")
    ax.plot(t, geo, label="geometric")
    ax.plot(t, ratio, label="speed w_e/w_s")
    ax.set_ylim(0.5, 3.5); ax.legend()
    ax = panel(1, 2, "sheave travel", ylabel="mm")
    ax.plot(t, _col(res, CH, "x_p") * 1e3, label="primary closing")
    ax.plot(t, _col(res, CH, "x_s") * 1e3, label="secondary opening"); ax.legend()
    ax = panel(2, 0, "belt span tension", ylabel="N")
    ax.plot(t, _col(res, CH, "T_tight"), label="tight")
    ax.plot(t, _col(res, CH, "T_slack"), label="slack"); ax.legend()
    ax = panel(2, 1, "belt slip", ylabel="m/s")
    ax.plot(t, _col(res, CH, "slip_p"), label="primary")
    ax.plot(t, _col(res, CH, "slip_s"), label="secondary"); ax.legend()
    ax = panel(2, 2, "bearing radial load", ylabel="N")
    for k in sh["bearing_cols"]:
        ax.plot(t, sh[k], lw=0.9, label=k.replace("F_bearing_", ""))
    ax.legend(ncol=2)
    ax = panel(3, 0, "shaft torque · bending", ylabel="N·m")
    ax.plot(t, np.abs(sh["T_prim"]), label="|T| primary")
    ax.plot(t, sh["M_prim"], label="M primary")
    ax.plot(t, np.abs(sh["T_sec"]), label="|T| secondary")
    ax.plot(t, sh["M_sec"], label="M secondary"); ax.legend(ncol=2)
    ax = panel(3, 1, "peak von mises (kt)", ylabel="MPa")
    ax.plot(t, sh["vm_prim"] / 1e6, label="primary shoulder")
    ax.plot(t, sh["vm_sec"] / 1e6, label="secondary shoulder"); ax.legend()
    ax = panel(3, 2, "energy closure", ylabel="J")
    E, Win, Wd = _col(res, CH, "E_mech"), _col(res, CH, "W_in"), _col(res, CH, "W_diss")
    ax.plot(t, Win / 1e3, label="input work [kJ]")
    ax.plot(t, ((E - E[0]) - (Win - Wd)), label="residual [J]"); ax.legend()
    return fig


def sweep(cases, CH, ref_name, metric_rows):
    """Overlay of variant runs (like comparing training runs) plus a convergence table."""
    S.use()
    fig = plt.figure(figsize=(16, 9))
    S.header(fig, "cvtsim / sweep", "convergence (time step, flank stiffness) and sensitivity (damping, flank angle)")
    gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.3, left=0.05, right=0.985, top=0.86, bottom=0.07)
    chans = [("w_e", "engine speed", "rpm", 30 / np.pi), ("x_p", "primary travel", "mm", 1e3),
             ("T_tight", "tight-span tension", "N", 1.0), ("Fa_s", "secondary thrust", "N", 1.0)]
    for i, (k, name, unit, scale) in enumerate(chans):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        S.title(ax, name)
        for j, (cname, res) in enumerate(cases):
            ref = cname == ref_name
            ax.plot(res[:, CH.index("t")], res[:, CH.index(k)] * scale,
                    color=S.WHITE if ref else [S.RED, "#9a9a9a", "#555555", "#3a3a3a"][j % 4],
                    lw=1.6 if ref else 0.8, alpha=1 if ref else 0.85, label=cname, zorder=3 if ref else 2)
        ax.set(xlabel="t [s]", ylabel=unit)
        if i == 0:
            ax.legend(ncol=2)
    ax = fig.add_subplot(gs[:, 2])
    ax.set_axis_off()
    S.title(ax, "end-of-window deltas vs baseline")
    base = dict(metric_rows)[ref_name]
    y = 0.95
    for cname, m in metric_rows:
        if cname == ref_name:
            continue
        ax.text(0.0, y, cname.upper(), fontsize=8, color=S.WHITE, transform=ax.transAxes)
        y -= 0.04
        for key in ("v_end", "x_s_end", "T_tight_max"):
            d = (m[key] - base[key]) / abs(base[key])
            ax.text(0.04, y, f"{key:<12s}{d:+.2%}", fontsize=7.5,
                    color=S.RED if abs(d) > 0.05 else S.MUTED, transform=ax.transAxes)
            y -= 0.032
        y -= 0.02
    return fig
