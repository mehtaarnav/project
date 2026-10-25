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

    fig = plt.figure(figsize=(16, 13))
    event = p.get("scenario_name", "launch")
    S.header(fig, f"cvtsim / {event}", f"{run_name}   ·   dt {p['sim']['dt']:.1e} s   ·   {int(p['sim']['belt_nodes'])} belt nodes   ·   * depends on shaft geometry, materials and bearing ratings: check python -m cvtsim.inputs --open")
    gs = fig.add_gridspec(5, 6, height_ratios=[0.42, 0.42, 1, 1, 1], hspace=0.55, wspace=0.45,
                          left=0.045, right=0.985, top=0.9, bottom=0.05)
    cards = [("engagement", f"{metrics['engage_rpm']:.0f}", "rpm", False),
             ("shift complete", "—" if np.isnan(metrics["t_shift_s"]) else f"{metrics['t_shift_s']:.2f}", "s", False),
             (f"speed @ {t[-1]:.1f} s", f"{metrics['v_end_kmh']:.1f}", "km/h", False),
             ("peak bearing*", f"{metrics['bearing_peak_N']:.0f}", "N", True),
             ("peak von mises*", f"{metrics['vm_peak_MPa']:.0f}", "MPa", True),
             ("energy residual / w_in", f"{metrics['energy_resid']:.0e}", "", False)]
    for i, (lab, val, unit, acc) in enumerate(cards):
        S.card(fig.add_subplot(gs[0, i]), lab, val, unit, acc)
    g = lambda k, fmt: fmt.format(metrics[k]) if np.isfinite(metrics.get(k, np.nan)) else "—"
    dur = [("PTO side load", g("pto_side_peak_N", "{:.0f}"), "N", True),
           ("PTO face moment", g("pto_face_moment_Nm", "{:.0f}"), "N·m", True),
           ("Goodman, crank PTO*", g("n_goodman_prim", "{:.2f}"), "", True),
           ("Goodman, jackshaft*", g("n_goodman_sec", "{:.2f}"), "", True),
           ("PTO main bearing P_eq, belt only", g("P_eq_prim2_N", "{:.0f}"), "N", True),
           ("jackshaft bearing C needed*", g("C_req_sec1_N", "{:.0f}"), "N", True)]
    for i, (lab, val, unit, acc) in enumerate(dur):
        S.card(fig.add_subplot(gs[1, i]), lab, val, unit, acc)

    def panel(r, c, name, xlabel="t [s]", ylabel=""):
        ax = fig.add_subplot(gs[r + 1, 2 * c:2 * c + 2])
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
    ax = panel(3, 1, "peak von mises (kt)*", ylabel="MPa")
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


# ---------------------------------------------------------------- engineering pages
# Channel contract: docs/theory/09_tuning_vehicle_belt.md section 6. A panel whose channels are
# absent is drawn as an empty frame naming them, so every page has a fixed layout.

FT = 0.3048
KMH = 3.6


def _get(res, CH, k):
    return res[:, CH.index(k)] if k in CH else None


def _have(CH, *keys):
    return all(k in CH for k in keys)


def _empty(ax, keys):
    ax.text(0.5, 0.5, "no channel: " + ", ".join(keys), ha="center", va="center",
            color=S.MUTED, fontsize=8, transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])


def _speed(res, CH, p):
    """Vehicle speed [m/s]: channel v, else the rigid-wheel value from w_s."""
    v = _get(res, CH, "v")
    if v is not None:
        return v
    return _col(res, CH, "w_s") * p["vehicle"]["tire_radius"] / p["vehicle"]["final_drive"]


def _distance(t, v):
    return np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(t))])


def _first_cross(t, y, level):
    """Time at which y first reaches level (linear interpolation); nan if never."""
    i = np.flatnonzero(y >= level)
    if i.size == 0:
        return np.nan
    j = i[0]
    return t[0] if j == 0 else t[j - 1] + (level - y[j - 1]) / (y[j] - y[j - 1]) * (t[j] - t[j - 1])


def _events(res, CH, p):
    """Indices of engagement (Fa_p > 1 N) and shift complete (x_p > 95 % of stop); None if absent."""
    eng = shift = None
    if "Fa_p" in CH and (_col(res, CH, "Fa_p") > 1.0).any():
        eng = int(np.argmax(_col(res, CH, "Fa_p") > 1.0))
    xs = p.get("derived", {}).get("geometry", {}).get("x_p_stop")
    if "x_p" in CH and xs and (_col(res, CH, "x_p") > 0.95 * xs).any():
        shift = int(np.argmax(_col(res, CH, "x_p") > 0.95 * xs))
    return eng, shift


def _margin(res, CH):
    """Belt torque over capacity at the primary and secondary sheaves."""
    tq_s = _col(res, CH, "tq_s_fix") + _col(res, CH, "tq_s_mov")
    cp, cs = _col(res, CH, "cap_p"), _col(res, CH, "cap_s")
    with np.errstate(divide="ignore", invalid="ignore"):
        mp = np.abs(_col(res, CH, "tq_p")) / np.where(cp > 1e-6, cp, np.nan)
        ms = np.abs(tq_s) / np.where(cs > 1e-6, cs, np.nan)
    return mp, ms


def key_metrics(res, CH, p):
    """Performance and clutch metrics for the cards and the comparison table (nan if unavailable)."""
    t = _col(res, CH, "t")
    v = _speed(res, CH, p)
    x = _distance(t, v)
    eng, shift = _events(res, CH, p)
    a = _get(res, CH, "a_x")
    m = dict(t_20=_first_cross(t, v, 20 / KMH), t_40=_first_cross(t, v, 40 / KMH),
             t_100ft=_first_cross(t, x, 100 * FT), t_150ft=_first_cross(t, x, 150 * FT),
             v_end_kmh=v[-1] * KMH, v_peak_kmh=v.max() * KMH,
             engage_rpm=_col(res, CH, "w_e")[eng] * 30 / np.pi if eng is not None else np.nan,
             t_shift=t[shift] if shift is not None else np.nan,
             a_peak=np.nanmax(np.gradient(v, t) if a is None else a))
    if _have(CH, "tq_p", "cap_p", "tq_s_fix", "tq_s_mov", "cap_s"):
        m["margin_max"] = float(np.nanmax(_margin(res, CH)))
    return m


def _grid(fig, nr, nc, **kw):
    kw = dict(hspace=0.55, wspace=0.38, left=0.06, right=0.975, top=0.88, bottom=0.07) | kw
    return fig.add_gridspec(nr, nc, **kw)


def _panel(fig, spec, name, xlabel="", ylabel=""):
    ax = fig.add_subplot(spec)
    S.title(ax, name)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    return ax


def engine(res, CH, p):
    """Engine operating trajectory on the full-throttle torque and power curves, coloured by time."""
    from cvtsim.drivetrain.engine import EngineCurve, RPM_TO_RAD
    S.use()
    t = _col(res, CH, "t")
    rpm = _col(res, CH, "w_e") / RPM_TO_RAD
    Te = _get(res, CH, "T_eng")
    curve = EngineCurve.from_params(p)
    r = np.linspace(min(rpm.min(), curve.rpm_min) * 0.95, curve.rpm_cut, 400)
    Tc = curve.torque(r * RPM_TO_RAD)
    eng, shift = _events(res, CH, p)

    fig = plt.figure(figsize=(16, 6.5))
    S.header(fig, "cvtsim / engine", "operating trajectory on the full-throttle curve   ·   colour = time")
    gs = _grid(fig, 1, 3, top=0.8, bottom=0.12)
    axes = [_panel(fig, gs[0, 0], "torque", "engine [rpm]", "N·m"),
            _panel(fig, gs[0, 1], "power", "engine [rpm]", "kW"),
            _panel(fig, gs[0, 2], "engine speed", "t [s]", "rpm")]
    lc = None
    for ax, y, yc in ((axes[0], Te, Tc),
                      (axes[1], None if Te is None else Te * rpm * RPM_TO_RAD / 1e3, Tc * r * RPM_TO_RAD / 1e3)):
        ax.plot(r, yc, color=S.MUTED, lw=1.0, ls=(0, (4, 2)), label="full throttle")
        ax.axvline(curve.rpm_cut, color=S.RED, lw=0.8, ls=(0, (2, 3)))
        ax.text(curve.rpm_cut, 0.02, "governor ", color=S.RED, fontsize=7, va="bottom", ha="right",
                rotation=90, transform=ax.get_xaxis_transform())
        if y is None:
            _empty(ax, ["T_eng"])
            continue
        lc = S.gradient_line(ax, rpm, y, t, lw=2.4)
        for idx, lab in ((eng, "engage"), (shift, "shift")):
            if idx is not None:
                ax.plot(rpm[idx], y[idx], "o", ms=6, color=S.WHITE, mec=S.RED, mew=1.2, zorder=5)
                ax.annotate(lab, (rpm[idx], y[idx]), xytext=(6, 8), textcoords="offset points",
                            fontsize=7.5, color=S.WHITE)
        ax.legend(loc="lower left")
    ax = axes[2]
    ax.plot(t, rpm)
    ax.axhline(curve.rpm_cut, color=S.RED, lw=0.8, ls=(0, (2, 3)))
    for idx in (eng, shift):
        if idx is not None:
            ax.axvline(t[idx], color=S.MUTED, lw=0.8, ls=(0, (2, 3)))
    if lc is not None:
        S.colorbar(fig, lc, axes[1], "t [s]")
    return fig


def clutch_forces(res, CH, p):
    """Axial force balance of the primary and secondary, against time and against travel."""
    S.use()
    t = _col(res, CH, "t")
    fig = plt.figure(figsize=(16, 9))
    S.header(fig, "cvtsim / clutch forces", "axial force balance   ·   belt thrust = clamping force on the belt")
    gs = _grid(fig, 2, 2)
    sets = {"p": [("F_fw", "flyweight", S.WHITE), ("F_spr_p", "spring", S.MUTED), ("Fa_p", "belt thrust", S.RED)],
            "s": [("F_spr_s", "linear spring", S.WHITE), ("F_tors", "torsion", "#9a9a9a"),
                  ("F_cam", "cam", S.MUTED), ("Fa_s", "belt thrust", S.RED)]}
    travel = {"p": "x_p", "s": "x_s"}
    for r, (key, name) in enumerate((("p", "primary"), ("s", "secondary"))):
        for c in (0, 1):
            xk = "t" if c == 0 else travel[key]
            ax = _panel(fig, gs[r, c], f"{name} force vs {'time' if c == 0 else 'travel'}",
                        "t [s]" if c == 0 else f"{name} travel [mm]", "N")
            drawn = [(lab, col, _get(res, CH, k)) for k, lab, col in sets[key] if k in CH]
            if xk not in CH or not drawn:
                _empty(ax, [xk] if xk not in CH else [k for k, _, _ in sets[key]])
                continue
            x = t if c == 0 else _col(res, CH, xk) * 1e3
            for lab, col, y in drawn:
                solid = lab in ("belt thrust", "flyweight", "linear spring")
                ax.plot(x, y, color=col, lw=1.8 if lab == "belt thrust" else 1.2,
                        ls="-" if solid else (0, (4, 2)), label=lab)
            ax.legend(ncol=2)
    return fig


def traction(res, CH, p):
    """Belt torque against traction capacity, slip power, belt temperature, tire force against mu N."""
    S.use()
    t = _col(res, CH, "t")
    fig = plt.figure(figsize=(16, 9))
    S.header(fig, "cvtsim / traction", "capacity 2 μ r F_a / cos β   ·   margin = torque / capacity, slip at 1")
    gs = _grid(fig, 2, 3)
    pk = ["tq_p", "cap_p"]
    sk = ["tq_s_fix", "tq_s_mov", "cap_s"]
    ax = _panel(fig, gs[0, 0], "primary belt torque", "t [s]", "N·m")
    if _have(CH, *pk):
        ax.plot(t, np.abs(_col(res, CH, "tq_p")), label="torque")
        ax.plot(t, _col(res, CH, "cap_p"), label="capacity"); ax.legend()
    else:
        _empty(ax, pk)
    ax = _panel(fig, gs[0, 1], "secondary belt torque", "t [s]", "N·m")
    if _have(CH, *sk):
        ax.plot(t, np.abs(_col(res, CH, "tq_s_fix") + _col(res, CH, "tq_s_mov")), label="torque")
        ax.plot(t, _col(res, CH, "cap_s"), label="capacity"); ax.legend()
    else:
        _empty(ax, sk)
    ax = _panel(fig, gs[0, 2], "traction margin", "t [s]", "torque / capacity")
    if _have(CH, *pk, *sk):
        mp, ms = _margin(res, CH)
        ax.plot(t, mp, label="primary"); ax.plot(t, ms, label="secondary")
        ax.axhline(1.0, color=S.RED, lw=0.8, ls=(0, (2, 3)))
        ax.set_ylim(0, max(1.2, min(2.0, np.nanmax([mp, ms]) * 1.1))); ax.legend()
    else:
        _empty(ax, pk + sk)
    ax = _panel(fig, gs[1, 0], "slip power", "t [s]", "W")
    if _have(CH, "P_slip_p", "P_slip_s"):
        ax.plot(t, _col(res, CH, "P_slip_p"), label="primary")
        ax.plot(t, _col(res, CH, "P_slip_s"), label="secondary"); ax.legend()
    else:
        _empty(ax, ["P_slip_p", "P_slip_s"])
    ax = _panel(fig, gs[1, 1], "belt temperature", "t [s]", "°C")
    if "T_belt" in CH:
        ax.plot(t, _col(res, CH, "T_belt") - 273.15)  # channel is in K
    else:
        _empty(ax, ["T_belt"])
    ax = _panel(fig, gs[1, 2], "tire force vs limit", "t [s]", "N")
    if _have(CH, "F_tire", "N_rear"):
        mu = p["vehicle"].get("mu_t", 0.7)  # ASSUMED, docs/theory/09 section 3
        ax.plot(t, _col(res, CH, "F_tire"), label="F_tire")
        ax.plot(t, mu * _col(res, CH, "N_rear"), color=S.RED, ls=(0, (4, 2)), label=f"μ N  (μ={mu:g})")
        ax.legend()
    else:
        _empty(ax, ["F_tire", "N_rear"])
    return fig


def performance(res, CH, p):
    """Speed and acceleration against time and distance, launch times, power to ground."""
    S.use()
    t = _col(res, CH, "t")
    v = _speed(res, CH, p)
    x = _distance(t, v)
    a = _get(res, CH, "a_x")
    a = np.gradient(v, t) if a is None else a
    m = key_metrics(res, CH, p)

    fig = plt.figure(figsize=(16, 10.5))
    S.header(fig, "cvtsim / performance", "full-throttle launch   ·   distance from integrated speed")
    gs = fig.add_gridspec(3, 3, height_ratios=[0.42, 1, 1], hspace=0.6, wspace=0.38,
                          left=0.06, right=0.975, top=0.88, bottom=0.07)

    def f(val, fmt="{:.2f}"):
        return "—" if np.isnan(val) else fmt.format(val)
    cards = [("0–20 km/h", f(m["t_20"]), "s"), ("0–40 km/h", f(m["t_40"]), "s"),
             ("100 ft", f(m["t_100ft"]), "s"), ("150 ft", f(m["t_150ft"]), "s"),
             ("peak speed", f(m["v_peak_kmh"], "{:.1f}"), "km/h"), ("peak accel", f(m["a_peak"], "{:.2f}"), "m/s²")]
    cg = gs[0, :].subgridspec(1, 6, wspace=0.1)
    for i, (lab, val, unit) in enumerate(cards):
        S.card(fig.add_subplot(cg[0, i]), lab, val, unit, accent=i in (1, 3))

    ax = _panel(fig, gs[1, 0], "speed vs time", "t [s]", "km/h"); ax.plot(t, v * KMH)
    ax = _panel(fig, gs[1, 1], "speed vs distance", "distance [m]", "km/h"); ax.plot(x, v * KMH)
    for d in (100 * FT, 150 * FT):
        if x[-1] >= d:
            ax.axvline(d, color=S.MUTED, lw=0.8, ls=(0, (2, 3)))
    ax = _panel(fig, gs[1, 2], "acceleration vs time", "t [s]", "m/s²")
    ax.plot(t, a)
    ax.axhline(0, color=S.GRID, lw=0.8)
    ax = _panel(fig, gs[2, 0], "acceleration vs distance", "distance [m]", "m/s²"); ax.plot(x, a)
    ax = _panel(fig, gs[2, 1], "power to ground", "t [s]", "kW")
    Te = _get(res, CH, "T_eng")
    Pe = None if Te is None else Te * _col(res, CH, "w_e") / 1e3
    Pg = _get(res, CH, "P_ground")
    if Pe is not None:
        ax.plot(t, Pe, color=S.MUTED, label="engine")
    if Pg is not None:
        ax.plot(t, Pg / 1e3, label="ground")
    if Pe is None and Pg is None:
        _empty(ax, ["T_eng", "P_ground"])
    else:
        ax.legend()
    ax = _panel(fig, gs[2, 2], "drivetrain efficiency", "t [s]", "P_ground / P_engine")
    if Pe is not None and Pg is not None:
        ok = Pe > 0.05 * Pe.max()
        ax.plot(t[ok], Pg[ok] / 1e3 / Pe[ok])
        ax.set_ylim(0, 1.1)
    else:
        _empty(ax, ["T_eng", "P_ground"])
    return fig


def compare(runs, CH, p):
    """Overlay of runs [(name, res)]: shift curves, speed, ratio, plus a key-metric table."""
    S.use()
    fig = plt.figure(figsize=(16, 9))
    S.header(fig, "cvtsim / compare", f"{len(runs)} runs   ·   first run drawn heavy   ·   best value per row in red")
    gs = _grid(fig, 2, 3, top=0.86, wspace=0.32, hspace=0.5, height_ratios=[1, 1.05])
    axs = [_panel(fig, gs[0, 0], "shift curve", "vehicle speed [km/h]", "engine [rpm]"),
           _panel(fig, gs[0, 1], "speed", "t [s]", "km/h"),
           _panel(fig, gs[0, 2], "ratio", "t [s]", "r_s / r_p")]
    dashes = ["-", (0, (4, 2)), (0, (1, 2))]
    mets = []
    for i, (name, res) in enumerate(runs):
        kw = dict(color=S.SERIES[i % len(S.SERIES)], ls=dashes[(i // len(S.SERIES)) % 3],
                  lw=1.9 if i == 0 else 1.1, label=name)
        t = _col(res, CH, "t")
        v = _speed(res, CH, p)
        axs[0].plot(v * KMH, _col(res, CH, "w_e") * 30 / np.pi, **kw)
        axs[1].plot(t, v * KMH, **kw)
        if _have(CH, "rstar_p", "rstar_s"):
            axs[2].plot(t, _col(res, CH, "rstar_s") / _col(res, CH, "rstar_p"), **kw)
        mets.append(key_metrics(res, CH, p))
    if not _have(CH, "rstar_p", "rstar_s"):
        _empty(axs[2], ["rstar_p", "rstar_s"])
    axs[0].legend(ncol=2)

    ax = fig.add_subplot(gs[1, :])
    ax.set_axis_off()
    S.title(ax, "key metrics")
    rows = [("engage [rpm]", "engage_rpm", "{:.0f}", None), ("shift complete [s]", "t_shift", "{:.2f}", "min"),
            ("0–20 km/h [s]", "t_20", "{:.2f}", "min"), ("0–40 km/h [s]", "t_40", "{:.2f}", "min"),
            ("100 ft [s]", "t_100ft", "{:.2f}", "min"), ("150 ft [s]", "t_150ft", "{:.2f}", "min"),
            ("end speed [km/h]", "v_end_kmh", "{:.1f}", "max"), ("peak accel [m/s²]", "a_peak", "{:.2f}", "max"),
            ("max torque/capacity", "margin_max", "{:.2f}", "min")]
    x0, dx = 0.22, 0.78 / max(len(runs), 1)
    y = 0.9
    for j, (name, _) in enumerate(runs):
        ax.text(x0 + dx * j, y, name, fontsize=8, color=S.MUTED, transform=ax.transAxes)
    step = 0.85 / (len(rows) + 1)
    for lab, key, fmt, best in rows:
        y -= step
        vals = np.array([m.get(key, np.nan) for m in mets], dtype=float)
        ax.text(0.0, y, lab, fontsize=8, color=S.MUTED, transform=ax.transAxes)
        ib = None
        if best and len(runs) > 1 and np.isfinite(vals).any():
            ib = int(np.nanargmin(vals) if best == "min" else np.nanargmax(vals))
        for j, val in enumerate(vals):
            ax.text(x0 + dx * j, y, "—" if np.isnan(val) else fmt.format(val), fontsize=9,
                    color=S.RED if j == ib else S.WHITE, transform=ax.transAxes)
        ax.plot([0, 1], [y - step * 0.35] * 2, color=S.GRID, lw=0.6, transform=ax.transAxes)
    return fig
