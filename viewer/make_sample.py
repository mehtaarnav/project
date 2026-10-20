"""Synthetic, physically consistent sample run for the viewer (not a simulation result).

Kinematics follow the repo geometry exactly (belt path from the open-belt formulas of
cvtsim/cvt/geometry.py for a prescribed smooth x_p sweep, ratio 3.03 -> 0.80). Forces use the
spec formulas (flyweight force, spring, helix, traction capacity) at the prescribed state;
engine torque is the brief's chart. The system is NOT integrated, so the dynamics are not
equilibrium-consistent. Writes viewer/sample.json.

    python viewer/make_sample.py
"""
import json
import pathlib
import shutil
import tempfile

import numpy as np

from cvtsim.cvt.geometry import (ratio_from_xp, r_star_primary, secondary_radius_from_xp,
                                 belt_length)
from cvtsim.cvt.system import CH as CH0
from cvtsim.io import params
from cvtsim.viz.export import write_viewer

HERE = pathlib.Path(__file__).parent
T_END, F, N, NS = 8.0, 100, 160, 800
CH = CH0 + ["theta_e", "theta_s", "r_w", "F_fw", "cap_p", "cap_s", "T_belt", "v"]


def smooth(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def belt_path(rp, rs, C):
    """Arc-length parametrised closed belt path (ccw from the primary's upper tangent point).

    Returns f(s) -> (points (n,2), segment id 0..3: primary wrap, lower span, secondary wrap,
    upper span, fraction along the segment)."""
    phi = np.arcsin((rs - rp) / C)
    nu = np.array([-np.sin(phi), np.cos(phi)])
    nl = np.array([-np.sin(phi), -np.cos(phi)])
    Lsp = C * np.cos(phi)
    L1, L2, L3 = rp * (np.pi - 2 * phi), Lsp, rs * (np.pi + 2 * phi)
    cuts = np.cumsum([L1, L2, L3, Lsp])
    a0 = np.pi / 2 + phi
    ctr2 = np.array([C, 0.0])
    P1l, P2l = rp * nl, ctr2 + rs * nl
    P2u, P1u = ctr2 + rs * nu, rp * nu

    def f(s):
        s = np.mod(s, cuts[-1])
        seg = np.searchsorted(cuts, s, side="right")
        start = np.concatenate([[0.0], cuts[:-1]])[seg]
        ln = np.array([L1, L2, L3, Lsp])[seg]
        u = (s - start) / ln
        pts = np.empty((len(s), 2))
        m = seg == 0
        a = a0 + u[m] * (np.pi - 2 * phi)
        pts[m] = rp * np.column_stack([np.cos(a), np.sin(a)])
        m = seg == 1
        pts[m] = P1l + u[m, None] * (P2l - P1l)
        m = seg == 2
        a = -(np.pi / 2 + phi) + u[m] * (np.pi + 2 * phi)
        pts[m] = ctr2 + rs * np.column_stack([np.cos(a), np.sin(a)])
        m = seg == 3
        pts[m] = P2u + u[m, None] * (P1u - P2u)
        return pts, seg, u

    return f, cuts[-1], phi


def main():
    p = params.load()
    g = p["derived"]["geometry"]
    d = p["derived"]
    pr, sc, ve, bt = p["primary"], p["secondary"], p["vehicle"], p["belt"]
    tb, beta, L, C = g["tan_beta"], g["beta"], g["L"], g["C"]
    mu = bt["mu"]
    r_cam, p_h = sc["r_cam"], d["p_h"]

    def sweep(t):
        return g["x_p_stop"] * smooth((t - 0.6) / 6.4)

    tf = np.linspace(0.0, T_END, 4001)
    xp_f = sweep(tf)
    rp_f = r_star_primary(xp_f, g)
    rs_f = np.array([secondary_radius_from_xp(x, g) for x in xp_f])
    rpm_f = 1800.0 + 1300.0 * smooth(tf / 1.2) + 500.0 * smooth((tf - 1.2) / 5.0)
    rpm_f = np.minimum(rpm_f, 3600.0)
    w_e_f = rpm_f * 2 * np.pi / 60
    w_s_f = w_e_f * rp_f / rs_f
    theta_e_f = np.concatenate([[0], np.cumsum(0.5 * (w_e_f[1:] + w_e_f[:-1]) * np.diff(tf))])
    theta_s_f = np.concatenate([[0], np.cumsum(0.5 * (w_s_f[1:] + w_s_f[:-1]) * np.diff(tf))])
    s_belt_f = np.concatenate([[0], np.cumsum(0.5 * (w_e_f * rp_f)[1:] * np.diff(tf)
                                              + 0.5 * (w_e_f * rp_f)[:-1] * np.diff(tf))])

    ts = np.linspace(0.0, T_END, NS)
    I = lambda a: np.interp(ts, tf, a)
    xp, rp, rs, w_e, w_s = I(xp_f), I(rp_f), I(rs_f), I(w_e_f), I(w_s_f)
    xs = 2 * tb * (p["geometry"]["rs_mech_max"] - rs)
    u = xp + pr["engagement_gap"]
    ramp_k = pr["kappa"]
    r_w = pr["r_w0"] + ramp_k * u
    F_fw = pr["n_weights"] * d["m_w"] * r_w * ramp_k * w_e ** 2
    F_spr = pr["F_preload"] + pr["k_spring"] * u
    Fa_p = np.maximum(F_fw - F_spr, 0.0)
    cap_p = 2 * mu * rp * Fa_p / np.cos(beta)
    T_eng = np.interp(w_e * 60 / (2 * np.pi), p["engine"]["rpm"], p["engine"]["torque"])
    tq_p = np.minimum(T_eng, 0.9 * cap_p)
    tq_s = tq_p * rs / rp
    h = xs / p_h
    hp = 1.0 / p_h
    Fa_s = d["F_l0"] + sc["k_linear"] * xs + (tq_s / 2 + sc["k_torsion"] * (d["theta_t0"] + h)) * hp
    cap_s = 2 * mu * rs * Fa_s / np.cos(beta)
    n_arc_p = np.pi - 2 * np.arcsin((rs - rp) / C)
    R = np.exp(0.35 * mu * n_arc_p / np.sin(beta))
    Tc = bt["rho"] * bt["area"] * (w_e * rp) ** 2
    delta = tq_p / rp
    T_slack = delta / (R - 1) + Tc
    T_tight = R * delta / (R - 1) + Tc
    P_in = T_eng * w_e
    W_in = np.concatenate([[0], np.cumsum(0.5 * (P_in[1:] + P_in[:-1]) * np.diff(ts))])
    W_diss = 0.06 * W_in
    phi = np.arcsin((rs - rp) / C)
    d_low = np.column_stack([np.cos(phi), -np.sin(phi)])
    d_up = np.column_stack([np.cos(phi), np.sin(phi)])
    Fp = T_tight[:, None] * d_low + T_slack[:, None] * d_up
    Fs = -Fp
    cols = dict(
        t=ts, w_e=w_e, w_s=w_s, x_p=xp, xd_p=np.gradient(xp, ts), x_s=xs, xd_s=np.gradient(xs, ts),
        rstar_p=rp, rstar_s=rs, Fa_p=Fa_p, Fa_s=Fa_s, tq_p=tq_p, tq_s_fix=tq_s / 2,
        tq_s_mov=tq_s / 2, T_eng=T_eng, T_load=tq_s, Fx_p=Fp[:, 0], Fy_p=Fp[:, 1],
        Fx_s=Fs[:, 0], Fy_s=Fs[:, 1], T_tight=T_tight, T_slack=T_slack, T_max=T_tight * 1.04,
        r_eff_p=rp, r_eff_s=rs, slip_p=0.002 + 0.02 * (tq_p / cap_p) ** 8,
        slip_s=0.002 + 0.02 * (tq_s / cap_s) ** 8, E_mech=W_in - W_diss, W_in=W_in,
        W_diss=W_diss, Fstop_p=0 * ts, Fstop_s=0 * ts, Fcam=tq_s / 2 / r_cam,
        theta_e=I(theta_e_f), theta_s=I(theta_s_f), r_w=r_w, F_fw=F_fw, cap_p=cap_p,
        cap_s=cap_s, T_belt=30.0 + 28.0 * (1 - np.exp(-ts / 3.0)),
        v=w_s * ve["tire_radius"] / ve["final_drive"])
    res = np.column_stack([cols[c] for c in CH])

    # frames
    tfr = np.linspace(0.0, T_END, F)
    J = lambda a: np.interp(tfr, tf, a)
    xp_fr = sweep(tfr)
    rp_fr = r_star_primary(xp_fr, g)
    rs_fr = np.array([secondary_radius_from_xp(x, g) for x in xp_fr])
    belt = np.empty((F, N, 2))
    tension = np.empty((F, N))
    for i in range(F):
        path, Lp, ph = belt_path(rp_fr[i], rs_fr[i], C)
        assert abs(Lp - L) < 1e-9
        s = np.arange(N) * L / N + J(s_belt_f)[i]
        belt[i] = path(s)[0]
        pm, seg, uu = path(s + 0.5 * L / N)
        ii = np.argmin(np.abs(ts - tfr[i]))
        lnT = np.log(R[ii])
        Tt = np.where(seg == 0, uu, np.where(seg == 1, 1.0, np.where(seg == 2, 1 - uu, 0.0)))
        tension[i] = (T_slack[ii] - Tc[ii]) * np.exp(lnT * Tt) + Tc[ii]
    snaps = dict(t=tfr, belt=belt, tension=tension, x_p=xp_fr, x_s=2 * tb * (p["geometry"]["rs_mech_max"] - rs_fr),
                 theta_e=J(theta_e_f), theta_s=J(theta_s_f), r_w=pr["r_w0"] + ramp_k * (xp_fr + pr["engagement_gap"]))

    # bearing reactions from the overhung / between-bearing statics (07 A7.4 layout)
    shp, shs = p["shafts"]["primary"], p["shafts"]["secondary"]
    z1, z2 = shp["bearings"]
    zs = shp["sheave_z"]
    zb1, zb2 = shs["bearings"]
    zss = shs["sheave_z"]
    Fpm, Fsm = np.hypot(Fp[:, 0], Fp[:, 1]), np.hypot(Fs[:, 0], Fs[:, 1])
    sh = dict(t=ts, brg_p1=Fpm * abs(zs - z2) / (z2 - z1), brg_p2=Fpm * abs(zs - z1) / (z2 - z1),
              brg_s1=Fsm * abs(zb2 - zss) / (zb2 - zb1), brg_s2=Fsm * abs(zss - zb1) / (zb2 - zb1))

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp) / "sample"
        tmp.mkdir()
        out = write_viewer(tmp, res, CH, snaps, p, sh)
        shutil.copy(out, HERE / "sample.json")
    print("wrote", HERE / "sample.json", (HERE / "sample.json").stat().st_size // 1024, "kB")


if __name__ == "__main__":
    main()
