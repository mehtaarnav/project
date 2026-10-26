"""Fast CVT model: belt closure, clamp-force balance and a stick-slip belt clutch (docs/theory/15).

Uses the same parameters, parts, engine, tire, vehicle and scenario schedule as the segmented-belt
model, so a tune or scenario runs in either. About 1e4 times faster; peaks of shaft torque need
the full model (A15.5).
"""
from __future__ import annotations

import numpy as np

from cvtsim.core.jit import maybe_njit
from cvtsim.cvt import system
from cvtsim.cvt.geometry import _rs_closure
from cvtsim.cvt.profiles import pp_eval
from cvtsim.drivetrain.vehicle import G_ACC, driveline_torque, road_load_force, tire_force

# calibration constants (15 section 5); see cvtsim/fast/calibration.yaml for fitted values
CAL_DEFAULT = dict(K=2.0e-5, g1=0.0, g2=0.0, c_mu=1.0, s_h=0.5, c_s=0.0, c_s2=0.0)
CAL_KEYS = ("K", "g1", "g2", "c_mu", "s_h", "c_s", "c_s2")

FCH = ["t", "w_e", "w_s", "x_p", "x_s", "r_p", "r_s", "F_fw", "F_spr", "F_p", "F_s", "F_pstar",
       "T_mean", "F_cap", "F_b", "slip", "stick", "v", "w_wheel", "F_tire", "T_eng", "T_load",
       "P_slip", "T_belt", "E_kin", "W_in", "W_diss", "throttle", "grade", "tau", "s", "v_target"]
NG = 4001  # geometry table points (15 section 2)


def load_calibration():
    """Fitted constants from cvtsim/fast/calibration.yaml (15 section 5), or {} if absent."""
    import pathlib
    import yaml
    q = pathlib.Path(__file__).with_name("calibration.yaml")
    return yaml.safe_load(q.read_text())["calibration"] if q.exists() else {}


def geometry_table(g):
    """r_s, phi_p, phi_s on a uniform grid of x_p in [0, x_p_stop] (15 section 2)."""
    x = np.linspace(0.0, g["x_p_stop"], NG)
    rp = g["rp_min"] + x / (2.0 * g["tan_beta"])
    rs = np.array([_rs_closure(r, g["L"], g["C"]) for r in rp])
    php = np.pi - 2.0 * np.arcsin((rs - rp) / g["C"])
    return x, rs, php, 2.0 * np.pi - php


def build_fast(p, scenario=None, cal=None):
    """Arrays and scalars for run_fast from params p (and a scenario, 14)."""
    full = system.build(p, scenario=scenario)          # shared parameter packing
    v = dict(zip(system.PRM, full["prm"]))
    g = p["derived"]["geometry"]
    xg, rsg, phpg, phsg = geometry_table(g)
    b = p["belt"]
    c = dict(CAL_DEFAULT, **(load_calibration() if cal is None else cal))
    r_w0 = float(pp_eval(full["rk"], full["rc"], 0.0)[0])
    sc = dict(
        J_e=v["J_f"] + v["J_c"] + v["n_w"] * v["m_w"] * r_w0 ** 2, J_s=v["J_s"] + v["J_m"],
        n_w=v["n_w"], m_w=v["m_w"], k_p=v["k_p"], F_p0=v["F_p0"], g_p=v["g_p"], x_stop=g["x_p_stop"],
        k_l=v["k_l"], F_l0=v["F_l0"], k_t=v["k_t"], th_t0=v["th_t0"], rp_low=g["rp_min"], rs_max=g["rs_max"],
        tb=g["tan_beta"], sinb=np.sin(g["beta"]), mu=b["mu"], rhoA=b["rho"] * b["area"],
        mu_r=v["mu_r"], mu_b=v["mu_b"], m_veh=v["m_veh"], r_t=v["r_t"], G=v["G"], Crr=v["Crr"],
        c_a=v["c_a"], eta_d=v["eta_d"], mu_t=v["mu_t"], s0=v["s0"], v_min=v["v_min"], b_L=v["b_L"],
        h_L=v["h_L"], rpm_cut=v["rpm_cut"], gov_band=v["gov_band"], T_mot=v["T_mot"],
        C_b=v["C_b"], hA=v["hA"], f_b=v["f_b"], T_amb=v["T_amb"],
        w0=full["st"][0], m_ph=v["m_ph"], c_p=v["c_p"], **{k: c[k] for k in CAL_KEYS})
    return dict(sc=sc, xg=xg, rsg=rsg, phpg=phpg, phsg=phsg, rk=full["rk"], rc=full["rc"],
                hk=full["hk"], hc=full["hc"], eng_w=full["eng_w"], eng_T=full["eng_T"],
                sched=full["sched"], p=p, cal=c)


SC_ORDER = ("J_e", "J_s", "n_w", "m_w", "k_p", "F_p0", "g_p", "x_stop", "k_l", "F_l0", "k_t", "th_t0",
            "rp_low", "rs_max", "tb", "sinb", "mu", "rhoA", "mu_r", "mu_b", "m_veh", "r_t", "G", "Crr",
            "c_a", "eta_d", "mu_t", "s0", "v_min", "b_L", "h_L", "rpm_cut", "gov_band", "T_mot",
            "C_b", "hA", "f_b", "T_amb", "w0", "m_ph", "c_p", "K", "g1", "g2", "c_mu", "s_h", "c_s", "c_s2")


@maybe_njit
def _interp(xg, yg, x):
    n = xg.shape[0]
    f = (x - xg[0]) / (xg[n - 1] - xg[0]) * (n - 1)
    if f <= 0.0:
        return yg[0]
    if f >= n - 1:
        return yg[n - 1]
    i = int(f)
    a = f - i
    return yg[i] * (1.0 - a) + yg[i + 1] * a


@maybe_njit
def capacity(Tm, phi_min, mu_eff):
    """Belt pull at gross slip, F_cap = 2 T' tanh(mu' phi_min / 2) (15 section 3)."""
    return 2.0 * Tm * np.tanh(mu_eff * phi_min / 2.0)


@maybe_njit
def stick_force(r_p, r_s, T_e, T_L, J_e, J_s, rpd, rsd, w_e, w_s):
    """Belt pull that keeps r_p w_e = r_s w_s (15 section 4)."""
    return ((r_p * T_e / J_e + r_s * T_L / J_s + rpd * w_e - rsd * w_s)
            / (r_p * r_p / J_e + r_s * r_s / J_s))


@maybe_njit
def run_fast(sv, xg, rsg, phpg, phsg, rk, rc, hk, hc, eng_w, eng_T, sched, h, nsteps, nout, trk, drv):
    """Time scenario (sched) or, when drv[0] > 0, a track driven by distance (16):
    trk rows [s, grade, traction, v_target]; drv = [mode, gain, throttle lag, s_end]."""
    (J_e, J_s, n_w, m_w, k_p, F_p0, g_p, x_stop, k_l, F_l0, k_t, th_t0, rp_low, rs_max, tb, sinb,
     mu, rhoA, mu_r, mu_b, m_veh, r_t, G, Crr, c_a, eta_d, mu_t, s0, v_min, b_L, h_L, rpm_cut,
     gov_band, T_mot, C_b, hA, f_b, T_amb, w0, m_ph, c_p, K, g1, g2, c_mu, s_h, c_s, c_s2) = sv
    stride = max(nsteps // max(nout - 1, 1), 1)
    out = np.zeros((nsteps // stride + 1, 32))
    w_e = w0; w_s = 0.0; v = 0.0; x_p = -g_p; T_b = T_amb; a_prev = 0.0
    stick = False; F_b = 0.0; tau = 0.0; T_s_prev = 0.0; xd_p = 0.0
    W_in = 0.0; W_d = 0.0; t = 0.0; k = 0; ks = 0
    s_pos = 0.0; js = 0; thr_f = 1.0; v_tgt = 0.0; track = drv[0] > 0.0
    mup = c_mu * mu / sinb
    for n in range(nsteps + 1):
        # scenario schedule (14 section 1)
        while ks < sched.shape[0] - 2 and t >= sched[ks + 1, 0]:
            ks += 1
        if t <= sched[0, 0]:
            thr = sched[0, 1]; grd = sched[0, 2]; trac = sched[0, 3]
        elif t >= sched[-1, 0]:
            thr = sched[-1, 1]; grd = sched[-1, 2]; trac = sched[-1, 3]
        else:
            fz = (t - sched[ks, 0]) / (sched[ks + 1, 0] - sched[ks, 0])
            thr = sched[ks, 1] + fz * (sched[ks + 1, 1] - sched[ks, 1])
            grd = sched[ks, 2] + fz * (sched[ks + 1, 2] - sched[ks, 2])
            trac = sched[ks, 3] + fz * (sched[ks + 1, 3] - sched[ks, 3])
        if track:   # course and driver by distance (16 sections 1-2)
            while js < trk.shape[0] - 2 and s_pos >= trk[js + 1, 0]:
                js += 1
            fz = (s_pos - trk[js, 0]) / (trk[js + 1, 0] - trk[js, 0])
            fz = min(max(fz, 0.0), 1.0)
            grd = trk[js, 1] + fz * (trk[js + 1, 1] - trk[js, 1])
            trac = trk[js, 2] + fz * (trk[js + 1, 2] - trk[js, 2])
            v_tgt = trk[js, 3] + fz * (trk[js + 1, 3] - trk[js, 3])
            cmd = min(max(drv[1] * (v_tgt - v), 0.0), 1.0)
            thr_f += h * (cmd - thr_f) / drv[2]
            thr = thr_f
        T_e, _ = system.engine_torque(w_e, eng_w, eng_T, rpm_cut, gov_band)
        if thr != 1.0:
            T_e = thr * T_e - (1.0 - thr) * T_mot * np.tanh(w_e / 10.0)
        # geometry (15 section 2)
        engaged = x_p >= 0.0
        xq = max(x_p, 0.0)
        r_p = rp_low + xq / (2.0 * tb)
        r_s = _interp(xg, rsg, xq)
        php = _interp(xg, phpg, xq)
        phs = _interp(xg, phsg, xq)
        x_s = 2.0 * tb * (rs_max - r_s)
        # clamp forces (15 section 3)
        r_w, kap, _ = pp_eval(rk, rc, x_p)
        F_fw = n_w * m_w * r_w * kap * w_e * w_e
        F_spr = F_p0 + k_p * x_p
        F_p = F_fw - F_spr
        hx, hp, _ = pp_eval(hk, hc, x_s)
        F_s = F_l0 + k_l * x_s + k_t * (th_t0 + hx) * hp + s_h * hp * T_s_prev
        if F_s < 0.0:
            F_s = 0.0
        gfac = 1.0 + g1 * tau + g2 * tau * tau
        T_c = rhoA * (r_p * w_e) ** 2
        Tm = 0.0; F_cap = 0.0; F_ps = 0.0
        tw = tb + mu * (c_s + c_s2 * xq / x_stop)   # effective wedge with radial friction (15 section 3)
        if engaged:
            Tm = min(2.0 * tw * F_s / phs, 2.0 * tw * max(F_p, 0.0) / (gfac * php))
            F_ps = gfac * F_s * php / phs
            F_cap = capacity(Tm, min(php, phs), mup)
        # vehicle and tire (09, 14)
        w_w = w_s / G
        F_x, N_r, _ = tire_force(w_w, v, a_prev, m_veh, r_t, mu_t, s0, v_min, b_L, h_L)
        F_x *= trac
        F_road = road_load_force(v, m_veh, Crr, c_a) + m_veh * G_ACC * np.sin(grd)
        T_L = driveline_torque(F_x, w_s, r_t, G, eta_d)
        # belt clutch (15 section 4)
        rpd = xd_p / (2.0 * tb)
        rsd = (_interp(xg, rsg, min(xq + 1e-6, x_stop)) - r_s) / 1e-6 * xd_p
        vrel = r_p * w_e - r_s * w_s
        if engaged and F_cap > 0.0:
            Fst = stick_force(r_p, r_s, T_e, T_L, J_e, J_s, rpd, rsd, w_e, w_s)
            if stick and abs(Fst) <= F_cap:
                F_b = Fst
            else:
                stick = False
                sgn = 1.0 if vrel > 0.0 else (-1.0 if vrel < 0.0 else (1.0 if Fst >= 0.0 else -1.0))
                F_b = F_cap * sgn
        else:
            stick = False
            F_b = 0.0
        tau = abs(F_b) / F_cap if F_cap > 0.0 else 0.0
        T_s_prev = F_b * r_s
        # shift law (15 section 4)
        if engaged:
            F_f = mu_r * F_fw + mu_b * 2.0 * Tm
            dF = F_p - F_ps
            dF = dF - F_f if dF > F_f else (dF + F_f if dF < -F_f else 0.0)
            xd_p = K * w_e * dF / max(F_s, 1.0)
            if F_p < 0.0:            # flyweights below the spring: release into the gap
                xd_p = min(xd_p, -1e-3)
        else:   # engagement gap: movable half with its mass (05 P2', no belt)
            M_ax = m_ph + n_w * m_w * kap * kap
            Fg = F_p - c_p * xd_p
            Fg = Fg - mu_r * F_fw if Fg > mu_r * F_fw else (Fg + mu_r * F_fw if Fg < -mu_r * F_fw else 0.0)
            vg = xd_p + h * Fg / M_ax
            if x_p <= -g_p and vg < 0.0:
                vg = 0.0
        P_slip = 0.0 if stick else abs(F_b * vrel)
        P_tire = F_x * (w_w * r_t - v)
        P_drv = T_L * w_s - F_x * w_w * r_t
        E = 0.5 * J_e * w_e * w_e + 0.5 * J_s * w_s * w_s + 0.5 * m_veh * v * v
        if n % stride == 0 and k < out.shape[0]:
            o = out[k]
            o[0] = t; o[1] = w_e; o[2] = w_s; o[3] = x_p; o[4] = x_s; o[5] = r_p; o[6] = r_s
            o[7] = F_fw; o[8] = F_spr; o[9] = F_p; o[10] = F_s; o[11] = F_ps; o[12] = Tm
            o[13] = F_cap; o[14] = F_b; o[15] = vrel; o[16] = 1.0 if stick else 0.0; o[17] = v
            o[18] = w_w; o[19] = F_x; o[20] = T_e; o[21] = T_L; o[22] = P_slip; o[23] = T_b
            o[24] = E; o[25] = W_in; o[26] = W_d; o[27] = thr; o[28] = grd; o[29] = tau
            o[30] = s_pos; o[31] = v_tgt
            k += 1
        # integrate (explicit, 15 section 4)
        wd_e = (T_e - F_b * r_p) / J_e
        wd_s = (F_b * r_s - T_L) / J_s
        a_x = (F_x - F_road) / m_veh
        W_in += T_e * w_e * h
        W_d += (P_slip + P_tire + P_drv + F_road * v) * h
        w_e += h * wd_e
        w_s += h * wd_s
        v += h * a_x
        a_prev = a_x
        if not engaged:
            xd_p = vg
        x_p = min(max(x_p + h * xd_p, -g_p), x_stop)
        if x_p >= x_stop and xd_p > 0.0:
            xd_p = 0.0
        T_b += h * (f_b * P_slip - hA * (T_b - T_amb)) / C_b
        t += h
        s_pos += h * v
        if track and s_pos >= drv[3]:
            o = out[min(k, out.shape[0] - 1)]
            o[0] = t; o[1] = w_e; o[2] = w_s; o[3] = x_p; o[17] = v; o[25] = W_in; o[26] = W_d
            o[30] = s_pos; o[27] = thr; o[24] = 0.5 * J_e * w_e * w_e + 0.5 * J_s * w_s * w_s + 0.5 * m_veh * v * v
            k = min(k + 1, out.shape[0])
            break
        # re-stick when the slip speed reaches zero within capacity
        if engaged and not stick and F_cap > 0.0:
            vrel_new = (rp_low + max(x_p, 0.0) / (2.0 * tb)) * w_e - _interp(xg, rsg, max(x_p, 0.0)) * w_s
            if vrel_new * vrel <= 0.0:
                Fst = stick_force(r_p, r_s, T_e, T_L, J_e, J_s, 0.0, 0.0, w_e, w_s)
                if abs(Fst) <= F_cap:
                    stick = True
        if stick:   # hold the kinematic constraint exactly; its kinetic-energy change is booked
            r_pn = rp_low + max(x_p, 0.0) / (2.0 * tb)     # as dissipation (the belt grabbing, 15 s4)
            r_sn = _interp(xg, rsg, max(x_p, 0.0))
            # belt impulse P that restores r_p w_e = r_s w_s (loss P * v_rel / 2 >= 0)
            P = (r_pn * w_e - r_sn * w_s) / (r_pn * r_pn / J_e + r_sn * r_sn / J_s)
            E0 = 0.5 * J_e * w_e * w_e + 0.5 * J_s * w_s * w_s
            w_e -= P * r_pn / J_e
            w_s += P * r_sn / J_s
            W_d += E0 - (0.5 * J_e * w_e * w_e + 0.5 * J_s * w_s * w_s)
    return out[:k]


def simulate_fast(fm, t_end, h=2e-4, nout=None, track=None):
    """Run the fast model for t_end seconds (or until the end of `track`). Returns (n, len(FCH)).

    track: None (time scenario from build_fast) or dict(trk=(m, 4) array, gain, lag, s_end) (16).
    """
    nsteps = int(round(t_end / h))
    nout = nout or int(t_end / 1e-2) + 1
    sv = np.array([fm["sc"][k] for k in SC_ORDER], dtype=np.float64)
    if track is None:
        trk, drv = np.zeros((2, 4)), np.zeros(4)
        trk[1, 0] = 1.0
    else:
        trk = np.asarray(track["trk"], float)
        drv = np.array([1.0, track["gain"], track["lag"], track["s_end"]])
    return run_fast(sv, fm["xg"], fm["rsg"], fm["phpg"], fm["phsg"], fm["rk"], fm["rc"], fm["hk"],
                    fm["hc"], fm["eng_w"], fm["eng_T"], fm["sched"], h, nsteps, nout, trk, drv)
