"""Full CVT launch model: segmented belt + analytic cone sheaves + Lagrangian
primary/secondary with profiled ramp and helix + engine with governor + vehicle with
tire slip, driveline losses and aero + lumped belt temperature (docs/theory 03-06, 09).
Symplectic Euler, whole time loop jitted (numba) for speed."""
import numpy as np
from cvtsim.core.jit import maybe_njit
from cvtsim.belt.belt import belt_forces, belt_tensions, init_belt_on_path, flank_stiffness
from cvtsim.contact.sheave import sheave_contact
from cvtsim.cvt.profiles import pp_eval
from cvtsim.drivetrain.vehicle import G_ACC, road_load_force, tire_force, driveline_torque

# output channel indices
CH = ["t", "w_e", "w_s", "x_p", "xd_p", "x_s", "xd_s", "rstar_p", "rstar_s",
      "Fa_p", "Fa_s", "tq_p", "tq_s_fix", "tq_s_mov", "T_eng", "T_load",
      "Fx_p", "Fy_p", "Fx_s", "Fy_s", "T_tight", "T_slack", "T_max",
      "r_eff_p", "r_eff_s", "slip_p", "slip_s",
      "E_mech", "W_in", "W_diss", "Fstop_p", "Fstop_s", "Fcam",
      # 09 section 6
      "theta_e", "theta_s", "r_w", "F_fw", "F_spr_p", "F_spr_s", "F_tors", "F_cam",
      "cap_p", "cap_s", "P_slip_p", "P_slip_s", "T_belt", "v", "w_wheel", "F_tire",
      "N_rear", "P_ground", "a_x"]
IDX = {k: i for i, k in enumerate(CH)}

# packed scalar parameters (order of run_loop's unpacking)
PRM = ["m_node", "L0", "EA", "eta", "EI", "cb", "beta", "k_n", "c_n", "mu", "v_eps",
       "rp_min", "rs_max", "C", "rout_p", "rout_s",
       "J_e", "n_w", "m_w", "k_p", "F_p0", "m_ph", "c_p", "s_p",
       "J_s", "J_m", "m_m", "k_l", "F_l0", "k_t", "th_t0", "c_s", "s_s",
       "k_stop", "c_stop_p", "c_stop_s", "m_veh", "r_t", "G", "Crr", "rpm_cut", "gov_band",
       "g_p", "k_hub", "rigid", "eta_d", "c_a", "mu_t", "s0", "v_min", "b_L", "h_L",
       "C_b", "hA", "f_b", "T_amb"]
# state vector
ST = ["w_e", "x_p", "xd_p", "w_s", "x_s", "xd_s", "t", "theta_e", "theta_s", "v", "T_belt", "a_x"]


@maybe_njit
def _stop(x, xd, lo, hi, k, c):
    """Penalty end stop. Returns (force, energy, non-conservative power)."""
    pen = 0.0
    if x < lo:
        pen = x - lo
    elif x > hi:
        pen = x - hi
    if pen == 0.0:
        return 0.0, 0.0, 0.0
    f_el = -k * pen
    f = f_el - c * xd
    if (pen < 0.0 and f < 0.0) or (pen > 0.0 and f > 0.0):
        f = 0.0
    return f, 0.5 * k * pen * pen, (f_el - f) * xd


@maybe_njit
def _span_tensions(X, T, C):
    """Mean element tension in upper (y>0) and lower (y<0) free spans."""
    N = X.shape[0]
    su = 0.0; nu = 0; sl = 0.0; nl = 0
    for e in range(N):
        a = X[e]; b = X[(e + 1) % N]
        mx = 0.5 * (a[0] + b[0]); my = 0.5 * (a[1] + b[1])
        if 0.3 * C < mx < 0.7 * C:
            if my > 0:
                su += T[e]; nu += 1
            else:
                sl += T[e]; nl += 1
    return su / max(nu, 1), sl / max(nl, 1)


@maybe_njit
def _contact_stats(X, V, cx, w, N_fl):
    """Effective wrap radius (normal-force weighted) and mean tangential slip ratio."""
    sn = 0.0; sr = 0.0; ss = 0.0
    for i in range(X.shape[0]):
        if N_fl[i] > 0.0:
            dx = X[i, 0] - cx; dy = X[i, 1]
            r = np.sqrt(dx * dx + dy * dy)
            tx = -dy / r; ty = dx / r
            vt = V[i, 0] * tx + V[i, 1] * ty
            sn += N_fl[i]; sr += N_fl[i] * r
            ss += N_fl[i] * (vt - w * r)
    if sn > 0.0:
        return sr / sn, ss / sn
    return 0.0, 0.0


@maybe_njit
def engine_torque(w_e, eng_w, eng_T, rpm_cut, gov_band):
    """Net engine torque and governor torque (T_eng - tau_gov, tau_gov) at crank speed w_e.

    Chart interpolated in rpm, held at its end values outside; the governor removes
    tau_gov = T_chart * clip((rpm - (rpm_cut - band))/band, 0, 1) (09 section 3).
    """
    rpm = w_e * 30.0 / np.pi
    if rpm <= eng_w[0]:
        T_ch = eng_T[0]
    elif rpm >= eng_w[-1]:
        T_ch = eng_T[-1]
    else:
        T_ch = np.interp(rpm, eng_w, eng_T)
    T_e = T_ch * min(1.0, max(0.0, (rpm_cut - rpm) / gov_band))
    return T_e, T_ch - T_e


@maybe_njit
def primary_accel(u, ud, w, Q_w, Q_u, J_e, n_w, m_w, m_ph, rk, rc):
    """(P1'),(P2') with ramp profile (rk, rc). Returns (w_dot, u_ddot, r_w, kappa, J, M).

    Q_w = T_eng + tau_belt - tau_gov; Q_u = -F_spr - c_p u_dot - F_a,p + F_stop.
    """
    r_w, kap, dkap = pp_eval(rk, rc, u)
    J = J_e + n_w * m_w * r_w * r_w
    M = m_ph + n_w * m_w * kap * kap
    wd = (Q_w - 2 * n_w * m_w * r_w * kap * ud * w) / J
    udd = (n_w * m_w * kap * r_w * w * w + Q_u - n_w * m_w * kap * dkap * ud * ud) / M
    return wd, udd, r_w, kap, J, M


@maybe_njit
def secondary_accel(xd, Q_th, Q_x, J_s, J_m, m_m, hp, hpp):
    """Secondary 2x2 system (09 section 2). Returns (w_dot, x_ddot, a11, a12, a22).

    Q_th = tau_fix + tau_mov - T_load; Q_x = -h' tau_mov + F_a,s - F_spring - k_t(th0+h)h'
    - c_s x_dot + F_stop (all generalized forces except the h'' velocity terms).
    """
    a11 = J_s + J_m; a12 = -J_m * hp; a22 = m_m + J_m * hp * hp
    b1 = Q_th + J_m * hpp * xd * xd
    b2 = Q_x - J_m * hp * hpp * xd * xd
    det = a11 * a22 - a12 * a12
    return (a22 * b1 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det, a11, a12, a22


@maybe_njit
def belt_temperature_step(T_b, P, C_b, hA, T_amb, h):
    """Advance C_b dT/dt = P - hA (T - T_amb) by h, exact for constant P (09 section 4)."""
    T_inf = T_amb + P / hA
    return T_inf + (T_b - T_inf) * np.exp(-h * hA / C_b)


@maybe_njit
def run_loop(X, V, st, prm, eng_w, eng_T, rk, rc, hk, hc, h, nsteps, nout, nframes, lock_engine):
    """Integrate nsteps of h. st: state (see ST), prm: scalars (see PRM), (rk, rc) ramp and
    (hk, hc) helix piecewise polynomials. Returns (out, frame_t, frame_X, frame_T, frame_s)."""
    (m_node, L0, EA, eta, EI, cb, beta, k_n, c_n, mu, v_eps,
     rp_min, rs_max, C, rout_p, rout_s,
     J_e, n_w, m_w, k_p, F_p0, m_ph, c_p, s_p,
     J_s, J_m, m_m, k_l, F_l0, k_t, th_t0, c_s, s_s,
     k_stop, c_stop_p, c_stop_s, m_veh, r_t, G, Crr, rpm_cut, gov_band,
     g_p, k_hub, rigid, eta_d, c_a, mu_t, s0, v_min, b_L, h_L,
     C_b, hA, f_b, T_amb) = prm[:56]
    tb = np.tan(beta)
    w_e, x_p, xd_p, w_s, x_s, xd_s, t = st[0], st[1], st[2], st[3], st[4], st[5], st[6]
    th_e, th_s, v, T_b, a_prev = st[7], st[8], st[9], st[10], st[11]
    tire = rigid == 0.0
    stride = max(nsteps // nout, 1)
    out = np.zeros((nsteps // stride + 1, 52))
    fstride = max(nsteps // max(nframes - 1, 1), 1)
    nf = min(nsteps // fstride + 1, nframes) if nframes > 0 else 0
    N = X.shape[0]
    f_t = np.zeros(nf); f_X = np.zeros((nf, N, 2)); f_T = np.zeros((nf, N)); f_s = np.zeros((nf, 5))
    W_in = 0.0; W_diss = 0.0
    k = 0; kf = 0
    for n in range(nsteps + 1):
        rsp = rp_min + x_p / (2 * tb); rsp_d = xd_p / (2 * tb)
        rss = rs_max - x_s / (2 * tb); rss_d = -xd_s / (2 * tb)
        hx, hp, hpp = pp_eval(hk, hc, x_s)
        fb, Ub, Pb = belt_forces(X, V, L0, EA, eta, EI, cb, True)
        w_mov_s = w_s - hp * xd_s
        fp, Fap, tqpf, tqpm, Pfp, Upp, dlp, Nfp = sheave_contact(
            X, V, 0.0, 0.0, rsp, rsp_d, beta, k_n, c_n, mu, v_eps, rout_p, w_e, w_e)
        fs, Fas, tqsf, tqsm, Pfs, Ups, dls, Nfs = sheave_contact(
            X, V, C, 0.0, rss, rss_d, beta, k_n, c_n, mu, v_eps, rout_s, w_s, w_mov_s)
        # free-running hub bearing under the open primary: frictionless radial support at rp_min
        U_hub = 0.0
        for i in range(N):
            r = np.sqrt(X[i, 0] ** 2 + X[i, 1] ** 2)
            if r < rp_min:
                pen = rp_min - r
                fp[i, 0] += k_hub * pen * X[i, 0] / r
                fp[i, 1] += k_hub * pen * X[i, 1] / r
                U_hub += 0.5 * k_hub * pen * pen
        # normal non-conservative power: sum 2 (N - k delta) delta_dot
        Pn = 0.0
        for i in range(N):
            for c in range(2):
                if c == 0:
                    dx = X[i, 0]; dy = X[i, 1]; rs_ = rsp; rsd = rsp_d; Ni = Nfp[i]; di = dlp[i]
                else:
                    dx = X[i, 0] - C; dy = X[i, 1]; rs_ = rss; rsd = rss_d; Ni = Nfs[i]; di = dls[i]
                if di > 0.0:
                    r = np.sqrt(dx * dx + dy * dy)
                    rdot = (dx * V[i, 0] + dy * V[i, 1]) / r
                    dd = (rsd - rdot) * np.sin(beta)
                    Pn += 2.0 * (Ni - k_n * di) * dd
        T_e, T_gov = engine_torque(w_e, eng_w, eng_T, rpm_cut, gov_band)
        # primary (P1'),(P2'); travel from fully open u = x_p + g_p
        u = x_p + g_p
        Fst_p, Ust_p, Pst_p = _stop(x_p, xd_p, -g_p, s_p, k_stop, c_stop_p)
        F_spr_p = F_p0 + k_p * u
        wd_e, xdd_p, r_w, kap, Jp, Mx_p = primary_accel(
            u, xd_p, w_e, T_e + tqpf + tqpm, -F_spr_p - c_p * xd_p - Fap + Fst_p,
            J_e, n_w, m_w, m_ph, rk, rc)
        if lock_engine:
            wd_e = 0.0
        # vehicle and driveline
        w_w = w_s / G
        if tire:
            F_x, N_r, _ = tire_force(w_w, v, a_prev, m_veh, r_t, mu_t, s0, v_min, b_L, h_L)
            F_road = road_load_force(v, m_veh, Crr, c_a)
            T_load = driveline_torque(F_x, w_s, r_t, G, eta_d)
        else:
            v = w_s * r_t / G
            F_road = road_load_force(v, m_veh, Crr, c_a)
            T_load = driveline_torque(F_road, w_s, r_t, G, eta_d)
        # secondary
        Fst_s, Ust_s, Pst_s = _stop(x_s, xd_s, 0.0, s_s, k_stop, c_stop_s)
        Fcam = -hp * tqsm
        F_spr_s = F_l0 + k_l * x_s
        F_tors = k_t * (th_t0 + hx) * hp
        Qt = tqsf + tqsm - T_load
        Qx = Fcam + Fas - F_spr_s - F_tors - c_s * xd_s + Fst_s
        wd_s, xdd_s, a11, a12, a22 = secondary_accel(xd_s, Qt, Qx, J_s, J_m, m_m, hp, hpp)
        if tire:
            a_x = (F_x - F_road) / m_veh
            P_veh = T_load * w_s - F_x * v + F_road * v
            KEv = 0.5 * m_veh * v * v
        else:
            a_x = wd_s * r_t / G
            F_x = F_road + m_veh * a_x
            N_r = m_veh * G_ACC * b_L + m_veh * a_prev * h_L
            P_veh = T_load * w_s
            KEv = 0.0
        # energy bookkeeping (state at t_n)
        KEb = 0.5 * m_node * np.sum(V * V)
        KEp = 0.5 * Jp * w_e * w_e + 0.5 * Mx_p * xd_p * xd_p
        KEs = 0.5 * (a11 * w_s * w_s + 2 * a12 * w_s * xd_s + a22 * xd_s * xd_s)
        Usp = F_p0 * u + 0.5 * k_p * u ** 2
        Uss = F_l0 * x_s + 0.5 * k_l * x_s * x_s + 0.5 * k_t * (th_t0 + hx) ** 2
        E = KEb + KEp + KEs + KEv + Ub + Upp + Ups + Usp + Uss + Ust_p + Ust_s + U_hub
        Pst = Pst_p + Pst_s
        P_d = Pb + Pfp + Pfs + Pn + c_p * xd_p ** 2 + c_s * xd_s ** 2 + P_veh + Pst
        P_in = T_e * w_e if not lock_engine else 0.0
        if n % stride == 0 and k < out.shape[0]:
            T = belt_tensions(X, V, L0, EA, eta, True)
            Tt, Tsl = _span_tensions(X, T, C)
            rep, slp = _contact_stats(X, V, 0.0, w_e, Nfp)
            res, sls = _contact_stats(X - np.array([C, 0.0]), V, 0.0, w_s, Nfs)
            Fpx = 0.0; Fpy = 0.0; Fsx = 0.0; Fsy = 0.0
            for i in range(N):
                Fpx -= fp[i, 0]; Fpy -= fp[i, 1]; Fsx -= fs[i, 0]; Fsy -= fs[i, 1]
            o = out[k]
            o[0] = t; o[1] = w_e; o[2] = w_s; o[3] = x_p; o[4] = xd_p; o[5] = x_s; o[6] = xd_s
            o[7] = rsp; o[8] = rss; o[9] = Fap; o[10] = Fas; o[11] = tqpf + tqpm
            o[12] = tqsf; o[13] = tqsm; o[14] = T_e; o[15] = T_load
            o[16] = Fpx; o[17] = Fpy; o[18] = Fsx; o[19] = Fsy
            o[20] = Tt; o[21] = Tsl; o[22] = T.max()
            o[23] = rep; o[24] = res; o[25] = slp; o[26] = sls
            o[27] = E; o[28] = W_in; o[29] = W_diss; o[30] = Fst_p; o[31] = Fst_s; o[32] = Fcam
            o[33] = th_e; o[34] = th_s; o[35] = r_w; o[36] = n_w * m_w * r_w * kap * w_e * w_e
            o[37] = F_spr_p; o[38] = F_spr_s; o[39] = F_tors; o[40] = Fcam
            o[41] = 2 * mu * rep * Fap / np.cos(beta); o[42] = 2 * mu * res * Fas / np.cos(beta)
            o[43] = Pfp; o[44] = Pfs; o[45] = T_b; o[46] = v; o[47] = w_w; o[48] = F_x
            o[49] = N_r; o[50] = F_x * v; o[51] = a_x
            k += 1
        if nf > 0 and n % fstride == 0 and kf < nf:
            f_t[kf] = t
            f_X[kf] = X
            f_T[kf] = belt_tensions(X, V, L0, EA, eta, True)
            f_s[kf, 0] = x_p; f_s[kf, 1] = x_s; f_s[kf, 2] = th_e; f_s[kf, 3] = th_s; f_s[kf, 4] = r_w
            kf += 1
        if n == nsteps:
            break
        W_in += P_in * h
        W_diss += P_d * h
        T_b = belt_temperature_step(T_b, f_b * (Pfp + Pfs) + Pb, C_b, hA, T_amb, h)
        # symplectic Euler
        for i in range(N):
            V[i, 0] += h * (fb[i, 0] + fp[i, 0] + fs[i, 0]) / m_node
            V[i, 1] += h * (fb[i, 1] + fp[i, 1] + fs[i, 1]) / m_node
            X[i, 0] += h * V[i, 0]
            X[i, 1] += h * V[i, 1]
        w_e += h * wd_e
        th_e += h * w_e
        xd_p += h * xdd_p; x_p += h * xd_p
        w_s += h * wd_s
        th_s += h * w_s
        xd_s += h * xdd_s; x_s += h * xd_s
        if tire:
            v += h * a_x
        a_prev = a_x
        t += h
    st[0] = w_e; st[1] = x_p; st[2] = xd_p; st[3] = w_s; st[4] = x_s; st[5] = xd_s; st[6] = t
    st[7] = th_e; st[8] = th_s; st[9] = v; st[10] = T_b; st[11] = a_prev
    return out[:k], f_t[:kf], f_X[:kf], f_T[:kf], f_s[:kf]


def build(p, dt=None, kn_scale=1.0):
    """Assemble the launch model from a params dict (cvtsim.io.params.load())."""
    g = p["derived"]["geometry"]
    d = p["derived"]
    beta = d["beta"]
    b, pr, se, eng, veh = p["belt"], p["primary"], p["secondary"], p["engine"], p["vehicle"]
    ramp, helix = d["ramp"], d["helix"]
    N = int(p["sim"]["belt_nodes"])
    X, L0 = init_belt_on_path(g["rp_min"], g["rs_low"], g["C"], N)
    X = np.ascontiguousarray(X, dtype=np.float64)
    # stress-free rest length = mean initial chord
    dX = np.roll(X, -1, axis=0) - X
    L0 = float(np.mean(np.hypot(dX[:, 0], dX[:, 1])))
    m_node = b["rho"] * b["area"] * L0
    k_n = kn_scale * flank_stiffness(b["G_hot"], b["flank_height"], L0, b["top_width"])
    c_n = 2.0 * b["normal_damping_ratio"] * np.sqrt(k_n * m_node)
    margin = p["sheave"]["r_outer_margin"]
    m_ax_p = pr["m_axial"] + pr["n_weights"] * d["m_w"] * ramp(0.0)[1] ** 2
    c_stop_p = 2 * 0.05 * np.sqrt(p["sheave"]["stop_k"] * m_ax_p)
    c_stop_s = 2 * 0.05 * np.sqrt(p["sheave"]["stop_k"] * se["m_axial"])
    G = veh["final_drive"]
    rigid = veh["model"] == "rigid"
    if veh["model"] not in ("rigid", "tire"):
        raise ValueError(f"vehicle.model must be 'rigid' or 'tire', not {veh['model']!r}")
    J_s = se["J_shaft"] + (veh["mass"] * veh["tire_radius"] ** 2 / G ** 2 if rigid else 0.0) + veh["J_wheels"] / G ** 2
    L_wb = veh["wheelbase"]
    th = d["belt_thermal"]
    vals = dict(
        m_node=m_node, L0=L0, EA=b["EA"], eta=b["eta"], EI=b["EI"], cb=0.0, beta=beta, k_n=k_n, c_n=c_n,
        mu=b["mu"], v_eps=b["v_eps"], rp_min=g["rp_min"], rs_max=g["rs_max"], C=g["C"],
        rout_p=p["geometry"]["rp_mech_max"] + margin, rout_s=g["rs_max"] + margin,
        J_e=eng["J_e"], n_w=pr["n_weights"], m_w=d["m_w"], k_p=pr["k_spring"], F_p0=pr["F_preload"],
        m_ph=pr["m_axial"], c_p=pr["c_axial"], s_p=g["x_p_stop"],
        J_s=J_s, J_m=se["J_movable"], m_m=se["m_axial"], k_l=se["k_linear"], F_l0=d["F_l0"],
        k_t=se["k_torsion"], th_t0=d["theta_t0"], c_s=se["c_axial"], s_s=se["stroke"],
        k_stop=p["sheave"]["stop_k"], c_stop_p=c_stop_p, c_stop_s=c_stop_s,
        m_veh=veh["mass"], r_t=veh["tire_radius"], G=G, Crr=veh["Crr"],
        rpm_cut=eng["rpm_cut"], gov_band=eng["governor_band"], g_p=pr["engagement_gap"], k_hub=k_n,
        rigid=float(rigid), eta_d=veh["eta"], c_a=0.5 * veh["rho_air"] * veh["CdA"], mu_t=veh["mu_t"],
        s0=veh["s0"], v_min=veh["v_min"], b_L=veh["b_over_L"], h_L=veh["h_cg"] / L_wb,
        C_b=th["C_b"], hA=th["hA"], f_b=b["f_slip_heat"], T_amb=th["T_amb"])
    prm = np.array([vals[k] for k in PRM], dtype=np.float64)
    st = np.zeros(len(ST))
    st[0] = eng["idle_rpm"] * np.pi / 30
    st[1] = -pr["engagement_gap"]
    st[4] = g["x_s0"]
    st[10] = th["T_amb"]
    return dict(X=X, V=np.zeros_like(X), st=st, prm=prm,
                eng_w=np.asarray(eng["rpm"], float), eng_T=np.asarray(eng["torque"], float),
                rk=ramp.xk, rc=ramp.c, hk=helix.xk, hc=helix.c,
                h=float(dt or p["sim"]["dt"]), p=p)


def _loop(model, st, prm, nsteps, nout, nframes, lock):
    return run_loop(model["X"], model["V"], st, prm, model["eng_w"], model["eng_T"],
                    model["rk"], model["rc"], model["hk"], model["hc"], model["h"],
                    nsteps, nout, nframes, lock)


def relax(model, t_relax=0.02):
    """Settle the belt in the stationary sheaves (engine locked, nothing turning)."""
    st = model["st"].copy()
    w0 = st[0]
    st[0] = 0.0
    prm = model["prm"].copy()
    prm[PRM.index("eta")] = max(prm[PRM.index("eta")], 0.0)
    n = int(t_relax / model["h"])
    for _ in range(4):  # damp by zeroing velocity between chunks (kinetic relaxation)
        _loop(model, st, prm, n // 4, 2, 0, True)
        model["V"][:] = 0.0
    for key in ("xd_p", "xd_s", "w_s", "t", "theta_e", "theta_s", "v", "a_x"):
        st[ST.index(key)] = 0.0
    st[ST.index("T_belt")] = model["st"][ST.index("T_belt")]
    st[1] = min(st[1], -model["p"]["primary"]["engagement_gap"])
    st[0] = w0
    model["st"] = st


def simulate(model, t_end, nout=3000, n_frames=400):
    """Run t_end seconds. Returns (res, snaps) per docs/viewer_contract.md.

    res: (n_out, len(CH)); snaps: t (F,), belt (F, N, 2) node positions, tension (F, N)
    element tensions (element k joins nodes k and k+1), and x_p, x_s, theta_e, theta_s, r_w (F,).
    """
    nsteps = int(round(t_end / model["h"]))
    res, ft, fX, fT, fs = _loop(model, model["st"], model["prm"], nsteps, nout, n_frames, False)
    snaps = dict(t=ft, belt=fX, tension=fT)
    for j, key in enumerate(("x_p", "x_s", "theta_e", "theta_s", "r_w")):
        snaps[key] = fs[:, j]
    return res, snaps
