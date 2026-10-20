"""Full CVT launch model: segmented belt + analytic cone sheaves + Lagrangian
primary/secondary + engine + vehicle (docs/theory 03-06). Symplectic Euler,
whole time loop jitted (numba) for speed."""
import numpy as np
from cvtsim.core.jit import maybe_njit
from cvtsim.belt.belt import belt_forces, belt_tensions, init_belt_on_path, flank_stiffness
from cvtsim.contact.sheave import sheave_contact

# output channel indices
CH = ["t", "w_e", "w_s", "x_p", "xd_p", "x_s", "xd_s", "rstar_p", "rstar_s",
      "Fa_p", "Fa_s", "tq_p", "tq_s_fix", "tq_s_mov", "T_eng", "T_load",
      "Fx_p", "Fy_p", "Fx_s", "Fy_s", "T_tight", "T_slack", "T_max",
      "r_eff_p", "r_eff_s", "slip_p", "slip_s",
      "E_mech", "W_in", "W_diss", "Fstop_p", "Fstop_s", "Fcam"]
IDX = {k: i for i, k in enumerate(CH)}


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
def run_loop(X, V, st, prm, eng_w, eng_T, h, nsteps, nout, lock_engine):
    """st = [w_e, x_p, xd_p, w_s, x_s, xd_s, t]; prm: float array (see pack_params)."""
    (m_node, L0, EA, eta, EI, cb, beta, k_n, c_n, mu, v_eps,
     rp_min, rs_max, C, rout_p, rout_s,
     J_e, n_w, m_w, r_w0, kappa, k_p, F_p0, m_ph, c_p, s_p,
     J_s, J_m, m_m, p_h, sigma, k_l, F_l0, k_t, th_t0, c_s, s_s,
     k_stop, c_stop_p, c_stop_s, m_veh, r_t, G, Crr, v0, rpm_cut, g_p, k_hub) = prm[:48]
    tb = np.tan(beta)
    w_e, x_p, xd_p, w_s, x_s, xd_s, t = st[0], st[1], st[2], st[3], st[4], st[5], st[6]
    stride = max(nsteps // nout, 1)
    out = np.zeros((nsteps // stride + 1, 33))
    W_in = 0.0; W_diss = 0.0
    k = 0
    for n in range(nsteps + 1):
        rsp = rp_min + x_p / (2 * tb); rsp_d = xd_p / (2 * tb)
        rss = rs_max - x_s / (2 * tb); rss_d = -xd_s / (2 * tb)
        fb, Ub, Pb = belt_forces(X, V, L0, EA, eta, EI, cb, True)
        w_mov_s = w_s + sigma * xd_s / p_h
        fp, Fap, tqpf, tqpm, Pfp, Upp, dlp, Nfp = sheave_contact(
            X, V, 0.0, 0.0, rsp, rsp_d, beta, k_n, c_n, mu, v_eps, rout_p, w_e, w_e)
        fs, Fas, tqsf, tqsm, Pfs, Ups, dls, Nfs = sheave_contact(
            X, V, C, 0.0, rss, rss_d, beta, k_n, c_n, mu, v_eps, rout_s, w_s, w_mov_s)
        # free-running hub bearing under the open primary: frictionless radial support at rp_min
        U_hub = 0.0
        for i in range(X.shape[0]):
            r = np.sqrt(X[i, 0] ** 2 + X[i, 1] ** 2)
            if r < rp_min:
                pen = rp_min - r
                fp[i, 0] += k_hub * pen * X[i, 0] / r
                fp[i, 1] += k_hub * pen * X[i, 1] / r
                U_hub += 0.5 * k_hub * pen * pen
        # normal non-conservative power: sum 2 (N - k delta) delta_dot
        Pn = 0.0
        for i in range(X.shape[0]):
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
        # engine torque
        rpm = w_e * 30.0 / np.pi
        if rpm <= eng_w[0]:
            T_e = eng_T[0]
        elif rpm >= eng_w[-1]:
            T_e = eng_T[-1] * max(0.0, (rpm_cut - rpm) / (rpm_cut - eng_w[-1]))
        else:
            T_e = np.interp(rpm, eng_w, eng_T)
        # primary
        # travel from fully open: x_p + g_p (spring and flyweights referenced there)
        r_w = r_w0 + kappa * (x_p + g_p)
        Jp = J_e + n_w * m_w * r_w * r_w
        Mx_p = m_ph + n_w * m_w * kappa * kappa
        Fst_p, Ust_p, Pst_p = _stop(x_p, xd_p, -g_p, s_p, k_stop, c_stop_p)
        Fcf = n_w * m_w * kappa * r_w * w_e * w_e
        wd_e = (T_e + tqpf + tqpm - 2 * n_w * m_w * r_w * kappa * xd_p * w_e) / Jp
        xdd_p = (Fcf - (F_p0 + k_p * (x_p + g_p)) - c_p * xd_p - Fap + Fst_p) / Mx_p
        if lock_engine:
            wd_e = 0.0
        # secondary
        v = w_s * r_t / G
        T_load = Crr * m_veh * 9.81 * np.tanh(v / v0) * r_t / G
        Fst_s, Ust_s, Pst_s = _stop(x_s, xd_s, 0.0, s_s, k_stop, c_stop_s)
        Fcam = sigma * tqsm / p_h
        Qt = tqsf + tqsm - T_load
        Qx = Fcam + Fas - (F_l0 + k_l * x_s) - k_t * (th_t0 + x_s / p_h) / p_h - c_s * xd_s + Fst_s
        a11 = J_s + J_m; a12 = sigma * J_m / p_h; a22 = m_m + J_m / (p_h * p_h)
        det = a11 * a22 - a12 * a12
        wd_s = (a22 * Qt - a12 * Qx) / det
        xdd_s = (a11 * Qx - a12 * Qt) / det
        # energy bookkeeping (state at t_n)
        KEb = 0.5 * m_node * np.sum(V * V)
        KEp = 0.5 * Jp * w_e * w_e + 0.5 * Mx_p * xd_p * xd_p
        KEs = 0.5 * (a11 * w_s * w_s + 2 * a12 * w_s * xd_s + a22 * xd_s * xd_s)
        Usp = F_p0 * (x_p + g_p) + 0.5 * k_p * (x_p + g_p) ** 2
        Uss = F_l0 * x_s + 0.5 * k_l * x_s * x_s + 0.5 * k_t * (th_t0 + x_s / p_h) ** 2
        E = KEb + KEp + KEs + Ub + Upp + Ups + Usp + Uss + Ust_p + Ust_s + U_hub
        Pst = Pst_p + Pst_s
        P_d = Pb + Pfp + Pfs + Pn + c_p * xd_p ** 2 + c_s * xd_s ** 2 + T_load * w_s + Pst
        P_in = T_e * w_e if not lock_engine else 0.0
        if n % stride == 0 and k < out.shape[0]:
            T = belt_tensions(X, V, L0, EA, eta, True)
            Tt, Tsl = _span_tensions(X, T, C)
            rep, slp = _contact_stats(X, V, 0.0, w_e, Nfp)
            res, sls = _contact_stats(X - np.array([C, 0.0]), V, 0.0, w_s, Nfs)
            Fpx = 0.0; Fpy = 0.0; Fsx = 0.0; Fsy = 0.0
            for i in range(X.shape[0]):
                Fpx -= fp[i, 0]; Fpy -= fp[i, 1]; Fsx -= fs[i, 0]; Fsy -= fs[i, 1]
            o = out[k]
            o[0] = t; o[1] = w_e; o[2] = w_s; o[3] = x_p; o[4] = xd_p; o[5] = x_s; o[6] = xd_s
            o[7] = rsp; o[8] = rss; o[9] = Fap; o[10] = Fas; o[11] = tqpf + tqpm
            o[12] = tqsf; o[13] = tqsm; o[14] = T_e; o[15] = T_load
            o[16] = Fpx; o[17] = Fpy; o[18] = Fsx; o[19] = Fsy
            o[20] = Tt; o[21] = Tsl; o[22] = T.max()
            o[23] = rep; o[24] = res; o[25] = slp; o[26] = sls
            o[27] = E; o[28] = W_in; o[29] = W_diss; o[30] = Fst_p; o[31] = Fst_s; o[32] = Fcam
            k += 1
        if n == nsteps:
            break
        W_in += P_in * h
        W_diss += P_d * h
        # symplectic Euler
        for i in range(X.shape[0]):
            V[i, 0] += h * (fb[i, 0] + fp[i, 0] + fs[i, 0]) / m_node
            V[i, 1] += h * (fb[i, 1] + fp[i, 1] + fs[i, 1]) / m_node
            X[i, 0] += h * V[i, 0]
            X[i, 1] += h * V[i, 1]
        w_e += h * wd_e
        xd_p += h * xdd_p; x_p += h * xd_p
        w_s += h * wd_s
        xd_s += h * xdd_s; x_s += h * xd_s
        t += h
    st[0] = w_e; st[1] = x_p; st[2] = xd_p; st[3] = w_s; st[4] = x_s; st[5] = xd_s; st[6] = t
    return out[:k]


def build(p, dt=None, kn_scale=1.0):
    """Assemble the launch model from a params dict (cvtsim.io.params.load())."""
    g = p["derived"]["geometry"]
    beta = p["derived"]["beta"]
    b, pr, se, eng, veh = p["belt"], p["primary"], p["secondary"], p["engine"], p["vehicle"]
    N = int(p["sim"]["belt_nodes"])
    X, L0 = init_belt_on_path(g["rp_min"], g["rs_max"], g["C"], N)
    X = np.ascontiguousarray(X, dtype=np.float64)
    # stress-free rest length = mean initial chord
    d = np.roll(X, -1, axis=0) - X
    L0 = float(np.mean(np.hypot(d[:, 0], d[:, 1])))
    m_node = b["rho"] * b["area"] * L0
    k_n = kn_scale * flank_stiffness(b["G_hot"], b["flank_height"], L0, b["top_width"])
    c_n = 2.0 * b["normal_damping_ratio"] * np.sqrt(k_n * m_node)
    margin = p["sheave"]["r_outer_margin"]
    m_ax_p = pr["m_axial"] + pr["n_weights"] * p["derived"]["m_w"] * pr["kappa"] ** 2
    c_stop_p = 2 * 0.05 * np.sqrt(p["sheave"]["stop_k"] * m_ax_p)
    c_stop_s = 2 * 0.05 * np.sqrt(p["sheave"]["stop_k"] * se["m_axial"])
    J_s = se["J_shaft"] + veh["mass"] * veh["tire_radius"] ** 2 / veh["final_drive"] ** 2
    prm = np.array([
        m_node, L0, b["EA"], b["eta"], b["EI"], 0.0, beta, k_n, c_n, b["mu"], b["v_eps"],
        g["rp_min"], g["rs_max"], g["C"], g["rp_max"] + margin, g["rs_max"] + margin,
        eng["J_e"], pr["n_weights"], p["derived"]["m_w"], pr["r_w0"], pr["kappa"],
        pr["k_spring"], pr["F_preload"], pr["m_axial"], pr["c_axial"], pr["stroke"],
        J_s, se["J_movable"], se["m_axial"], p["derived"]["p_h"], -1.0,
        se["k_linear"], p["derived"]["F_l0"], se["k_torsion"], p["derived"]["theta_t0"],
        se["c_axial"], se["stroke"],
        p["sheave"]["stop_k"], c_stop_p, c_stop_s,
        veh["mass"], veh["tire_radius"], veh["final_drive"], veh["Crr"], 0.05, eng["rpm_cut"],
        pr["engagement_gap"], k_n],
        dtype=np.float64)
    st = np.array([eng["idle_rpm"] * np.pi / 30, -pr["engagement_gap"], 0.0, 0.0, 0.0, 0.0, 0.0])
    return dict(X=X, V=np.zeros_like(X), st=st, prm=prm,
                eng_w=np.asarray(eng["rpm"], float), eng_T=np.asarray(eng["torque"], float),
                h=float(dt or p["sim"]["dt"]), p=p)


def relax(model, t_relax=0.02):
    """Settle the belt in the stationary sheaves (engine locked, nothing turning)."""
    st = model["st"].copy()
    w0 = st[0]
    st[0] = 0.0
    prm = model["prm"].copy()
    prm[3] = max(prm[3], 0.0)
    X, V = model["X"], model["V"]
    n = int(t_relax / model["h"])
    for _ in range(4):  # damp by zeroing velocity between chunks (kinetic relaxation)
        run_loop(X, V, st, prm, model["eng_w"], model["eng_T"], model["h"], n // 4, 2, True)
        V[:] = 0.0
    st[2] = 0.0; st[5] = 0.0; st[3] = 0.0; st[6] = 0.0
    st[1] = min(st[1], -model["p"]["primary"]["engagement_gap"])
    st[0] = w0
    model["st"] = st


def simulate(model, t_end, nout=3000):
    nsteps = int(round(t_end / model["h"]))
    return run_loop(model["X"], model["V"], model["st"], model["prm"],
                    model["eng_w"], model["eng_T"], model["h"], nsteps, nout, False)
