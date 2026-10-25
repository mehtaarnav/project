"""Shaft/bearing load and stress histories from a launch result (docs/theory/07 A7.1:
quasi-static shaft response, linear superposition of unit-load solutions)."""
import numpy as np
from cvtsim.fe.shaft import Shaft
from cvtsim.fe.peterson import kt_shoulder_bending, kt_shoulder_torsion

K_TORSION_GROUND = 1e9  # N*m/rad, torque reaction station (rigid in practice)


def _build(sd, mat, extra_z):
    z = np.asarray(sd["stations"], float)
    zs = np.unique(np.concatenate([z, np.asarray(extra_z, float)]))
    # map diameters (per original element) onto the refined stations
    dia = []
    for a, b in zip(zs[:-1], zs[1:]):
        k = np.searchsorted(z, 0.5 * (a + b)) - 1
        dia.append(sd["diameters"][min(max(k, 0), len(sd["diameters"]) - 1)])
    return Shaft(zs, np.asarray(dia), mat["E"], mat["nu"])


def _unit_responses(shaft, load_z, crit_z, bearings, k_t_z):
    """Responses to unit Fx, Fy, Tz at each load station: bearing radial forces and section loads."""
    for zb in bearings:
        shaft.add_bearing(zb, k_r=1e12, k_a=0.0, k_t=K_TORSION_GROUND if np.isclose(zb, k_t_z) else 0.0)
    if not any(np.isclose(zb, k_t_z) for zb in bearings):
        shaft.add_bearing(k_t_z, k_r=0.0, k_a=0.0, k_t=K_TORSION_GROUND)
    out = {}
    for zl in load_z:
        for comp, vec in (("Fx", (1, 0, 0, 0, 0, 0)), ("Fy", (0, 1, 0, 0, 0, 0)), ("Tz", (0, 0, 0, 0, 0, 1))):
            r = shaft.solve([(zl, *vec)])
            reac = np.array([[rb[0], rb[1]] for rb in r.reactions[:len(bearings)]])  # (nb, 2)
            Vx, Vy, Mx, My, T, Nn = r.internal(crit_z)
            out[(zl, comp)] = (reac, np.array([Mx, My, T]))
    return out


def _history(res, CH, sd, mat, loads, bearings, k_t_z, crit):
    """loads: list of (z, Fx_series, Fy_series, Tz_series)."""
    shaft = _build(sd, mat, [l[0] for l in loads] + list(bearings) + [k_t_z, crit["z"]])
    U = _unit_responses(shaft, [l[0] for l in loads], crit["z"], bearings, k_t_z)
    n = res.shape[0]
    reac = np.zeros((n, len(bearings), 2))
    sec = np.zeros((n, 3))
    for zl, fx, fy, tz in loads:
        for comp, ser in (("Fx", fx), ("Fy", fy), ("Tz", tz)):
            R, S = U[(zl, comp)]
            reac += ser[:, None, None] * R[None]
            sec += ser[:, None] * S[None]
    M = np.hypot(sec[:, 0], sec[:, 1])
    T = sec[:, 2]
    d = crit["d"]
    ktb = kt_shoulder_bending(crit["D"], d, crit["r"])
    kts = kt_shoulder_torsion(crit["D"], d, crit["r"])
    sb = 32 * M / (np.pi * d ** 3)
    tau = 16 * T / (np.pi * d ** 3)
    vm = np.sqrt((ktb * sb) ** 2 + 3 * (kts * tau) ** 2)
    return np.linalg.norm(reac, axis=2), M, T, vm, (ktb, kts)


def shaft_histories(res, p, CH):
    c = lambda k: res[:, CH.index(k)]
    mat = p["shafts"]["material"]
    sp, ss = p["shafts"]["primary"], p["shafts"]["secondary"]
    zero = np.zeros(res.shape[0])
    # primary: belt hub force at the sheave; torsion is the PTO torque T_pto (11 section 2), which
    # carries the belt torque plus clutch inertia and the firing pulses. Applied at the sheave with
    # the sign of a load on the shaft (-T_pto, as tq_p ~ -T_pto in steady drive), reacted at the crank.
    # Crank-train gas and inertia forces on the main bearings are not included (belt load only).
    Rp, Mp, Tp, vmp, ktp = _history(res, CH, sp, mat,
                                    [(sp["sheave_z"], c("Fx_p"), c("Fy_p"), -c("T_pto"))],
                                    sp["bearings"], sp["torque_in_z"], sp["critical"])
    # secondary: belt hub force + belt torque at sheave; output gear force (PLACEHOLDER: along -Y,
    # magnitude T/r_pitch) and torque reaction at output station
    Ts = c("tq_s_fix") + c("tq_s_mov")
    Fout = Ts / ss["output_pitch_radius"]
    Rs, Ms, Tsec, vms, kts = _history(res, CH, ss, mat,
                                      [(ss["sheave_z"], c("Fx_s"), c("Fy_s"), Ts),
                                       (ss["output_z"], zero, -Fout, zero)],
                                      ss["bearings"], ss["output_z"], ss["critical"])
    cols, table = ["t"], [c("t")]
    bearing_cols = []
    out = {}
    for name, R in (("prim", Rp), ("sec", Rs)):
        for j in range(R.shape[1]):
            k = f"F_bearing_{name}{j+1}"
            out[k] = R[:, j]; bearing_cols.append(k); cols.append(k); table.append(R[:, j])
    for k, v in (("T_prim", Tp), ("M_prim", Mp), ("vm_prim", vmp),
                 ("T_sec", Tsec), ("M_sec", Ms), ("vm_sec", vms)):
        out[k] = v; cols.append(k); table.append(v)
    out.update(cols=cols, table=np.column_stack(table), bearing_cols=bearing_cols,
               kt=dict(primary=ktp, secondary=kts))
    return out
