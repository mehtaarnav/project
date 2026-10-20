"""Steady state and fixed-ratio launch simulation (docs/theory/06_drivetrain.md, sections 4 and 7;
09 section 3 for the vehicle with tire slip and driveline efficiency)."""
import numpy as np

from cvtsim.core.jit import maybe_njit
from cvtsim.drivetrain.vehicle import (G_ACC, reflected_inertia, road_load_torque,
                                       road_load_force, tire_force, driveline_torque)


def steady_speed_linear(a, b, R, G, r_t, Crr, m, c_a, eta=1.0):
    """Steady vehicle speed [m/s] for linear engine T = a - b omega at fixed ratio R (D3).

    Solves c_a v^2 + eta b k^2 v + (Crr m g - eta a k) = 0 with k = R G / r_t (driveline
    efficiency eta, no tire slip) and returns the positive root. Engine speed is k v.
    """
    k = R * G / r_t
    A = c_a
    B = eta * b * k * k
    C0 = Crr * m * G_ACC - eta * a * k
    if A == 0.0:
        return -C0 / B
    return (-B + np.sqrt(B * B - 4.0 * A * C0)) / (2.0 * A)


def simulate_fixed_ratio(engine_torque_fn, J_e, R, p, t_end, h, omega_s0=0.0):
    """Fixed-ratio launch with classical RK4.

    Solves (J_e R^2 + J_veh + J_sec) omega_s' = T_e(R omega_s) R - T_load(omega_s).
    engine_torque_fn takes the engine speed [rad/s]. Returns (t, omega_s).
    """
    v = p["vehicle"]
    G, r_t = v["final_drive"], v["tire_radius"]
    J_veh = reflected_inertia(v["mass"], r_t, G)
    J_sec = p["secondary"]["J_shaft"]
    J_eff = J_e * R * R + J_veh + J_sec

    def rhs(w):
        return (engine_torque_fn(R * w) * R - road_load_torque(w, p)) / J_eff

    n = int(round(t_end / h))
    t = np.arange(n + 1) * h
    w = np.empty(n + 1)
    w[0] = omega_s0
    for i in range(n):
        wi = w[i]
        k1 = rhs(wi)
        k2 = rhs(wi + 0.5 * h * k1)
        k3 = rhs(wi + 0.5 * h * k2)
        k4 = rhs(wi + h * k3)
        w[i + 1] = wi + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return t, w


@maybe_njit
def _fixed_ratio_tire(eng_w, eng_T, J, R, m, r_t, G, eta, Crr, c_a, mu_t, s0, v_min, b_L, h_L, h, n):
    out = np.zeros((n + 1, 6))
    w_s = 0.0; v = 0.0; a = 0.0
    for i in range(n + 1):
        F_x, N_r, s = tire_force(w_s / G, v, a, m, r_t, mu_t, s0, v_min, b_L, h_L)
        T_ld = driveline_torque(F_x, w_s, r_t, G, eta)
        o = out[i]
        o[0] = i * h; o[1] = w_s; o[2] = v; o[3] = F_x; o[4] = s; o[5] = N_r
        wd = (np.interp(R * w_s, eng_w, eng_T) * R - T_ld) / J
        a = (F_x - road_load_force(v, m, Crr, c_a)) / m
        w_s += h * wd
        v += h * a
    return out


def simulate_fixed_ratio_tire(eng_w, eng_T, J_e, R, p, t_end, h):
    """Fixed-ratio launch of the 09 vehicle: secondary speed and vehicle speed coupled by tire slip.

    (J_e R^2 + J_shaft + J_w/G^2) w_s' = T_e(R w_s) R - T_ld(F_x), m v' = F_x - F_roll - F_aero.
    eng_w [rad/s], eng_T [N m]: engine table (np.interp, clamped). Explicit Euler, whose fixed
    point is the exact steady state. Returns dict of t, w_s, v, F_x, s, N_r arrays.
    """
    v = p["vehicle"]
    G = v["final_drive"]
    J = J_e * R * R + p["secondary"]["J_shaft"] + v["J_wheels"] / G ** 2
    out = _fixed_ratio_tire(np.asarray(eng_w, float), np.asarray(eng_T, float), J, R, v["mass"],
                            v["tire_radius"], G, v["eta"], v["Crr"], 0.5 * v["rho_air"] * v["CdA"],
                            v["mu_t"], v["s0"], v["v_min"], v["b_over_L"], v["h_cg"] / v["wheelbase"],
                            h, int(round(t_end / h)))
    return {k: out[:, j] for j, k in enumerate(("t", "w_s", "v", "F_x", "s", "N_r"))}
