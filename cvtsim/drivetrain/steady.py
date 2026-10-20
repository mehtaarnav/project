"""Steady state and fixed-ratio launch simulation (docs/theory/06_drivetrain.md, sections 4 and 7)."""
import numpy as np

from cvtsim.drivetrain.vehicle import G_ACC, reflected_inertia, road_load_torque


def steady_speed_linear(a, b, R, G, r_t, Crr, m, c_a):
    """Steady vehicle speed [m/s] for linear engine T = a - b omega at fixed ratio R (D3).

    Solves c_a v^2 + b k^2 v + (Crr m g - a k) = 0 with k = R G / r_t and returns the
    positive root. Engine speed is k v.
    """
    k = R * G / r_t
    A = c_a
    B = b * k * k
    C0 = Crr * m * G_ACC - a * k
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
