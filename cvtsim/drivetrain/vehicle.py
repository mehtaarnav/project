"""Vehicle inertia and road load reflected to the secondary shaft (docs/theory/06_drivetrain.md,
section 3), and the tire/driveline kernels of the vehicle with its own speed state (09 section 3)."""
import numpy as np

from cvtsim.core.jit import maybe_njit

G_ACC = 9.81       # gravity [m/s^2]
V_SMOOTH = 0.05    # tanh smoothing speed for the rolling-resistance sign [m/s]


def reflected_inertia(m, r_t, G):
    """J_veh = m r_t^2 / G^2 at the secondary shaft (D1)."""
    return m * r_t ** 2 / G ** 2


def vehicle_speed(omega_s, p):
    """Vehicle speed v = omega_s r_t / G for secondary speed omega_s [rad/s]."""
    v = p["vehicle"]
    return omega_s * v["tire_radius"] / v["final_drive"]


def road_load_torque(omega_s, p):
    """Road-load torque at the secondary shaft [N m] (D2)."""
    v = p["vehicle"]
    vel = vehicle_speed(omega_s, p)
    F = (v["Crr"] * v["mass"] * G_ACC * np.tanh(vel / V_SMOOTH)
         + 0.5 * v["rho_air"] * v["CdA"] * vel ** 2)
    return F * v["tire_radius"] / v["final_drive"]


@maybe_njit
def road_load_force(v, m, Crr, c_a):
    """Rolling plus aerodynamic resistance F_roll + F_aero [N] at speed v (09 section 3)."""
    return Crr * m * G_ACC * np.tanh(v / V_SMOOTH) + c_a * v * abs(v)


@maybe_njit
def tire_force(w_w, v, a_prev, m, r_t, mu_t, s0, v_min, b_L, h_L):
    """Rear-tire force (F_x, N_r, s) for wheel speed w_w and vehicle speed v (09 section 3).

    N_r = m g b/L + m a h_cg/L with the acceleration a lagged one step, clipped at 0;
    s = (w_w r_t - v)/max(|v|, v_min); F_x = mu_t N_r tanh(s/s0).
    """
    N_r = m * G_ACC * b_L + m * a_prev * h_L
    if N_r < 0.0:
        N_r = 0.0
    s = (w_w * r_t - v) / max(abs(v), v_min)
    return mu_t * N_r * np.tanh(s / s0), N_r, s


@maybe_njit
def driveline_torque(F_x, w_s, r_t, G, eta):
    """Load torque on the secondary for wheel force F_x through efficiency eta (09 section 3)."""
    if F_x * w_s >= 0.0:
        return F_x * r_t / G / eta
    return eta * F_x * r_t / G
