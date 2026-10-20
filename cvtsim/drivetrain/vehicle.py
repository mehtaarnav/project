"""Vehicle inertia and road load reflected to the secondary shaft (docs/theory/06_drivetrain.md, section 3)."""
import numpy as np

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
