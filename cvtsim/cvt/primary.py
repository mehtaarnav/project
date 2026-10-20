"""CVT primary: centrifugal flyweights on straight ramps (docs/theory/05_cvt.md, section 3)."""
import numpy as np


def _pp(p):
    pr = p["primary"]
    return pr["n_weights"], p["derived"]["m_w"], pr["kappa"], pr["r_w0"]


def flyweight_radius(x_p, p):
    """Flyweight radius r_w(x_p) = r_w0 + kappa x_p."""
    _, _, kap, rw0 = _pp(p)
    return rw0 + kap * x_p


def flyweight_axial_force(x_p, omega, p):
    """Flyweight axial force F_fw = n m_w kappa r_w(x_p) omega^2 (pushes the primary closed)."""
    n, m_w, kap, _ = _pp(p)
    return n * m_w * kap * flyweight_radius(x_p, p) * omega ** 2


def spring_force(x_p, p):
    """Primary spring force F_p0 + k_p x_p (opposes closure)."""
    pr = p["primary"]
    return pr["F_preload"] + pr["k_spring"] * x_p


def engagement_omega(p):
    """Static engagement speed [rad/s] from F_fw(0, omega) = F_p0 (P3)."""
    n, m_w, kap, rw0 = _pp(p)
    return np.sqrt(p["primary"]["F_preload"] / (n * m_w * kap * rw0))


def negative_stiffness_omega(p):
    """Speed [rad/s] above which the net primary axial stiffness is negative (P4)."""
    n, m_w, kap, _ = _pp(p)
    return np.sqrt(p["primary"]["k_spring"] / (n * m_w * kap ** 2))


def primary_mass(x_p, p):
    """Return (J_rot, m_ax): rotary inertia about the shaft and axial mass at x_p."""
    n, m_w, kap, _ = _pp(p)
    rw = flyweight_radius(x_p, p)
    J_rot = p["engine"]["J_e"] + n * m_w * rw ** 2
    m_ax = p["primary"]["m_axial"] + n * m_w * kap ** 2
    return J_rot, m_ax


def primary_coriolis_torque(x_p, xdot_p, omega, p):
    """Left-hand term 2 n m_w r_w kappa xdot omega of (P1), i.e. (dJ/dt) omega.

    The speed equation is J omega_dot = -(this value) + T_eng + tau_belt - c omega.
    """
    n, m_w, kap, _ = _pp(p)
    return 2.0 * n * m_w * flyweight_radius(x_p, p) * kap * xdot_p * omega


def flyweight_force_profile(u, omega, p):
    """Static flyweight force F_fw(u, omega) = n m_w r_w(u) kappa(u) omega^2 with the ramp
    profile (09 section 1); u is travel from full open."""
    r_w, kap, _ = p["derived"]["ramp"](u)
    return p["primary"]["n_weights"] * p["derived"]["m_w"] * r_w * kap * omega ** 2


def engagement_omega_profile(p):
    """Static engagement speed [rad/s] from F_fw(0, omega) = F_p0 with the ramp profile."""
    return np.sqrt(p["primary"]["F_preload"] / flyweight_force_profile(0.0, 1.0, p))


def negative_stiffness_omega_profile(u, p):
    """Speed [rad/s] above which k_p - dF_fw/du < 0 at travel u, with
    dF_fw/du = n m_w (kappa^2 + r_w kappa') omega^2."""
    r_w, kap, dkap = p["derived"]["ramp"](u)
    n, m_w = p["primary"]["n_weights"], p["derived"]["m_w"]
    return np.sqrt(p["primary"]["k_spring"] / (n * m_w * (kap * kap + r_w * dkap)))
