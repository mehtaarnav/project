"""CVT secondary: helix cam, torsion and linear springs (docs/theory/05_cvt.md, section 4)."""
import numpy as np

SIGMA = -1.0  # opening rotates the movable half backward (torque feedback)


def secondary_static_force(x_s, Tq, p):
    """Static closing force F_s(x_s, Tq) (S2) for shaft torque Tq [N m]."""
    s, d = p["secondary"], p["derived"]
    ph, th0 = d["p_h"], d["theta_t0"]
    return (d["F_l0"] + s["k_linear"] * x_s
            + (Tq / 2.0 + s["k_torsion"] * (th0 + x_s / ph)) / ph)


def secondary_mass_matrix(p):
    """Mass matrix in (omega_s, xdot_s) coordinates (S1)."""
    s, d = p["secondary"], p["derived"]
    ph = d["p_h"]
    Js, Jm, mm = s["J_shaft"], s["J_movable"], s["m_axial"]
    return np.array([[Js + Jm, SIGMA * Jm / ph],
                     [SIGMA * Jm / ph, mm + Jm / ph ** 2]])


def secondary_spring_force(x_s, p):
    """Spring part of the generalized force Q_x (S1), as a contribution to Q_x.

    Equals -(F_l0 + k_l x_s) - k_t (theta_t0 + x_s/p_h)/p_h. Belt thrust, cam torque,
    damping and stop forces are excluded.
    """
    s, d = p["secondary"], p["derived"]
    ph, th0 = d["p_h"], d["theta_t0"]
    return -(d["F_l0"] + s["k_linear"] * x_s
             + s["k_torsion"] * (th0 + x_s / ph) / ph)
