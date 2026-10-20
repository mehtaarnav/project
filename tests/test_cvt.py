"""CVT numeric test targets T5.1-T5.10 (docs/theory/05_cvt.md, section 8)."""
import math

import numpy as np
import pytest

from cvtsim.io.params import load
from cvtsim.cvt.geometry import derive_geometry, ratio_from_xp
from cvtsim.cvt.primary import (engagement_omega, negative_stiffness_omega,
                                flyweight_axial_force, spring_force, primary_mass,
                                primary_coriolis_torque)
from cvtsim.cvt.secondary import secondary_static_force

RPM = 2.0 * math.pi / 60.0


def to_rpm(omega):
    return omega / RPM


@pytest.fixture(scope="module")
def p():
    return load()


def test_T5_1_derive_geometry_beta45():
    g = derive_geometry(np.radians(45.0), 0.02548, 0.02448, 3.0, 0.8)
    targets = dict(rp_min=0.0101964, rp_max=0.0229364, rs_max=0.0305891,
                   rs_min=0.0183491, C=0.251495, L=0.632776)
    for key, val in targets.items():
        assert g[key] == pytest.approx(val, rel=1e-5), key


def test_T5_2_engagement_rpm_kappa_half(p):
    assert to_rpm(engagement_omega(p)) == pytest.approx(1725.66, rel=1e-5)


def test_T5_3_engagement_rpm_kappa_one(p):
    p1 = load(overrides={"primary.kappa": 1.0})
    assert to_rpm(engagement_omega(p1)) == pytest.approx(1220.23, rel=1e-5)


def test_T5_4_negative_stiffness_onset_kappa_half(p):
    assert to_rpm(negative_stiffness_omega(p)) == pytest.approx(2262.57, rel=1e-5)


def net_restoring_stiffness(p, rpm_val, x0=0.01, h=1e-7):
    """k_net = d(F_spring - F_fw)/dx_p by central difference (P4).

    The restoring force on the movable half is F_spring - F_fw, so
    k_net = k_p - n m_w kappa^2 omega^2.
    """
    w = rpm_val * RPM

    def restoring(x):
        return spring_force(x, p) - flyweight_axial_force(x, w, p)

    return (restoring(x0 + h) - restoring(x0 - h)) / (2.0 * h)


def test_T5_5_k_net_sign(p):
    assert net_restoring_stiffness(p, 2300.0) < 0.0
    assert net_restoring_stiffness(p, 2200.0) > 0.0


def test_T5_6_secondary_torsion_axial_preload(p):
    F_t0 = p["secondary"]["k_torsion"] * p["derived"]["theta_t0"] / p["derived"]["p_h"]
    assert F_t0 == pytest.approx(54.009, rel=1e-4)
    assert p["derived"]["F_torsion0"] == pytest.approx(54.009, rel=1e-4)


def test_T5_7_secondary_static_force_torque_20(p):
    assert secondary_static_force(0.0, 20.0, p) == pytest.approx(698.097, rel=1e-5)


def test_T5_8_secondary_static_force_stroke_10mm(p):
    assert secondary_static_force(0.01, 0.0, p) == pytest.approx(496.235, rel=1e-4)


def test_T5_9_energy_conservation_free_primary(p):
    """Free primary (P1-P2, T_eng = 0, no belt, no damping, no stops) under RK4 for 1 s."""
    k, Fp0 = p["primary"]["k_spring"], p["primary"]["F_preload"]

    def energy(state):
        w, x, v = state
        _, m_ax = primary_mass(x, p)
        J, _ = primary_mass(x, p)
        return 0.5 * J * w * w + 0.5 * m_ax * v * v + Fp0 * x + 0.5 * k * x * x

    def rhs(state):
        w, x, v = state
        J, m_ax = primary_mass(x, p)
        coriolis = primary_coriolis_torque(x, v, w, p)  # (dJ/dt) omega
        w_dot = -coriolis / J
        x_ddot = (flyweight_axial_force(x, w, p) - spring_force(x, p)) / m_ax
        return np.array([w_dot, v, x_ddot])

    state = np.array([2000.0 * RPM, 0.01, 0.0])  # 2000 rpm, x_p inside the stroke
    E0 = energy(state)
    h, n_steps = 1e-4, 10000
    for _ in range(n_steps):
        k1 = rhs(state)
        k2 = rhs(state + 0.5 * h * k1)
        k3 = rhs(state + 0.5 * h * k2)
        k4 = rhs(state + h * k3)
        state = state + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    E1 = energy(state)
    assert abs(E1 - E0) / abs(E0) < 1e-6


def test_T5_10_ratio_limits(p):
    g = p["derived"]["geometry"]
    assert ratio_from_xp(0.0, g) == pytest.approx(3.0, rel=1e-4)
    assert ratio_from_xp(p["primary"]["stroke"], g) == pytest.approx(0.8, rel=1e-4)
