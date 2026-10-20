"""Drivetrain numeric test targets T6.1-T6.8 (docs/theory/06_drivetrain.md, section 8)."""
import pytest

from cvtsim.io.params import load
from cvtsim.drivetrain.engine import EngineCurve, RPM_TO_RAD
from cvtsim.drivetrain.vehicle import reflected_inertia, road_load_torque
from cvtsim.drivetrain.steady import steady_speed_linear, simulate_fixed_ratio


@pytest.fixture(scope="module")
def p():
    return load()


@pytest.fixture(scope="module")
def engine(p):
    return EngineCurve.from_params(p)


def rad(rpm):
    return rpm * RPM_TO_RAD


def test_T6_1_engine_torque_grid_points(engine):
    assert engine.torque(rad(2400.0)) == pytest.approx(25.08, abs=1e-12)
    assert engine.torque(rad(3000.0)) == pytest.approx(22.51, abs=1e-12)
    assert engine.torque(rad(3600.0)) == pytest.approx(18.30, abs=1e-12)


def test_T6_2_engine_torque_midpoint(engine):
    assert engine.torque(rad(2900.0)) == pytest.approx((23.59 + 22.51) / 2, abs=1e-12)


def test_T6_3_engine_torque_below_chart_flagged(engine):
    assert engine.torque(rad(2000.0)) == pytest.approx(25.08, abs=1e-12)
    assert engine.flag_unverified(rad(2000.0)) is True
    assert engine.flag_unverified(rad(2400.0)) is False


def test_T6_4_engine_torque_governor_taper(engine):
    assert engine.torque(rad(3700.0)) == pytest.approx(9.15, abs=1e-12)
    assert engine.torque(rad(3900.0)) == pytest.approx(0.0, abs=1e-12)


def test_T6_5_reflected_inertia():
    assert reflected_inertia(250.0, 0.2921, 8.1) == pytest.approx(0.325112, rel=1e-6)


def test_T6_6_steady_speed_linear_engine():
    v = steady_speed_linear(40.0, 0.05, 1.0, 8.1, 0.2921, 0.015, 250.0, 0.5)
    assert v == pytest.approx(21.744, rel=1e-4)
    engine_rpm = 1.0 * 8.1 * v / 0.2921 / RPM_TO_RAD
    assert engine_rpm == pytest.approx(5757.91, rel=1e-4)


def test_T6_7_simulate_fixed_ratio_reaches_steady_speed():
    # c_a = 0.5 = (1/2) rho_air CdA, so CdA = 2 c_a / rho_air
    p_c = load(overrides={"vehicle.CdA": 2.0 * 0.5 / 1.2, "vehicle.Crr": 0.015})  # spec 06 inputs
    t, w = simulate_fixed_ratio(lambda we: 40.0 - 0.05 * we, 0.05, 1.0, p_c,
                                60.0, 1e-3)
    v_final = w[-1] * p_c["vehicle"]["tire_radius"] / p_c["vehicle"]["final_drive"]
    assert v_final == pytest.approx(21.744, rel=1e-3)


def test_T6_8_rolling_resistance_torque_at_10_mps():
    p = load(overrides={"vehicle.Crr": 0.015, "vehicle.CdA": 0.0})  # spec 06: rolling resistance only
    omega_s = 10.0 * p["vehicle"]["final_drive"] / p["vehicle"]["tire_radius"]
    assert road_load_torque(omega_s, p) == pytest.approx(1.32662, rel=1e-6)
