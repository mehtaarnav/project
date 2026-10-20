"""Event scenarios: docs/theory/14_scenarios.md section 3 (T14.1-T14.4)."""
import numpy as np
import pytest

from cvtsim.cvt.system import CH, build, relax, schedule, simulate
from cvtsim.drivetrain.vehicle import G_ACC, V_SMOOTH
from cvtsim.io.params import load, load_scenario, scenario_names

c = lambda res, k: res[:, CH.index(k)]


def _run(knots, t_end=0.02, **ov):
    p = load(overrides=ov or None)
    m = build(p, scenario=dict(knots=knots) if knots is not None else None)
    relax(m)
    res, _ = simulate(m, t_end, nout=200, n_frames=2)
    return p, res


def test_t14_1_launch_scenario_is_the_launch_model():
    _, a = _run(None)
    _, b = _run(load_scenario("launch")["knots"])
    assert np.array_equal(a, b)


def test_t14_2_airborne_on_grade():
    p, res = _run([[0.0, 1.0, 10.0, 0.0], [1.0, 1.0, 10.0, 0.0]], t_end=0.05)
    v, ax = c(res, "v"), c(res, "a_x")
    veh = p["vehicle"]
    ca = 0.5 * veh["rho_air"] * veh["CdA"]
    m = veh["mass"]
    want = -(veh["Crr"] * m * G_ACC * np.tanh(v / V_SMOOTH) + ca * v * np.abs(v) + m * G_ACC * np.sin(np.radians(10))) / m
    assert np.allclose(ax[1:], want[1:], rtol=1e-9, atol=1e-12)
    assert np.all(c(res, "F_tire") == 0.0)


def test_t14_3_closed_throttle_motoring_torque():
    p, res = _run([[0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
    T_mot = p["engine"]["motoring_torque"]
    assert np.allclose(c(res, "T_eng"), -T_mot * np.tanh(c(res, "w_e") / 10.0), rtol=1e-9, atol=1e-12)


@pytest.mark.slow
def test_t14_4_energy_closure_with_schedule_changes():
    knots = [[0.0, 1.0, 0.0, 1.0], [0.3, 1.0, 0.0, 1.0], [0.35, 0.0, 5.0, 0.0], [0.5, 1.0, 5.0, 2.0], [0.6, 1.0, 0.0, 1.0]]
    _, res = _run(knots, t_end=0.6)
    E, Win, Wd = c(res, "E_mech"), c(res, "W_in"), c(res, "W_diss")
    assert np.max(np.abs((E - E[0]) - (Win - Wd))) / np.max(np.abs(Win)) < 1e-4


def test_scenario_files_are_valid():
    assert {"launch", "hill_climb", "backshift", "landing"} <= set(scenario_names())
    for name in scenario_names():
        k = schedule(name)
        assert k[0, 0] == 0.0 and k[-1, 0] >= load_scenario(name)["t_end"]


def test_bad_knots_rejected():
    with pytest.raises(ValueError):
        schedule(dict(knots=[[0.0, 1.5, 0.0, 1.0], [1.0, 1.0, 0.0, 1.0]]))
    with pytest.raises(ValueError):
        schedule(dict(knots=[[1.0, 1.0, 0.0, 1.0], [0.5, 1.0, 0.0, 1.0]]))
