"""Fast CVT model: docs/theory/15_fast_cvt.md section 6 (T15.1-T15.7)."""
import math

import numpy as np
import pytest

from cvtsim.fast.model import FCH, build_fast, capacity, geometry_table, simulate_fast, stick_force
from cvtsim.io.params import load


@pytest.fixture(scope="module")
def p():
    return load(tune="default")


def test_t15_1_geometry_ends(p):
    g = p["derived"]["geometry"]
    x, rs, php, phs = geometry_table(g)
    assert rs[0] == pytest.approx(g["rs_low"], rel=1e-9)
    assert rs[0] / g["rp_min"] == pytest.approx(g["ratio_low"], rel=1e-9)
    assert rs[-1] / g["rp_max"] == pytest.approx(p["ratio"]["high"], rel=1e-9)
    assert np.allclose(php + phs, 2 * np.pi)


def test_t15_2_equal_radii_wrap():
    C, r = 0.232, 0.06
    php = math.pi - 2 * math.asin((r - r) / C)
    assert php == math.pi


def test_t15_3_capacity():
    mu_eff = 0.4 / math.sin(math.radians(12))
    assert capacity(500.0, math.pi, mu_eff) == pytest.approx(995.2675501957627, rel=1e-12)


def test_t15_4_locked_belt_acceleration():
    r_p, r_s, T, J_e, J_s = 0.04, 0.09, 20.0, 0.05, 0.03
    F = stick_force(r_p, r_s, T, 0.0, J_e, J_s, 0.0, 0.0, 300.0, 300.0 * r_p / r_s)
    wd_e, wd_s = (T - F * r_p) / J_e, F * r_s / J_s
    assert r_p * wd_e == pytest.approx(r_s * wd_s, rel=1e-12)
    assert wd_e == pytest.approx(T * r_s ** 2 / (J_e * r_s ** 2 + J_s * r_p ** 2), rel=1e-12)


def test_t15_6_energy_closure(p):
    r = simulate_fast(build_fast(p), 8.0)
    c = lambda k: r[:, FCH.index(k)]
    E, Wi, Wd = c("E_kin"), c("W_in"), c("W_diss")
    assert np.max(np.abs((E - E[0]) - (Wi - Wd))) / Wi[-1] < 1e-3


def test_t15_7_step_convergence(p):
    fm = build_fast(p)
    a, b = simulate_fast(fm, 8.0, h=2e-4), simulate_fast(fm, 8.0, h=1e-4)
    for k in ("w_e", "v"):
        i = FCH.index(k)
        assert b[-1, i] == pytest.approx(a[-1, i], rel=5e-3)


def test_runs_every_scenario_quickly(p):
    import time
    from cvtsim.io.params import load_scenario, scenario_names
    for name in scenario_names():
        sc = load_scenario(name)
        t0 = time.time()
        r = simulate_fast(build_fast(p, scenario=sc), sc["t_end"])
        assert np.all(np.isfinite(r)) and time.time() - t0 < 5.0


FULL_RUNS = ["runs/20261020-211533-tune-default", "runs/20261020-213705-tune-default-hill_climb",
             "runs/20261025-201504-tune-default-backshift", "runs/20261025-201508-tune-cvtech_0651_1100"]


@pytest.mark.parametrize("run", FULL_RUNS[2:])
def test_calibration_validates_on_unseen_runs(run):
    """Section 5 acceptance checks on full-model runs not used in the fit (needs their series)."""
    from pathlib import Path
    from cvtsim.fast.calibrate import checks, load_full, run_on
    from cvtsim.fast.model import load_calibration
    if not (Path(run) / "series.csv.gz").exists():
        pytest.skip("full-model series not present (git-ignored)")
    full = load_full(run)
    c = checks(full, run_on(full, load_calibration()))
    assert c["rpm_ok"] and c["xp_ok"] and c["v_ok"] and c["engage_ok"], c
