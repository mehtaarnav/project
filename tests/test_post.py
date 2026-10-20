"""Post-processing tests: docs/theory/13 section 6 (T13.1-T13.12).

Targets are hand closed forms on synthetic arrays, written from the spec expressions and not
from the code under test. Tolerance rel 1e-6 unless stated.
"""
import math

import numpy as np
import pytest

from cvtsim.cvt.system import CH
from cvtsim.io.params import load
from cvtsim.post.belt_life import (
    belt_damage,
    bending_strain,
    damage_rate,
    friction_fade,
    inputs,
    temperature_factor,
    wrap_angle_and_length,
)
from cvtsim.post.efficiency import cvt_efficiency, efficiency_map, power_terms, quasi_steady_mask

REL = 1e-6


def synth(n=1, **cols):
    """Result array with the named channels set (scalars broadcast)."""
    res = np.zeros((n, len(CH)))
    for name, val in cols.items():
        res[:, CH.index(name)] = val
    return res


def base_p(**kw):
    """Synthetic input dict (same keys as inputs(params)); tests override what they vary."""
    p = dict(P_min=1.0, ratio_edges=[0.0, 2.0, 4.0], rpm_edges=[0.0, 3000.0, 6000.0],
             N_ref=1.0e6, eps_ref=0.05, k=5.0, duty=1.0, T_ref=313.15, dT_half=10.0,
             c_lo=0.001, c_hi=0.003, C=0.2, mu0=0.3, h_n=0.0025)
    p.update(kw)
    return p


# T13.1 equal sheaves: theta = pi, L = 0.05 pi
def test_t13_1_wrap_equal_sheaves():
    theta, L = wrap_angle_and_length(0.05, 0.05, 0.2)
    assert theta == pytest.approx(math.pi, rel=REL)
    assert L == pytest.approx(0.05 * math.pi, rel=REL)


# T13.2 small sheave: theta = pi - 2 asin(0.15), L = 0.02 theta
def test_t13_2_wrap_small_sheave():
    theta, L = wrap_angle_and_length(0.02, 0.05, 0.2)
    assert theta == pytest.approx(math.pi - 2.0 * math.asin(0.15), rel=REL)
    assert theta == pytest.approx(2.8404561080, rel=1e-9)
    assert L == pytest.approx(0.02 * theta, rel=REL)


# T13.3 bending strain eps = h_n / r = 0.006 / 0.02
def test_t13_3_bending_strain():
    assert bending_strain(0.006, 0.02) == pytest.approx(0.3, rel=REL)


# T13.4 damage rate, two equal sheaves, eps = 0.1, eps_ref = 0.05, k = 5, T = T_ref.
# n_dot = (40/pi) each (v = 2 m/s, L = 0.05 pi); per cycle 2^5 / N_ref.
# Ddot = 2 (40/pi) 32 / N_ref = 2560 / (pi N_ref)
def test_t13_4_damage_rate_single_strain():
    p = base_p(h_n=0.005)
    res = synth(r_eff_p=0.05, r_eff_s=0.05, w_e=40.0, w_s=40.0, T_belt=p["T_ref"])
    assert damage_rate(res, CH, p)[0] == pytest.approx(2560.0 / (math.pi * 1.0e6), rel=REL)


# T13.5 temperature factor: 2, 0.5, 1 at +10 K, -10 K, 0 K
def test_t13_5_temperature_factor():
    T_ref = 313.15
    assert temperature_factor(T_ref + 10.0, T_ref, 10.0) == pytest.approx(2.0, rel=REL)
    assert temperature_factor(T_ref - 10.0, T_ref, 10.0) == pytest.approx(0.5, rel=REL)
    assert temperature_factor(T_ref, T_ref, 10.0) == pytest.approx(1.0, rel=REL)


# T13.6 power accounting: w_e = 100, tq_p = -10, T_eng = 12, tq_s = 300 - 200, w_s = 5,
# P_slip = 150. P_eng = 1200, P_in = 1000, P_out = 500, eta = 0.5, R_bal = 350
def test_t13_6_power_and_efficiency():
    res = synth(w_e=100.0, tq_p=-10.0, T_eng=12.0, tq_s_fix=300.0, tq_s_mov=-200.0,
                w_s=5.0, P_slip_p=100.0, P_slip_s=50.0)
    t = power_terms(res, CH)
    assert t["P_eng"][0] == pytest.approx(1200.0, rel=REL)
    assert t["P_belt_in"][0] == pytest.approx(1000.0, rel=REL)
    assert t["P_out"][0] == pytest.approx(500.0, rel=REL)
    assert t["P_slip"][0] == pytest.approx(150.0, rel=REL)
    assert t["R_bal"][0] == pytest.approx(350.0, rel=REL)
    assert cvt_efficiency(res, CH, base_p())[0] == pytest.approx(0.5, rel=REL)


# T13.7 threshold: P_belt_in = 0.5 W < P_min = 1 W gives NaN
def test_t13_7_efficiency_undefined_below_threshold():
    res = synth(w_e=1.0, tq_p=-0.5, w_s=1.0, tq_s_fix=0.2, tq_s_mov=0.0)
    assert np.isnan(cvt_efficiency(res, CH, base_p())[0])


# T13.13 quasi-steady mask: P_in = 1000 W; P_slip = 1500 W excluded, P_slip = 500 W kept
def test_t13_13_quasi_steady_mask():
    res = synth(n=2, w_e=100.0, tq_p=-10.0, w_s=5.0, tq_s_fix=100.0, tq_s_mov=0.0,
                P_slip_p=np.array([1500.0, 300.0]), P_slip_s=np.array([0.0, 200.0]))
    m = quasi_steady_mask(res, CH, base_p())
    assert m.tolist() == [False, True]


# T13.8 map: instants (g = 1.0, 2000 rpm, eta = 0.8), (g = 1.5, 2500 rpm, eta = 0.6),
# (g = 3.0, 4500 rpm, eta = 0.9); edges g [0,2,4], rpm [0,3000,6000].
# Bin (0,0) = 0.7 (two instants), bin (1,1) = 0.9, bins (0,1), (1,0) empty -> NaN
def test_t13_8_efficiency_map_bins():
    rpm = np.array([2000.0, 2500.0, 4500.0])
    g = np.array([1.0, 1.5, 3.0])
    eta = np.array([0.8, 0.6, 0.9])
    w_e = rpm * math.pi / 30.0
    P_in = 10.0 * w_e
    res = synth(n=3, w_e=w_e, tq_p=-10.0, w_s=10.0, rstar_p=0.02, rstar_s=0.02 * g,
                tq_s_fix=eta * P_in / 10.0, tq_s_mov=0.0)
    mean, counts, _, _ = efficiency_map(res, CH, base_p())
    assert mean[0, 0] == pytest.approx(0.7, rel=REL)
    assert counts[0, 0] == 2
    assert mean[1, 1] == pytest.approx(0.9, rel=REL)
    assert counts[1, 1] == 1
    assert np.isnan(mean[0, 1]) and np.isnan(mean[1, 0])


# T13.9 friction fade: mu0 = 0.3, c = 0.002/K, T - T_ref = 50 K -> 0.27
def test_t13_9_friction_fade():
    p = base_p(c_lo=0.002, c_hi=0.002)
    res = synth(T_belt=p["T_ref"] + 50.0)
    f = friction_fade(res, CH, p)
    assert f["mu_lo"][0] == pytest.approx(0.27, rel=REL)
    assert f["mu_hi"][0] == pytest.approx(0.27, rel=REL)


# T13.10 constant rate, r = 0.05, w = 40 rad/s (v = 2), h_n = 0.0025 (eps = eps_ref), T = T_ref,
# t in [0, 2], N_ref = 1e6, duty = 1.
# Ddot = 80 / (pi 1e6); D = 2 Ddot; life_h = 2 / (3600 D) = 10.908308
def test_t13_10_belt_damage_life():
    p = base_p()
    t = np.array([0.0, 1.0, 2.0])
    res = synth(n=3, t=t, r_eff_p=0.05, r_eff_s=0.05, w_e=40.0, w_s=40.0, T_belt=p["T_ref"])
    out = belt_damage(res, CH, p)
    D = 2.0 * 80.0 / (math.pi * 1.0e6)
    assert out["T_run"] == pytest.approx(2.0, rel=REL)
    assert out["D"] == pytest.approx(D, rel=REL)
    assert out["life_h"] == pytest.approx(2.0 / (3600.0 * D), rel=REL)
    assert out["life_h"] == pytest.approx(10.908308, rel=1e-6)


# T13.11 duty 0.5 halves the projected life of T13.10 (life doubles: 21.816616 h)
def test_t13_11_duty_cycle_scales_life():
    p = base_p(duty=0.5)
    t = np.array([0.0, 1.0, 2.0])
    res = synth(n=3, t=t, r_eff_p=0.05, r_eff_s=0.05, w_e=40.0, w_s=40.0, T_belt=p["T_ref"])
    out = belt_damage(res, CH, p)
    assert out["life_h"] == pytest.approx(21.816616, rel=1e-6)


# T13.12 smoke on a real short launch: functions run, eta within [0, 1] where defined
@pytest.fixture(scope="module")
def launch():
    from cvtsim.cvt.system import build, relax, simulate
    params = load()
    m = build(params)
    relax(m)
    res, _ = simulate(m, 0.15, nout=100)   # belt engages within this window
    return params, res


def test_t13_12_smoke_real_launch(launch):
    params, res = launch
    p = inputs(params)
    assert np.all(np.isfinite(res))
    eta = cvt_efficiency(res, CH, p)
    ok = quasi_steady_mask(res, CH, p)
    assert ok.sum() > 0, "no quasi-steady efficiency: the run never reached belt transmission"
    assert np.all((eta[ok] >= 0.0) & (eta[ok] <= 1.0))
    mean, counts, _, _ = efficiency_map(res, CH, p)
    assert mean.shape == (len(p["ratio_edges"]) - 1, len(p["rpm_edges"]) - 1)
    out = belt_damage(res, CH, p)
    assert np.isfinite(out["D"]) and out["D"] >= 0.0
    f = friction_fade(res, CH, p)
    assert np.all(np.isfinite(f["mu_lo"])) and np.all(np.isfinite(f["mu_hi"]))
