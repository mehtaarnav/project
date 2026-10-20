"""Launch durability: docs/theory/12_fatigue_bearings.md section 7 (T12.15-T12.19)."""
import math

import numpy as np
import pytest

from cvtsim.post import durability as du

SU, SE = 600e6, 300e6


def test_t12_15_damage_by_revolutions(monkeypatch):
    t = np.linspace(0.0, 1.0, 1001)
    w = np.full_like(t, 2 * np.pi * 10)
    d = 0.025
    M = np.full_like(t, 400e6 * math.pi * d ** 3 / 32)  # sigma_a = 400 MPa
    monkeypatch.setattr(du, "endurance_limit", lambda *a, **k: SE)
    f = du.shaft_fatigue(t, w, M, np.zeros_like(t), dict(d=d, r=0.002), (1.0, 1.0), SU, 500e6, "ground")
    assert f["Kf"] == 1.0
    assert f["N_rev"] == pytest.approx(10.0, rel=1e-12)
    assert f["D"] == pytest.approx(10 / 34017.4381512233, rel=1e-6)
    assert f["D_tors"] == 0.0


def test_t12_16_launches_to_l10():
    t = np.linspace(0.0, 1.0, 101)
    w = np.full_like(t, 2 * np.pi * 10)
    b = du.bearing_life(t, w, np.full_like(t, 1000.0), C=10000.0)
    assert b["P_eq_N"] == pytest.approx(1000.0, rel=1e-12)
    assert b["launches_to_L10"] == pytest.approx(1e8, rel=1e-9)


def test_t12_17_required_rating():
    t = np.linspace(0.0, 60.0, 61)
    w = np.full_like(t, 1000 * np.pi / 30)
    b = du.bearing_life(t, w, np.full_like(t, 1000.0), target_hours=100.0, target_launches=6e6 / 1000.0)
    assert b["rpm_mean"] == pytest.approx(1000.0, rel=1e-12)
    assert b["C_req_hours_N"] == pytest.approx(1817.1205928321397, rel=1e-6)
    assert b["C_req_N"] == pytest.approx(1817.1205928321397, rel=1e-6)


def test_t12_18_equivalent_load_by_revolutions():
    t = np.array([0.0, 1.0, 1.0 + 1e-12, 2.0])
    w = np.full_like(t, 2 * np.pi)
    R = np.array([1000.0, 1000.0, 2000.0, 2000.0])
    b = du.bearing_life(t, w, R)
    assert b["P_eq_N"] == pytest.approx(1650.9636244473127, rel=1e-6)


def test_t12_19_side_load():
    s = du.side_load(np.full(5, 300.0), np.full(5, 400.0), 0.06009, F_allow=1000.0)
    assert s["F_hub_peak_N"] == pytest.approx(500.0)
    assert s["M_face_peak_Nm"] == pytest.approx(30.045, rel=1e-9)
    assert s["side_load_ratio"] == pytest.approx(0.5)


def test_sampling_flag():
    t = np.arange(0, 1, 1e-4)
    assert du.sampling_ok(t, 3600.0, f_pto=187.0)
    assert not du.sampling_ok(np.arange(0, 1, 2.7e-3), 3600.0)


def test_torsion_family_counts_rainflow():
    t = np.linspace(0.0, 1.0, 7)
    T = np.array([0.0, 50.0, -50.0, 50.0, -50.0, 0.0, 0.0])
    f = du.shaft_fatigue(t, np.zeros_like(t), np.zeros_like(t), T, dict(d=0.025, r=0.002), (2.0, 1.5),
                         SU, 500e6, "ground")
    assert f["N_rev"] == 0.0 and f["D_bend"] == 0.0
    zt = 16 / (math.pi * 0.025 ** 3)
    assert f["sa_env_MPa"] == pytest.approx(math.sqrt(3) * 1.5 * 50 * zt / 1e6, rel=1e-12)  # 7.3, M = 0
