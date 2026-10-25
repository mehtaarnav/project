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
                         SU, 500e6, "ground", win=3)
    assert f["N_rev"] == 0.0 and f["D_bend"] == 0.0


def test_goodman_on_simultaneous_stresses():
    """7.3: per sample, bending at that sample; torque half-range and mean over the window."""
    rng = np.random.default_rng(4)
    t = np.linspace(0.0, 1.0, 41)
    M, T = 60 * rng.random(41), 80 * rng.random(41) - 20
    d, win = 0.025, 5
    f = du.shaft_fatigue(t, np.full_like(t, 50.0), M, T, dict(d=d, r=0.002), (2.0, 1.5), SU, 500e6, "ground", win)
    zb, zt = 32 / (math.pi * d ** 3), 16 / (math.pi * d ** 3)
    Se, Kf = f["Se_MPa"] * 1e6, f["Kf"]
    best = math.inf
    for i in range(41):
        seg = T[max(0, i - 2):i + 3]
        pad = [T[0]] * max(0, 2 - i) + list(seg) + [T[-1]] * max(0, i + 3 - 41)
        ta, tm = 0.5 * (max(pad) - min(pad)), sum(pad) / 5
        sa = math.sqrt((Kf * zb * M[i]) ** 2 + 3 * (1.5 * zt * ta) ** 2)
        sm = math.sqrt(3) * 1.5 * zt * abs(tm)
        best = min(best, 1 / (sa / Se + sm / SU))
    assert f["n_goodman"] == pytest.approx(best, rel=1e-9)
    vm = np.sqrt((zb * M) ** 2 + 3 * (zt * T) ** 2).max()
    assert f["n_yield"] == pytest.approx(500e6 / vm, rel=1e-12)


def test_primary_section_carries_pto_torque():
    """The crank PTO section is loaded by T_pto (belt torque + clutch inertia + pulses), not tq_p."""
    from cvtsim.cvt.system import CH
    from cvtsim.fe.shaft_loads import shaft_histories
    from cvtsim.io.params import load
    res = np.zeros((3, len(CH)))
    res[:, CH.index("t")] = [0.0, 0.1, 0.2]
    res[:, CH.index("T_pto")] = [40.0, 250.0, 40.0]
    res[:, CH.index("tq_p")] = -40.0
    sh = shaft_histories(res, load(), CH)
    assert np.allclose(np.abs(sh["T_prim"]), [40.0, 250.0, 40.0])
