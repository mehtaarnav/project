"""Fatigue and bearing tests: docs/theory/12_fatigue_bearings.md section 4 (T12.1-T12.14).

Targets are hand closed forms; tolerances are as stated in the spec (rel 1e-5 unless noted).
"""
import math
from pathlib import Path

import pytest
import yaml

from cvtsim.fe.bearings import equivalent_load, l10_hours, l10_revolutions, pto_side_load_ratio
from cvtsim.fe.fatigue import (
    basquin_constants,
    bending_stress_amplitude,
    cycles_to_failure,
    endurance_limit,
    fatigue_kf,
    goodman_n,
    marin_ka,
    marin_kb,
    miner_damage,
    notch_sensitivity,
    peterson_a_mm,
    rainflow,
    shaft_cycle_stresses,
    yield_n,
)

ROOT = Path(__file__).resolve().parents[1]
REL = 1e-5


def _close(a, b, rel=REL):
    assert math.isclose(a, b, rel_tol=rel, abs_tol=0.0), f"{a} vs {b}"


def test_t12_1_rainflow_six_point():
    cyc = sorted(rainflow([0, 10, -10, 10, -10, 0]))
    expected = sorted([(20.0, 0.0, 1.0), (10.0, 5.0, 0.5), (20.0, 0.0, 0.5), (10.0, -5.0, 0.5)])
    assert len(cyc) == len(expected)
    for got, exp in zip(cyc, expected):
        assert got == pytest.approx(exp, abs=0)


def test_t12_2_rainflow_sine_three_cycles():
    cyc = rainflow([10, -10, 10, -10, 10, -10, 10])
    assert sum(n for _, _, n in cyc) == pytest.approx(3.0, abs=0)
    assert all(r == pytest.approx(20.0, abs=0) for r, _, _ in cyc)
    assert all(m == pytest.approx(0.0, abs=0) for _, m, _ in cyc)


def test_t12_3_bending_stress_matches_shaft_target():
    _close(bending_stress_amplitude(100.0, 0.025, 1.0) / 1e6, 65.1899, rel=1e-5)


def test_t12_4_peterson_r1_su1000():
    _close(peterson_a_mm(1000.0), 0.0940980, rel=1e-5)
    q = notch_sensitivity(1.0, 1000.0)
    _close(q, 0.913995, rel=1e-5)
    _close(fatigue_kf(2.0, q), 1.913995, rel=1e-5)


def test_t12_5_peterson_r2_su1000():
    q = notch_sensitivity(2.0, 1000.0)
    _close(q, 0.955065, rel=1e-5)
    _close(fatigue_kf(2.0, q), 1.955065, rel=1e-5)


def test_t12_6_goodman():
    _close(goodman_n(100.0, 50.0, 300.0, 600.0), 2.4, rel=1e-12)


def test_t12_7_yield():
    _close(yield_n(500.0, 100.0, 50.0), 3.333333333, rel=1e-8)


def test_t12_8_l10_constant_load():
    _close(l10_revolutions(10.0, 1.0, 3.0), 1e9, rel=1e-12)
    _close(l10_hours(10.0, 1.0, [1000.0] * 10, 3.0), 16666.666667, rel=1e-8)


def test_t12_9_equivalent_load_two_level():
    _close(equivalent_load([1000.0, 2000.0], [0.75, 0.25], 3.0), 1401.0197, rel=1e-6)


def test_t12_10_marin_and_endurance():
    _close(marin_ka(600.0, "machined"), 0.827878, rel=1e-5)
    _close(marin_kb(25.0), 0.878703, rel=1e-5)
    _close(endurance_limit(600.0, "machined", 25.0) / 1e6, 218.2377, rel=1e-5)


def test_t12_11_basquin_constants():
    A, b = basquin_constants(600e6, 300e6)
    _close(A / 1e6, 972.0, rel=1e-9)
    _close(b, -0.0850908, rel=1e-6)
    _close(cycles_to_failure(300e6, 600e6, 300e6), 1e6, rel=1e-9)


def test_t12_12_miner_single_level():
    D = miner_damage([(400e6, 0.0, 1000.0)], 600e6, 300e6)
    _close(D, 0.0293967, rel=1e-5)
    _close(cycles_to_failure(400e6, 600e6, 300e6), 34017.44, rel=1e-5)


def test_t12_12_miner_endurance_cutoff():
    assert miner_damage([(250e6, 0.0, 1e9)], 600e6, 300e6) == 0.0


def test_t12_13_pto_side_load_ratio():
    assert pto_side_load_ratio([100.0, -500.0, 300.0], 1000.0) == pytest.approx(0.5, abs=0)


def test_t12_14_shaft_full_cycle_von_mises():
    cycles = shaft_cycle_stresses(100.0, [0, 50, -50, 50, -50, 0], 0.025, 1.0, 1.0)
    full = [c for c in cycles if c[2] == 1.0]
    assert len(full) == 1
    sa, sm, _ = full[0]
    _close(sa / 1e6, 71.0390, rel=1e-5)
    assert sm == pytest.approx(0.0, abs=1e-9)


def test_rainflow_monotone_residue():
    assert sum(n for _, _, n in rainflow([0, 1, 2, 3])) == pytest.approx(0.5, abs=0)


def test_materials_yaml_has_sources():
    data = yaml.safe_load((ROOT / "cvtsim" / "io" / "materials.yaml").read_text())
    tags = ("BRIEF", "measured", "DERIVED", "ASSUMED", "UNVERIFIED", "PLACEHOLDER")
    for group in ("materials", "bearings"):
        for name, entry in data[group].items():
            for key, val in entry.items():
                if isinstance(val, dict) and "source" in val:
                    assert any(val["source"].startswith(t) for t in tags), f"{group}.{name}.{key}"
    assert data["pto"]["allowable_side_load_N"]["source"].startswith("PLACEHOLDER")
