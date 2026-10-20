"""Shaft tests: docs/theory/07_shafts.md section 7 (T7.1-T7.11).

Rigid supports: the spec's closed forms assume rigid bearings. Bearing springs
at the default k_r = 1e12 N/m add ~2.9e-5 relative tip error in T7.5 (rigid
spring compliance), so the tests pass k_r = RIGID = 1e15 N/m explicitly.
"""
import math

import numpy as np
import pytest

from cvtsim.fe.peterson import kt_hole_plate_tension, kt_shoulder_bending, kt_shoulder_torsion
from cvtsim.fe.shaft import Shaft, cowper_kappa, first_bending_frequency, influence_matrix

E = 200e9
NU = 0.3
D = 0.025
RIGID = 1e15
G = E / (2.0 * (1.0 + NU))
KAPPA = 6.0 * (1.0 + NU) / (7.0 + 6.0 * NU)  # 0.886364 (Cowper)


def _ss_shaft(kappa=None):
    s = Shaft([0.0, 0.15, 0.3], [D, D], E, NU, kappa=kappa)
    s.add_bearing(0.0, k_r=RIGID)
    s.add_bearing(0.3, k_r=RIGID)
    return s


def test_T7_1_ss_beam_reactions():
    # P = 1000 N directed -y at midspan -> bearings push +y: 500 N each.
    r = _ss_shaft().solve([(0.15, 0.0, -1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert r.reactions[0, 1] == pytest.approx(500.0, rel=1e-9)
    assert r.reactions[1, 1] == pytest.approx(500.0, rel=1e-9)


def test_T7_2_ss_midspan_deflection_timoshenko():
    r = _ss_shaft().solve([(0.15, 0.0, -1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert abs(r.u[1, 1]) == pytest.approx(1.48918e-4, rel=1e-5)


def test_T7_3_ss_midspan_deflection_euler_bernoulli_limit():
    r = _ss_shaft(kappa=1e9).solve([(0.15, 0.0, -1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert abs(r.u[1, 1]) == pytest.approx(1.46677e-4, rel=1e-5)


def test_T7_4_overhung_reactions():
    s = Shaft([0.0, 0.2, 0.25], [D, D], E, NU)
    s.add_bearing(0.0, k_r=RIGID)
    s.add_bearing(0.2, k_r=RIGID)
    r = s.solve([(0.25, 0.0, 1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert r.reactions[0, 1] == pytest.approx(250.0, rel=1e-9)
    assert r.reactions[1, 1] == pytest.approx(-1250.0, rel=1e-9)


def test_T7_5_overhung_tip_deflection_timoshenko():
    s = Shaft([0.0, 0.2, 0.25], [D, D], E, NU)
    s.add_bearing(0.0, k_r=RIGID)
    s.add_bearing(0.2, k_r=RIGID)
    r = s.solve([(0.25, 0.0, 1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert r.u[2, 1] == pytest.approx(5.61923e-5, rel=1e-5)


def test_T7_6_torsional_twist():
    s = Shaft([0.0, 0.3], [D], E, NU)
    s.add_bearing(0.0)  # torsional ground added at this bearing by the solver
    r = s.solve([(0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 50.0)])
    twist = r.u[1, 5] - r.u[0, 5]
    assert twist == pytest.approx(5.08481e-3, rel=1e-6)


def test_T7_7_stresses_M100_T50():
    # Bearings at 0 and 0.3 carry F/2 each; F = 4000/3 N gives M = 100 N.m at z = 0.15.
    # Torque 50 N.m applied at the far end, reacted at the z = 0 bearing.
    s = Shaft([0.0, 0.3], [D], E, NU)
    s.add_bearing(0.0, k_r=RIGID)
    s.add_bearing(0.3, k_r=RIGID)
    r = s.solve([(0.15, 0.0, -4000.0 / 3.0, 0.0, 0.0, 0.0, 0.0),
                 (0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 50.0)])
    sig_b, tau, sig_a, sig_vm = r.stress(0.15)
    assert sig_b == pytest.approx(65.1899e6, rel=1e-5)
    assert tau == pytest.approx(16.2975e6, rel=1e-5)
    assert sig_a == pytest.approx(0.0, abs=1e-6)
    assert sig_vm == pytest.approx(71.0390e6, rel=1e-5)


def test_T7_8_shoulder_bending_fit_value():
    assert kt_shoulder_bending(1.2, 1.0, 0.1) == pytest.approx(1.66575, abs=1e-4)


def test_T7_9_hole_plate_limit_is_3():
    assert kt_hole_plate_tension(1e-6, 1.0) == pytest.approx(3.004, abs=1e-3)


def test_T7_10_shoulder_factors_decrease_with_radius_and_exceed_one():
    r_b = np.linspace(0.05, 0.5, 46)                 # bending valid for t/r in [0.1, 2]
    kb = np.array([kt_shoulder_bending(1.2, 1.0, r) for r in r_b])
    assert np.all(np.diff(kb) < 0) and np.all(kb > 1.0)
    r_t = np.linspace(0.05, 0.4, 36)                 # torsion valid for t/r in [0.25, 4]
    kt = np.array([kt_shoulder_torsion(1.2, 1.0, r) for r in r_t])
    assert np.all(np.diff(kt) < 0) and np.all(kt > 1.0)


def test_T7_10_out_of_range_raises():
    with pytest.raises(ValueError):
        kt_shoulder_torsion(1.2, 1.0, 0.5)           # t/r = 0.2 < 0.25
    with pytest.raises(ValueError):
        kt_shoulder_bending(1.2, 1.0, 0.01)          # t/r = 10 > 2


def test_T7_11_equilibrium_random_load_set():
    rng = np.random.default_rng(12345)
    zb = [0.0, 0.11, 0.37]
    kz = [0.0, 0.05, 0.1, 0.2, 0.29, 0.37, 0.4]
    s = Shaft([0.0, 0.2, 0.4], [0.02, 0.025], E, NU)
    s.add_bearing(zb[0], k_r=RIGID, k_a=RIGID)     # thrust bearing
    s.add_bearing(zb[1], k_r=RIGID)
    s.add_bearing(zb[2], k_r=RIGID)
    # torsional ground at the first bearing is added automatically
    loads = []
    for z in kz:
        loads.append((z, *rng.uniform(-1000, 1000, 3), *rng.uniform(-100, 100, 3)))
    r = s.solve(loads)
    L = np.array(loads)
    R = np.zeros((len(zb), 7))
    R[:, 0] = zb
    R[:, 1:4] = r.reactions[:, :3]
    R[:, 6] = r.reactions[:, 3]
    A = np.vstack([L, R])                      # columns: z, Fx, Fy, Fz, Mx, My, Tz
    z, Fx, Fy, Fz, Mx, My, Tz = A.T
    res_F = np.array([Fx.sum(), Fy.sum(), Fz.sum()])
    # moment about the origin of (Mx, My, Tz) + r x F with r = (0, 0, z):
    # r x F = (-z Fy, z Fx, 0)
    res_M = np.array([np.sum(Mx - z * Fy), np.sum(My + z * Fx), Tz.sum()])
    scale_F = np.abs(A[:, 1:4]).sum()
    scale_M = np.abs(A[:, 4:7]).sum() + np.sum(np.abs(z) * (np.abs(Fx) + np.abs(Fy)))
    assert np.all(np.abs(res_F) <= 1e-9 * scale_F)
    assert np.all(np.abs(res_M) <= 1e-9 * scale_M)


def test_influence_matrix_superposition_matches_solve():
    s = _ss_shaft()
    # columns: load z=0.15 with components Fx, Fy, Fz, Mx, My, Tz (unit each)
    # rows: reactions (4 per bearing: Fx, Fy, Fz, Tz) then internal at z=0.1
    Gm = influence_matrix(s, load_z=[0.15], eval_z=[0.1])
    col_Fy = Gm[:, 1]
    r = s.solve([(0.15, 0.0, -1000.0, 0.0, 0.0, 0.0, 0.0)])
    assert col_Fy[1] * (-1000.0) == pytest.approx(r.reactions[0, 1], rel=1e-9)
    assert col_Fy[5] * (-1000.0) == pytest.approx(r.reactions[1, 1], rel=1e-9)


def test_first_bending_frequency_ss_beam_sanity():
    # Closed form for a simply supported uniform beam: omega1 = pi^2 sqrt(EI/(rho A L^4)).
    # Lumped-mass model with 20 elements (not a spec row; sanity check, 1% tolerance).
    s = Shaft(np.linspace(0.0, 0.3, 21), np.full(20, D), E, NU)
    s.add_bearing(0.0, k_r=RIGID)
    s.add_bearing(0.3, k_r=RIGID)
    rho = 7850.0
    mpl = rho * math.pi * D**2 / 4.0
    I = math.pi * D**4 / 64.0
    f_exact = math.pi**2 * math.sqrt(E * I / (mpl * 0.3**4)) / (2.0 * math.pi)
    f = first_bending_frequency(s, mpl)
    assert f == pytest.approx(f_exact, rel=1e-2)
