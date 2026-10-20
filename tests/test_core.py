"""Numeric tests for cvtsim.core (docs/theory/01_rigid_body.md, sec. 7: T1.1-T1.13).

Targets and tolerances are taken verbatim from the spec.

Pendulum setup: the revolute axis is +Z, so gravity must lie in the swing plane;
the pendulum tests use g = (0, -9.81, 0). The CoM is at (L sin th, -L cos th, 0).
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from cvtsim.core.joints import Prismatic, Revolute, Screw
from cvtsim.core.mbs import System
from cvtsim.core.rigid import RigidBody, quat_exp

G = 9.81
L = 1.0
PENDULUM_GRAVITY = (0.0, -G, 0.0)


def _pendulum_system(theta0_deg, baumgarte=(0.0, 0.0), project=False):
    th0 = np.radians(theta0_deg)
    body = RigidBody(1.0, 1e-6 * np.eye(3),
                     r=[L * np.sin(th0), -L * np.cos(th0), 0.0],
                     q=quat_exp([0.0, 0.0, th0]),   # Rz(theta0): (0, L, 0) -> pivot offset
                     v=[0.0, 0.0, 0.0], w=[0.0, 0.0, 0.0])
    joint = Revolute(None, 0, si=[0.0, 0.0, 0.0], sj=[0.0, L, 0.0],
                     ai=[0.0, 0.0, 1.0], aj=[0.0, 0.0, 1.0])
    return body, System([body], [joint], gravity=PENDULUM_GRAVITY,
                        baumgarte=baumgarte, project=project)


def _period(t, th):
    """Mean period from same-direction (negative -> positive) zero crossings, linearly interpolated."""
    idx = np.where((th[:-1] < 0.0) & (th[1:] >= 0.0))[0]
    tc = t[idx] + (0.0 - th[idx]) * (t[idx + 1] - t[idx]) / (th[idx + 1] - th[idx])
    assert tc.size >= 4, "need at least 3 periods"
    return float(np.mean(np.diff(tc)))


@functools.lru_cache(maxsize=None)
def _period_run(theta0_deg, h, method, t_end, project):
    body, s = _pendulum_system(theta0_deg, project=project)
    nsteps = int(round(t_end / h))
    t = np.arange(nsteps + 1) * h
    th = np.empty(nsteps + 1)
    th[0] = np.radians(theta0_deg)
    for k in range(nsteps):
        s.step(h, method=method)
        th[k + 1] = np.arctan2(body.r[0], -body.r[1])
    return _period(t, th)


@functools.lru_cache(maxsize=None)
def _shared_run(baumgarte_on, project, method='symplectic_euler', h=1e-3, t_end=10.0):
    """Pendulum at 60 deg. Returns (energy error, constraint violation) histories."""
    a = b = 1.0 / (10.0 * h) if baumgarte_on else 0.0
    body, s = _pendulum_system(60.0, baumgarte=(a, b), project=project)
    E0 = s.energy()
    n = int(round(t_end / h))
    dE = np.empty(n)
    cv = np.empty(n)
    for k in range(n):
        s.step(h, method=method)
        dE[k] = abs(s.energy() - E0) / abs(E0)
        cv[k] = s.constraint_violation()
    return dE, cv


# ---- T1.1 - T1.3: pendulum period -------------------------------------------

def test_T1_1_period_10deg():
    # h = 1e-4, symplectic Euler with projection; run to 3.85 periods (>= 3 periods averaged)
    T = 2.00989
    Tm = _period_run(10.0, 1e-4, 'symplectic_euler', 7.7, True)
    assert abs(Tm - T) / T < 1e-3, (Tm, T)


def test_T1_2_period_60deg():
    # RK4 with h = 5e-3 (step size is not fixed by the spec for this row)
    T = 2.15287
    Tm = _period_run(60.0, 5e-3, 'rk4', 8.3, False)
    assert abs(Tm - T) / T < 1e-3, (Tm, T)


def test_T1_3_period_120deg():
    # RK4 with h = 5e-3 (step size is not fixed by the spec for this row)
    T = 2.75409
    Tm = _period_run(120.0, 5e-3, 'rk4', 10.6, False)
    assert abs(Tm - T) / T < 1e-3, (Tm, T)


# ---- T1.4 - T1.7: energy and constraint drift --------------------------------

def test_T1_4_energy_rk4_projection():
    dE, _ = _shared_run(False, True, method='rk4')
    assert dE.max() < 1e-6, dE.max()


def test_T1_4b_symplectic_euler_energy_drift_first_order():
    # Projected symplectic Euler is not symplectic on the constraint manifold: the energy
    # error drifts linearly, at a rate proportional to h (01 spec, T1.4b).
    d1 = _shared_run(False, True, h=1e-3, t_end=2.0)[0][-1]
    d2 = _shared_run(False, True, h=5e-4, t_end=2.0)[0][-1]
    assert 1.7 <= d1 / d2 <= 2.3, (d1, d2)


def test_T1_5_constraint_drift_baumgarte():
    _, cv = _shared_run(True, False)
    h = 1e-3
    bound = 10.0 * G * h * h
    assert cv.max() <= bound, (cv.max(), bound)


def test_T1_6_drift_with_projection():
    _, cv = _shared_run(False, True)
    assert cv.max() <= 1e-10, cv.max()


def test_T1_7_drift_without_stabilization_is_larger():
    _, cv_unstab = _shared_run(False, False)
    _, cv_proj = _shared_run(False, True)
    assert cv_unstab.max() > cv_proj.max(), (cv_unstab.max(), cv_proj.max())


# ---- T1.8 - T1.10: free symmetric top (rk4, h = 1e-3, 10 s) -----------------

@functools.lru_cache(maxsize=None)
def _top():
    body = RigidBody(1.0, np.diag([1.0, 1.0, 2.0]), [0, 0, 0], [1, 0, 0, 0],
                     [0, 0, 0], [1.0, 0.0, 3.0])
    s = System([body], [], gravity=(0, 0, 0))
    h = 1e-3
    n = 10000
    t = np.arange(n + 1) * h
    wb = np.empty((n + 1, 3))
    e3 = np.empty((n + 1, 3))
    Lw = np.empty((n + 1, 3))

    def record(k):
        R = body.R
        wb[k] = R.T @ body.w
        e3[k] = R[:, 2]
        Lw[k] = body.Jw @ body.w

    record(0)
    for k in range(n):
        s.step(h, method='rk4')
        record(k + 1)
    return t, wb, e3, Lw


def test_T1_8_body_frame_precession_rate():
    t, wb, _, _ = _top()
    phi = np.unwrap(np.arctan2(wb[:, 1], wb[:, 0]))
    rate = np.polyfit(t, phi, 1)[0]
    target = (2.0 - 1.0) / 1.0 * 3.0  # (I3 - I1)/I1 * w3 = 3.000
    assert abs(rate - target) / target < 1e-3, (rate, target)


def test_T1_9_space_precession_about_L():
    t, _, e3, Lw = _top()
    Lvec = Lw[0]
    Lhat = Lvec / np.linalg.norm(Lvec)
    b1 = np.cross(Lhat, [1.0, 0.0, 0.0])
    if np.linalg.norm(b1) < 1e-6:
        b1 = np.cross(Lhat, [0.0, 1.0, 0.0])
    b1 /= np.linalg.norm(b1)
    b2 = np.cross(Lhat, b1)
    p = e3 - np.outer(e3 @ Lhat, Lhat)
    phi = np.unwrap(np.arctan2(p @ b2, p @ b1))
    rate = np.polyfit(t, phi, 1)[0]
    target = np.linalg.norm(Lvec) / 1.0  # |L|/I1 = 6.08276
    assert target == pytest.approx(6.08276, rel=1e-5)
    # rate is the signed precession about +Lhat; the spec row is the magnitude
    assert abs(abs(rate) - target) / target < 1e-3, (rate, target)


def test_T1_10_angular_momentum_conserved():
    _, _, _, Lw = _top()
    Lmag = np.linalg.norm(Lw, axis=1)
    drift = np.max(np.abs(Lmag - Lmag[0])) / Lmag[0]
    assert drift < 1e-6, drift


# ---- T1.11 - T1.12: screw joint ---------------------------------------------

def _screw_system():
    body = RigidBody(2.0, 0.01 * np.eye(3), [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0, 0, 0])
    joint = Screw(None, 0, [0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 1], 0.01)
    s = System([body], [joint], gravity=(0, 0, 0), project=True)
    return body, joint, s


def test_T1_11_screw_axial_acceleration():
    body, _, s = _screw_system()
    body.force = np.array([0.0, 0.0, 10.0])
    acc = s.accelerations()
    target = 10.0 / (2.0 + 0.01 / 0.01 ** 2)  # F/(m + I/p^2) = 0.0980392
    assert acc[2] == pytest.approx(target, rel=1e-6), (acc[2], target)


def test_T1_12_screw_constraint_holds():
    body, joint, s = _screw_system()
    body.force = np.array([0.0, 0.0, 10.0])
    h = 1e-3
    worst = 0.0
    for _ in range(1000):
        s.step(h, method='symplectic_euler')
        worst = max(worst, abs(body.r[2] - joint.p * joint.theta(s)))
    assert worst <= 1e-9, worst


# ---- T1.13: prismatic block on incline --------------------------------------

def test_T1_13_prismatic_incline_acceleration():
    th = np.radians(30.0)
    t_dir = np.array([np.cos(th), 0.0, -np.sin(th)])   # downhill direction
    body = RigidBody(1.0, np.eye(3), [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0, 0, 0])
    joint = Prismatic(None, 0, [0, 0, 0], [0, 0, 0], t_dir, t_dir)
    s = System([body], [joint], gravity=(0.0, 0.0, -G))
    acc = s.accelerations()
    along = float(acc[0:3] @ t_dir)
    target = G * np.sin(th)  # 4.905
    assert along == pytest.approx(target, rel=1e-9), (along, target)
