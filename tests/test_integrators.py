"""Numeric tests for cvtsim.core.integrators (docs/theory/02_integrators.md, sec. 8: T2.1-T2.10)."""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from cvtsim.core.integrators import (
    GeneralizedAlpha,
    genalpha_amplification,
    genalpha_params,
    symplectic_euler,
    velocity_verlet,
)

OMEGA = 2.0 * np.pi


def _genalpha_sdof(rho, h, t_end=1.0):
    k = OMEGA ** 2
    ga = GeneralizedAlpha(M=np.eye(1), f_int=lambda x, v, t: k * x, f_ext=None,
                          rho_inf=rho, K=np.array([[k]]), Cd=np.zeros((1, 1)))
    nsteps = int(round(t_end / h))
    t, X, _, _ = ga.run([1.0], [0.0], 0.0, h, nsteps)
    return np.max(np.abs(X[:, 0] - np.cos(OMEGA * t)))


def _verlet_sdof(h, t_end=1.0):
    k = OMEGA ** 2
    t, X, _ = velocity_verlet(lambda x, v, t: -k * x, [1.0], [0.0], 0.0, h, int(round(t_end / h)))
    return np.max(np.abs(X[:, 0] - np.cos(OMEGA * t)))


def _euler_sdof(h, t_end=1.0):
    k = OMEGA ** 2
    t, X, _ = symplectic_euler(lambda x, v, t: -k * x, [1.0], [0.0], 0.0, h, int(round(t_end / h)))
    return np.max(np.abs(X[:, 0] - np.cos(OMEGA * t)))


# ---- T2.1 ------------------------------------------------------------------

def test_T2_1_genalpha_params():
    am, af, gam, bet = genalpha_params(0.8)
    np.testing.assert_allclose([am, af, gam, bet],
                               [0.333333, 0.444444, 0.611111, 0.308642], atol=1e-6)


# ---- T2.2 : gen-alpha order ---------------------------------------------------

def test_T2_2_genalpha_order_two():
    e = {}
    for n in (100, 200, 400):
        e[n] = _genalpha_sdof(0.8, 1.0 / n)
    o1 = np.log2(e[100] / e[200])
    o2 = np.log2(e[200] / e[400])
    assert abs(o1 - 2.0) <= 0.15, (o1, e)
    assert abs(o2 - 2.0) <= 0.15, (o2, e)


# ---- T2.3 - T2.5 : amplification matrix --------------------------------------

def test_T2_3_spectral_radius_at_high_frequency():
    for rho in (0.0, 0.5, 0.8, 1.0):
        G = genalpha_amplification(1e5, rho)
        sr = np.max(np.abs(np.linalg.eigvals(G)))
        assert abs(sr - rho) <= 2e-3, (rho, sr)


def test_T2_4_unconditional_stability():
    for rho in (0.0, 0.5, 1.0):
        for z in np.logspace(-3, 4, 200):
            sr = np.max(np.abs(np.linalg.eigvals(genalpha_amplification(z, rho))))
            assert sr <= 1.0 + 1e-12, (rho, z, sr)


def test_T2_5_no_dissipation_at_rho_one():
    for z in np.logspace(-3, 4, 200):
        sr = np.max(np.abs(np.linalg.eigvals(genalpha_amplification(z, 1.0))))
        assert abs(sr - 1.0) <= 1e-10, (z, sr)


# ---- T2.6 - T2.7 : Verlet ------------------------------------------------------

def test_T2_6_verlet_order_two():
    e = {}
    for n in (100, 200, 400):
        e[n] = _verlet_sdof(1.0 / n)
    o1 = np.log2(e[100] / e[200])
    o2 = np.log2(e[200] / e[400])
    assert abs(o1 - 2.0) <= 0.15, (o1, e)
    assert abs(o2 - 2.0) <= 0.15, (o2, e)


def test_T2_7_verlet_stability_threshold():
    def max_abs(wh):
        h = wh / 1.0  # omega = 1
        _, X, _ = velocity_verlet(lambda x, v, t: -x, [1.0], [0.0], 0.0, h, 10000)
        # treat overflow/NaN after blow-up as unbounded
        return float(np.max(np.nan_to_num(np.abs(X[:, 0]), nan=np.inf, posinf=np.inf)))

    assert max_abs(1.99) <= 20.0
    assert max_abs(2.01) > 1e3


# ---- T2.8 : Verlet energy, pendulum ----------------------------------------

def test_T2_8_verlet_pendulum_energy():
    g_over_L = 9.81
    h = 0.01
    nsteps = 100000

    def accel(x, v, t):
        return -g_over_L * np.sin(x)

    _, X, V = velocity_verlet(accel, [1.0], [0.0], 0.0, h, nsteps)
    th = X[:, 0]
    w = V[:, 0]
    E = 0.5 * w ** 2 - g_over_L * np.cos(th)
    E0 = E[0]
    err = np.abs(E - E0) / np.abs(E0)
    assert err.max() < 1e-3, err.max()
    n10 = (nsteps + 1) // 10
    assert err[-n10:].max() <= 1.5 * err[:n10].max(), (err[-n10:].max(), err[:n10].max())


# ---- T2.9 : nonlinear gen-alpha vs reference -------------------------------

def test_T2_9_duffing_genalpha_matches_reference():
    ga = GeneralizedAlpha(M=np.eye(1), f_int=lambda x, v, t: x + x ** 3, f_ext=None, rho_inf=0.9)
    _, X, _, _ = ga.run([1.0], [0.0], 0.0, 1e-3, 5000)
    x_ga = X[-1, 0]

    def rhs(t, y):
        return [y[1], -y[0] - y[0] ** 3]

    sol = solve_ivp(rhs, (0.0, 5.0), [1.0, 0.0], method='RK45', rtol=1e-10, atol=1e-12)
    x_ref = sol.y[0, -1]
    assert abs(x_ga - x_ref) < 1e-4, (x_ga, x_ref)


# ---- T2.10 : symplectic Euler order ----------------------------------------

def test_T2_10_symplectic_euler_order_one():
    e = {}
    for n in (100, 200, 400):
        e[n] = _euler_sdof(1.0 / n)
    o1 = np.log2(e[100] / e[200])
    o2 = np.log2(e[200] / e[400])
    assert abs(o1 - 1.0) <= 0.15, (o1, e)
    assert abs(o2 - 1.0) <= 0.15, (o2, e)
