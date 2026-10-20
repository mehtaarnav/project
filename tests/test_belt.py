"""Segmented belt tests: NUMERIC TEST TARGETS T3.1-T3.8 (docs/theory/03_belt.md).

Discretization and pre-stretch convention used throughout: nodes are laid out at
the stretched spacing l = span/N and the rest length is L0 = l/(1 + eps0), with
eps0 = T0/EA. Lumped mass is m = rhoA * L0 (rest length, as in belt_forces).
"""
import numpy as np
import pytest

from cvtsim.belt.belt import (
    belt_forces,
    belt_length,
    belt_tensions,
    init_belt_on_path,
    node_mass,
)

EA = 243.5e3
RHOA = 0.345  # kg/m (task spec, T3.2)


# ---------------------------------------------------------------- helpers
def _straight_chain(span, N, T0, EA_=EA):
    """Open chain of N segments along x, stretched spacing span/N, pre-tension T0."""
    l = span / N
    eps0 = T0 / EA_
    L0 = l / (1.0 + eps0)
    x = np.zeros((N + 1, 2))
    x[:, 0] = np.arange(N + 1) * l
    return x, L0


def _newton_static(x, free, ext, EA_, L0, tol=1e-10, maxit=50):
    """Solve f_int(x) + ext = 0 on free nodes (open chain, EI = 0, eta = 0)."""
    x = x.copy()
    zero_v = np.zeros_like(x)
    idx = np.flatnonzero(free)
    for _ in range(maxit):
        f, _, _ = belt_forces(x, zero_v, L0, EA_, 0.0, 0.0, 0.0, False)
        R = (f + ext)[idx].ravel()
        if np.max(np.abs(R)) < tol:
            break
        h = 1e-7
        J = np.zeros((R.size, R.size))
        for j, (k, c) in enumerate([(k, c) for k in idx for c in (0, 1)]):
            xp = x.copy()
            xp[k, c] += h
            xm = x.copy()
            xm[k, c] -= h
            fp, _, _ = belt_forces(xp, zero_v, L0, EA_, 0.0, 0.0, 0.0, False)
            fm, _, _ = belt_forces(xm, zero_v, L0, EA_, 0.0, 0.0, 0.0, False)
            J[:, j] = ((fp - fm) / (2 * h))[idx].ravel()
        dx = np.linalg.solve(J, -R)
        x[idx] += dx.reshape(-1, 2)
    return x


def _transverse_jacobian(x, L0, EA_, EI, N_nodes_free_idx, h=1e-7):
    """Numerical d(f_y)/d(y) on the free interior nodes (EI given, eta = 0)."""
    zero_v = np.zeros_like(x)
    n = len(N_nodes_free_idx)
    K = np.zeros((n, n))
    for j, k in enumerate(N_nodes_free_idx):
        xp = x.copy()
        xp[k, 1] += h
        xm = x.copy()
        xm[k, 1] -= h
        fp, _, _ = belt_forces(xp, zero_v, L0, EA_, 0.0, EI, 0.0, False)
        fm, _, _ = belt_forces(xm, zero_v, L0, EA_, 0.0, EI, 0.0, False)
        col = -(fp[N_nodes_free_idx, 1] - fm[N_nodes_free_idx, 1]) / (2 * h)
        K[:, j] = col
    return K


def _lowest_transverse_hz(span, N, T0, EI, rhoA=RHOA):
    x, L0 = _straight_chain(span, N, T0)
    m = node_mass(rhoA, L0)
    free = list(range(1, N))
    K = _transverse_jacobian(x, L0, EA, EI, free)
    K = 0.5 * (K + K.T)
    lam = np.linalg.eigvalsh(K / m)
    return np.sqrt(lam[0]) / (2 * np.pi)


# ---------------------------------------------------------------- T3.1
def test_T3_1_static_stretch():
    """Open chain L=1 m, N=50 nodes, P=500 N at the free end, node 0 fixed.
    Target delta = P L / EA (2.0533881e-3 m), rel 1e-9."""
    N_nodes = 50
    L = 1.0
    P = 500.0
    L0 = L / (N_nodes - 1)
    x = np.zeros((N_nodes, 2))
    x[:, 0] = np.arange(N_nodes) * L0
    ext = np.zeros_like(x)
    ext[-1, 0] = P
    free = np.ones(N_nodes, dtype=bool)
    free[0] = False
    xs = _newton_static(x, free, ext, EA, L0)
    delta = xs[-1, 0] - x[-1, 0]
    target = P * L / EA
    assert delta == pytest.approx(target, rel=1e-9)
    # tension in every element equals P
    T = belt_tensions(xs, np.zeros_like(xs), L0, EA, 0.0, False)
    assert np.allclose(T, P, rtol=1e-9)


# ---------------------------------------------------------------- T3.2
def test_T3_2_string_discrete():
    """Span 0.2 m, N=20 segments, T0=100 N, rhoA=0.345 kg/m, EI=0.
    Lowest transverse eigenfrequency of the linearized system vs B2:
    42.5278 Hz (B2 with pre-stretch), rel 1e-4."""
    f1 = _lowest_transverse_hz(0.2, 20, 100.0, 0.0)
    assert f1 == pytest.approx(42.5278, rel=1e-4)


def test_T3_3_string_continuum():
    """Same setup, N=20 vs continuum 42.5628 Hz, rel 2e-3."""
    f1 = _lowest_transverse_hz(0.2, 20, 100.0, 0.0)
    assert f1 == pytest.approx(42.5628, rel=2e-3)


# ---------------------------------------------------------------- T3.4
def test_T3_4_string_with_bending():
    """N=100, EI=0.05, span 0.2 m, T0=100 N, compared with
    sqrt((T0/rhoA) k^2 + (EI/rhoA) k^4)/(2 pi), k = pi/span. rel 2e-2."""
    span, T0, EI = 0.2, 100.0, 0.05
    k = np.pi / span
    target = np.sqrt((T0 / RHOA) * k**2 + (EI / RHOA) * k**4) / (2 * np.pi)
    f1 = _lowest_transverse_hz(span, 100, T0, EI)
    assert f1 == pytest.approx(target, rel=2e-2)


# ---------------------------------------------------------------- T3.5
def test_T3_5_time_domain_string():
    """Velocity Verlet (h = 1e-6 s) of the T3.2 string, initial shape = discrete
    mode 1 with amplitude 1e-4 m. Measured period from zero crossings of the
    midspan node vs 1/42.5278 s, rel 1e-3."""
    N = 20
    span, T0 = 0.2, 100.0
    x0, L0 = _straight_chain(span, N, T0)
    m = node_mass(RHOA, L0)
    j = np.arange(N + 1)
    mode = np.sin(np.pi * j / N)
    mode[0] = mode[-1] = 0.0
    x = x0.copy()
    x[:, 1] = 1e-4 * mode
    v = np.zeros_like(x)
    free = np.ones(N + 1, dtype=bool)
    free[0] = free[-1] = False
    h = 1e-6
    n_steps = 120_000  # ~5 periods (T ~ 23.5 ms)
    mid = N // 2

    def accel(xx, vv):
        f, _, _ = belt_forces(xx, vv, L0, EA, 0.0, 0.0, 0.0, False)
        a = f / m
        a[~free] = 0.0
        return a

    a = accel(x, v)
    ts = []
    ys = []
    t = 0.0
    for _ in range(n_steps):
        v = v + 0.5 * h * a
        x = x + h * v
        a = accel(x, v)
        v = v + 0.5 * h * a
        ys.append(x[mid, 1])
        t += h
        ts.append(t)
    ys = np.asarray(ys)
    ts = np.asarray(ts)
    # zero crossings with linear interpolation
    s = np.sign(ys)
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    tc = ts[idx] - ys[idx] * (ts[idx + 1] - ts[idx]) / (ys[idx + 1] - ys[idx])
    # successive crossings are half a period apart
    half = np.diff(tc)
    period = 2.0 * np.mean(half)
    assert period == pytest.approx(1.0 / 42.5278, rel=1e-3)


# ---------------------------------------------------------------- T3.6
def test_T3_6_energy_closed_belt():
    """Closed belt, EI > 0, eta = 0, cb = 0, small random velocities,
    velocity Verlet 1e4 steps: max|dE|/E0 < 1e-3."""
    N = 100
    R = 0.05
    th = 2 * np.pi * np.arange(N) / N
    x = np.stack([R * np.cos(th), R * np.sin(th)], axis=1)
    L0 = 2 * R * np.sin(np.pi / N)  # zero axial strain
    m = node_mass(RHOA, L0)
    EI = 0.05
    rng = np.random.default_rng(12345)
    v = 1e-3 * rng.standard_normal(x.shape)
    h = 1e-6

    def energy(xx, vv):
        _, U, _ = belt_forces(xx, vv, L0, EA, 0.0, EI, 0.0, True)
        return 0.5 * m * np.sum(vv * vv) + U

    def accel(xx, vv):
        f, _, _ = belt_forces(xx, vv, L0, EA, 0.0, EI, 0.0, True)
        return f / m

    E0 = energy(x, v)
    a = accel(x, v)
    max_dE = 0.0
    for _ in range(10_000):
        v = v + 0.5 * h * a
        x = x + h * v
        a = accel(x, v)
        v = v + 0.5 * h * a
        max_dE = max(max_dE, abs(energy(x, v) - E0))
    assert max_dE / E0 < 1e-3


# ---------------------------------------------------------------- T3.7
def test_T3_7_rigid_rotation_no_damping_force():
    """Closed belt on a circle rotating rigidly at omega, eta > 0: the axial
    damping force is zero. Checked as f(eta) - f(0) at identical state, which
    isolates the damping contribution; target |f_damp| < 1e-9. (Spec setup has
    no bending damper; see the cb test below.)"""
    N = 120
    R = 0.05
    th = 2 * np.pi * np.arange(N) / N
    x = np.stack([R * np.cos(th), R * np.sin(th)], axis=1)
    L0 = 2 * R * np.sin(np.pi / N)
    omega = 37.0
    v = omega * np.stack([-x[:, 1], x[:, 0]], axis=1)  # zhat x x
    f_d, _, P_d = belt_forces(x, v, L0, EA, 1e-4, 0.05, 0.0, True)
    f_0, _, P_0 = belt_forces(x, v, L0, EA, 0.0, 0.05, 0.0, True)
    assert np.max(np.abs(f_d - f_0)) < 1e-9
    assert P_d < 1e-12  # no dissipation for rigid rotation


def test_T3_7_supplement_bending_damper_rigid_rotation():
    """Supplement to T3.7: with the bending damper cb > 0, the rigid-rotation
    damping force is also zero (phi is rotation invariant). Tolerance 1e-6 N
    here, because cb*|grad phi|^2 amplifies round-off of v = omega x r by ~1e5
    for this geometry, so a 1e-9 absolute bound is below the floating-point floor."""
    N = 120
    R = 0.05
    th = 2 * np.pi * np.arange(N) / N
    x = np.stack([R * np.cos(th), R * np.sin(th)], axis=1)
    L0 = 2 * R * np.sin(np.pi / N)
    omega = 37.0
    v = omega * np.stack([-x[:, 1], x[:, 0]], axis=1)
    f_d, _, _ = belt_forces(x, v, L0, EA, 1e-4, 0.05, 0.2, True)
    f_0, _, _ = belt_forces(x, v, L0, EA, 0.0, 0.05, 0.0, True)
    f_c, _, _ = belt_forces(x, v, L0, EA, 0.0, 0.05, 0.2, True)
    assert np.max(np.abs(f_d - f_0)) < 1e-6
    assert np.max(np.abs(f_c - belt_forces(x, v, L0, EA, 0.0, 0.05, 0.0, True)[0])) < 1e-6


# ---------------------------------------------------------------- T3.8
def test_T3_8_path_length():
    """init_belt_on_path(0.0101964, 0.0305891, 0.251495, 200): total polyline
    length vs belt_length (chord error), rel 1e-3."""
    rp, rs, C, N = 0.0101964, 0.0305891, 0.251495, 200
    x, L0 = init_belt_on_path(rp, rs, C, N)
    chords = np.linalg.norm(np.roll(x, -1, axis=0) - x, axis=1)  # closed loop
    exact = belt_length(rp, rs, C)
    assert chords.sum() == pytest.approx(exact, rel=1e-3)
    assert L0 * N == pytest.approx(exact, rel=1e-12)
    assert chords.min() > 0.99 * L0  # no coincident nodes
