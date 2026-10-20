"""V-belt contact tests: NUMERIC TEST TARGETS T4.1-T4.9 (docs/theory/04_contact.md).

Rig solves use Newton on the static equations (belt nodes at rest, v = 0), so the
sliding friction is evaluated at the sheave surface speed (slip = -u). The
anchor (pinned end) reaction is the force the belt receives from the anchor.
"""
import numpy as np
import pytest

from cvtsim.belt.belt import belt_forces, node_mass
from cvtsim.contact.sheave import (
    axial_force_incipient,
    capstan_ratio,
    sheave_contact,
)

EA = 2.435e5
BETA = np.pi / 4
MU = 0.3
R_STAR = 0.05
V_EPS = 1e-3
R_OUT = R_STAR + 5e-3


# ---------------------------------------------------------------- closed forms
def test_T4_1_capstan_ratio():
    """T4.1: capstan_ratio(0.3, pi, pi/4) = 3.79189 (abs tol 1e-5; exact 3.7918912)."""
    assert capstan_ratio(0.3, np.pi, np.pi / 4) == pytest.approx(3.79189, abs=1e-5)


def test_T4_4_axial_force_incipient():
    """T4.4: axial_force_incipient(10, 0.05, pi/4, 0.3) = 235.702 N, abs tol 1e-3."""
    assert axial_force_incipient(10.0, 0.05, np.pi / 4, 0.3) == pytest.approx(235.702, abs=1e-3)


# ---------------------------------------------------------------- rig helpers
def _capstan_path(N):
    """Belt path from end A (right leg bottom) up the right leg, over the upper
    half of the pulley (radius R_STAR, centre origin), down the left leg to end B.
    Legs are 0.1 m. Nodes equally spaced by arc length. Wrap nodes are placed
    1e-6 m inside the free radius (small initial penetration).

    Geometry choice: the spec fixes that A is the tight side (T1 = e^{mu pi/sin b} T2)
    but not which leg A is on. With the pulley turning +Z the friction on the belt
    at the top points toward -X, so the tight side is the +X leg; A is put there."""
    leg = 0.1
    arc = np.pi * R_STAR
    Ltot = 2 * leg + arc
    s = np.linspace(0.0, Ltot, N)
    x = np.zeros((N, 2))
    for k, sk in enumerate(s):
        if sk <= leg:  # right leg, going up from A
            x[k] = (R_STAR, -leg + sk)
        elif sk <= leg + arc:  # over the top, counter-clockwise from angle 0 to pi
            th = (sk - leg) / R_STAR
            r = R_STAR - 1e-6
            x[k] = (r * np.cos(th), r * np.sin(th))
        else:  # left leg, going down to B
            x[k] = (-R_STAR, -(sk - leg - arc))
    return x, s[1] - s[0]


def _rig_residual(x, L0, kn, mu, ext):
    zero_v = np.zeros_like(x)
    f_b, _, _ = belt_forces(x, zero_v, L0, EA, 0.0, 0.0, 0.0, False)
    fc = sheave_contact(x, zero_v, 0.0, 0.0, R_STAR, 0.0, BETA, kn, 0.0, mu, V_EPS,
                        R_OUT, 10.0, 10.0)
    return f_b + fc[0] + ext, fc


def _rig_jacobian(x, L0, kn, ext, idx, h=1e-9):
    """Jacobian of the rig residual over free nodes idx, by colored finite
    differences. Row node i depends only on nodes i-1, i, i+1, so within
    {i-1, i, i+1} exactly one node has each residue mod 3: perturbing all nodes of
    one residue class at once recovers every entry with 6 residual pairs."""
    N = x.shape[0]
    pos = {k: n for n, k in enumerate(idx)}
    J = np.zeros((2 * len(idx), 2 * len(idx)))
    for c in (0, 1):
        for g in range(3):
            nodes = [k for k in idx if k % 3 == g]
            xp = x.copy()
            xm = x.copy()
            xp[nodes, c] += h
            xm[nodes, c] -= h
            Rp, _ = _rig_residual(xp, L0, kn, MU, ext)
            Rm, _ = _rig_residual(xm, L0, kn, MU, ext)
            D = (Rp - Rm) / (2 * h)
            for i in idx:
                for k in (i - 1, i, i + 1):
                    if 1 <= k <= N - 1 and k % 3 == g:
                        for d in (0, 1):
                            J[2 * pos[i] + d, 2 * pos[k] + c] = D[i, d]
    return J


def _solve_capstan_rig(kn, N=150, T2=100.0, tol=1e-9, maxit=60):
    """Newton solve for statics of the capstan rig. Returns (x, F_anchor_on_belt,
    contact outputs at solution). Node 0 = end A (pinned), node N-1 = end B
    (constant force T2 downward)."""
    x, ds = _capstan_path(N)
    eps0 = 1.5e-3  # initial tension guess ~ 365 N, keeps the chain stiff
    L0 = ds / (1.0 + eps0)
    ext = np.zeros_like(x)
    ext[-1] = (0.0, -T2)
    idx = np.arange(1, N)  # free nodes
    for _ in range(maxit):
        R, _ = _rig_residual(x, L0, kn, MU, ext)
        res = np.max(np.abs(R[idx]))
        if res < tol:
            break
        J = _rig_jacobian(x, L0, kn, ext, idx)
        du = np.linalg.solve(J, -R[idx].ravel()).reshape(-1, 2)
        # backtracking line search on the residual norm
        alpha = 1.0
        r0 = np.linalg.norm(R[idx])
        for _ls in range(30):
            xt = x.copy()
            xt[idx] += alpha * du
            Rt, _ = _rig_residual(xt, L0, kn, MU, ext)
            if np.linalg.norm(Rt[idx]) < (1.0 - 1e-4 * alpha) * r0 or alpha < 1e-6:
                break
            alpha *= 0.5
        x = xt
    R, fc = _rig_residual(x, L0, kn, MU, ext)
    assert np.max(np.abs(R[idx])) < 1e-6, "rig Newton did not converge"
    # Pinned end A carries no external load, so the anchor reaction on the belt is
    # minus the total belt+contact force on node 0 (its residual is what the anchor takes).
    f_b, _, _ = belt_forces(x, np.zeros_like(x), L0, EA, 0.0, 0.0, 0.0, False)
    anchor_reaction = -(f_b[0] + fc[0][0])
    return x, anchor_reaction, fc, N


def _rig_ratio_and_thrust(kn):
    x, anchor, fc, N = _solve_capstan_rig(kn)
    T1 = np.linalg.norm(anchor)
    T2 = 100.0
    Fa = fc[1]
    wrap_nodes = int(np.count_nonzero(fc[7] > 0.0))
    return T1 / T2, Fa, T1, T2, wrap_nodes, float(np.max(fc[6])), anchor


# ---------------------------------------------------------------- T4.2 / T4.3
def test_T4_2_capstan_rig_anchor_ratio():
    """T4.2: anchor force / T2 vs 3.79189, rel 2e-2, k_n = 1e7, >= 60 nodes in wrap."""
    ratio, Fa, T1, T2, wrap, dmax, anchor = _rig_ratio_and_thrust(1e7)
    assert wrap >= 60
    assert ratio == pytest.approx(3.79189, rel=2e-2)


def test_T4_3_capstan_rig_axial_thrust():
    """T4.3: F_a = cos(beta)(T1 - T2)/(2 mu) with measured T1, T2, rel 2e-2."""
    ratio, Fa, T1, T2, wrap, dmax, anchor = _rig_ratio_and_thrust(1e7)
    Fa_C5 = np.cos(BETA) * (T1 - T2) / (2 * MU)
    assert Fa == pytest.approx(Fa_C5, rel=2e-2)


# ---------------------------------------------------------------- T4.5
def test_T4_5_penalty_convergence():
    """T4.5: max penetration ~ 1/k_n (successive ratio in [5, 20]); error of T4.2
    non-increasing as k_n rises over {1e5, 1e6, 1e7}."""
    kns = [1e5, 1e6, 1e7]
    pen = []
    err = []
    for kn in kns:
        ratio, Fa, T1, T2, wrap, dmax, anchor = _rig_ratio_and_thrust(kn)
        pen.append(dmax)
        err.append(abs(ratio - 3.79189))
    r1 = pen[0] / pen[1]
    r2 = pen[1] / pen[2]
    assert 5.0 <= r1 <= 20.0
    assert 5.0 <= r2 <= 20.0
    assert err[0] >= err[1] >= err[2]


# ---------------------------------------------------------------- T4.6
def _contact_state(seed, n=12):
    """Random nodes inside the sheave (penetrating by 1e-5..1e-4 m) and outside."""
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0.2, np.pi - 0.2, n)
    depth = rng.uniform(1e-5, 1e-4, n)
    r = R_STAR - depth
    x = np.stack([r * np.cos(ang), r * np.sin(ang)], axis=1)
    v = rng.normal(0.0, 0.2, (n, 2))
    return x, v


def test_T4_6_energy_consistency():
    """T4.6 (mu = 0, c_n = 0): central FD of U_pen w.r.t. node positions equals -f,
    and (dU/dr*)/(2 tan beta) = F_a, both rel 1e-6."""
    kn = 1e6
    x, v = _contact_state(7)
    cx = cy = 0.0
    f, Fa, _, _, _, U, _, _ = sheave_contact(x, v, cx, cy, R_STAR, 0.0, BETA, kn, 0.0,
                                              0.0, V_EPS, R_OUT, 0.0, 0.0)
    h = 1e-8
    grad = np.zeros_like(x)
    for k in range(x.shape[0]):
        for c in (0, 1):
            xp = x.copy()
            xp[k, c] += h
            xm = x.copy()
            xm[k, c] -= h
            Up = sheave_contact(xp, v, cx, cy, R_STAR, 0.0, BETA, kn, 0.0, 0.0, V_EPS,
                                R_OUT, 0.0, 0.0)[5]
            Um = sheave_contact(xm, v, cx, cy, R_STAR, 0.0, BETA, kn, 0.0, 0.0, V_EPS,
                                R_OUT, 0.0, 0.0)[5]
            grad[k, c] = (Up - Um) / (2 * h)
    assert np.allclose(-grad, f, rtol=1e-6, atol=1e-9)

    hr = 1e-8
    Up = sheave_contact(x, v, cx, cy, R_STAR + hr, 0.0, BETA, kn, 0.0, 0.0, V_EPS,
                        R_OUT + hr, 0.0, 0.0)[5]
    Um = sheave_contact(x, v, cx, cy, R_STAR - hr, 0.0, BETA, kn, 0.0, 0.0, V_EPS,
                        R_OUT - hr, 0.0, 0.0)[5]
    dU_drs = (Up - Um) / (2 * hr)
    assert dU_drs / (2 * np.tan(BETA)) == pytest.approx(Fa, rel=1e-6)


# ---------------------------------------------------------------- T4.7
def test_T4_7_friction_dissipates():
    """T4.7: P_diss >= 0 for random states (friction is dissipative)."""
    rng = np.random.default_rng(3)
    for trial in range(200):
        n = 20
        ang = rng.uniform(0.0, 2 * np.pi, n)
        r = R_STAR - rng.uniform(-2e-3, 3e-3, n)
        x = np.stack([r * np.cos(ang), r * np.sin(ang)], axis=1)
        v = rng.normal(0.0, 3.0, (n, 2))
        out = sheave_contact(x, v, 0.0, 0.0, R_STAR, rng.normal(0, 0.01), BETA,
                             1e6, 10.0, MU, V_EPS, R_OUT,
                             rng.normal(0, 50.0), rng.normal(0, 50.0))
        assert out[4] >= 0.0


# ---------------------------------------------------------------- T4.8
def test_T4_8_no_contact_outside():
    """T4.8: forces are exactly zero for r > r* and for r > r_outer."""
    x = np.array([
        [R_STAR + 1e-3, 0.0],          # beyond free radius
        [0.0, R_STAR + 2e-2],          # far outside
        [0.0, -(R_STAR + 1e-9)],       # just outside
    ])
    v = np.array([[1.0, 2.0], [0.5, -0.5], [0.0, 3.0]])
    out = sheave_contact(x, v, 0.0, 0.0, R_STAR, 0.0, BETA, 1e7, 100.0, MU, V_EPS,
                         R_OUT, 10.0, 10.0)
    f, Fa, tqf, tqm, Pd, U, delta, Nf = out
    assert np.array_equal(f, np.zeros_like(f))
    assert Fa == 0.0 and tqf == 0.0 and tqm == 0.0 and Pd == 0.0 and U == 0.0
    assert np.all(delta == 0.0) and np.all(Nf == 0.0)
    # r_outer below r*: node inside r* but beyond r_outer must also see no force
    x2 = np.array([[0.0, 0.0499], [0.0, 0.0495]])  # r = 0.0499, 0.0495
    out2 = sheave_contact(x2, np.zeros_like(x2), 0.0, 0.0, R_STAR, 0.0, BETA, 1e7, 0.0,
                          MU, V_EPS, 0.0495 - 1e-6, 10.0, 10.0)
    assert np.array_equal(out2[0], np.zeros_like(x2))


# ---------------------------------------------------------------- T4.9 (slow)
def _stadium_path(N, C, R):
    """Closed belt path for two equal pulleys at (0,0) and (C,0), CCW order."""
    L = 2 * C + 2 * np.pi * R
    s = np.arange(N) * (L / N)
    x = np.zeros((N, 2))
    a1 = C
    a2 = np.pi * R
    a3 = C
    for k, sk in enumerate(s):
        if sk < a1:  # bottom, +X
            x[k] = (sk, -R)
        elif sk < a1 + a2:  # right arc, CCW from -pi/2 to pi/2 about (C,0)
            th = -np.pi / 2 + (sk - a1) / R
            x[k] = (C + R * np.cos(th), R * np.sin(th))
        elif sk < a1 + a2 + a3:  # top, -X
            x[k] = (C - (sk - a1 - a2), R)
        else:  # left arc, CCW from pi/2 to 3pi/2 about (0,0)
            th = np.pi / 2 + (sk - a1 - a2 - a3) / R
            x[k] = (R * np.cos(th), R * np.sin(th))
    return x, L / N


def _driven_loss(Tq, t_end=0.3, h=1e-6):
    """Two-pulley closed-belt drive: driver at constant omega = 100 rad/s, driven
    pulley (I = 1e-3) with load torque Tq. Returns relative speed loss at steady state.

    Choices not fixed by the spec (see report): rhoA = 0.345 kg/m (T3.2 value);
    v_eps = 0.05 m/s (contact doc placeholder); eta = 0. With v_eps = 1e-3 the
    friction velocity-stiffness mu*N/v_eps exceeds the explicit Verlet limit at
    h = 1e-6, and with eta = 1e-4 the Kelvin-Voigt rate EA*eta/(L0*m) ~ 2e6 /s does.
    The steady sliding state does not depend on axial damping."""
    N, C, R = 120, 0.2, 0.05
    kn = 1e7
    rhoA = 0.345
    eta = 0.0
    v_eps = 0.05
    I_d = 1e-3
    w_in = 100.0
    x, ds = _stadium_path(N, C, R)
    L0 = 0.998 * ds  # rest length 0.2 % shorter than the path -> ~490 N pretension
    m = node_mass(rhoA, L0)
    cn = 0.1 * np.sqrt(kn * m)  # 5 % of critical normal damping
    v = np.zeros_like(x)
    w = 0.0

    def acc(xx, vv, ww):
        fb, _, _ = belt_forces(xx, vv, L0, EA, eta, 0.0, 0.0, True)
        f1 = sheave_contact(xx, vv, 0.0, 0.0, R, 0.0, BETA, kn, cn, MU, v_eps,
                            R + 5e-3, w_in, w_in)
        f2 = sheave_contact(xx, vv, C, 0.0, R, 0.0, BETA, kn, cn, MU, v_eps,
                            R + 5e-3, ww, ww)
        F = fb + f1[0] + f2[0]
        a = F / m
        alpha = (f2[2] + f2[3] - Tq) / I_d
        return a, alpha

    a, al = acc(x, v, w)
    n_steps = int(round(t_end / h))
    hist = []
    for step in range(n_steps):
        vm = v + 0.5 * h * a
        wm = w + 0.5 * h * al
        x = x + h * vm
        a, al = acc(x, vm, wm)
        v = vm + 0.5 * h * a
        w = wm + 0.5 * h * al
        if step >= int(0.75 * n_steps):
            hist.append(w)
    hist = np.asarray(hist)
    half = len(hist) // 2
    m1, m2 = hist[:half].mean(), hist[half:].mean()
    w_ss = float(hist.mean())
    # steadiness: the mean over the two halves of the averaging window agrees
    assert abs(m2 - m1) / abs(w_ss) < 1e-3, f"driven pulley not steady ({m1:.4f} vs {m2:.4f})"
    return 1.0 - w_ss / w_in


@pytest.mark.slow
def test_T4_9_slip_increases_with_load():
    """T4.9 (Gerbert, qualitative): relative speed loss 1 - w_out/w_in is strictly
    increasing in load torque Tq in {2, 5, 10, 15} N m."""
    losses = [_driven_loss(Tq) for Tq in (2.0, 5.0, 10.0, 15.0)]
    assert all(b > a for a, b in zip(losses, losses[1:])), losses
