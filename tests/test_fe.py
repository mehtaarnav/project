"""FE tests: docs/theory/08_fe_tet10.md section 7 (T8.1-T8.6)."""
import numpy as np
import pytest

from cvtsim.fe import mesh as M
from cvtsim.fe import tet10 as T

# ---------------------------------------------------------------- helpers

def _D_ref(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def _voigt_strain(A):
    """Engineering-strain Voigt vector of the linear field u = A x + c."""
    S = 0.5 * (A + A.T)
    return np.array([S[0, 0], S[1, 1], S[2, 2], 2 * S[1, 2], 2 * S[0, 2], 2 * S[0, 1]])


def kirsch_cart(x, y, a=1.0, s=1.0):
    """Cartesian stress (sxx, syy, sxy) of the Kirsch solution for uniaxial s along x."""
    r2 = x * x + y * y
    th = np.arctan2(y, x)
    c, sn = np.cos(th), np.sin(th)
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    ar2 = a * a / r2
    ar4 = ar2 * ar2
    srr = s / 2 * (1 - ar2) + s / 2 * (1 - 4 * ar2 + 3 * ar4) * c2
    stt = s / 2 * (1 + ar2) - s / 2 * (1 + 3 * ar4) * c2
    srt = -s / 2 * (1 + 2 * ar2 - 3 * ar4) * s2
    sxx = srr * c * c + stt * sn * sn - 2 * srt * sn * c
    syy = srr * sn * sn + stt * c * c + 2 * srt * sn * c
    sxy = (srr - stt) * sn * c + srt * (c * c - sn * sn)
    return sxx, syy, sxy


A_HOLE, W_PLATE, T_PLATE = 1.0, 5.0, 0.2


def _iso_grading(n_r, n_t, a=A_HOLE, span=W_PLATE - A_HOLE):
    """Radial grading ratio q such that the first radial layer at the hole is as thick
    as the hole's angular element arc (a * (pi/2) / (2 n_t)). Isotropic hole elements.
    Solves (q-1)/(q^n_r - 1) * span = arc by bisection."""
    arc = a * (np.pi / 2) / (2 * n_t)
    f = lambda q: (q - 1.0) / (q**n_r - 1.0) * span - arc  # noqa: E731
    lo, hi = 1.0 + 1e-12, 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _kirsch_solve(n_r, n_t, grading=None):
    """Quarter plate with hole, exact Kirsch tractions on x=W and y=W, nu=0, E=1."""
    if grading is None:
        grading = _iso_grading(n_r, n_t)
    m = M.plate_hole_quarter_mesh(A_HOLE, W_PLATE, T_PLATE, n_r, n_t, 1, grading=grading)
    faces = np.vstack([
        M.boundary_faces(m, lambda p: np.isclose(p[:, 0], W_PLATE, atol=1e-12)),
        M.boundary_faces(m, lambda p: np.isclose(p[:, 1], W_PLATE, atol=1e-12)),
    ])

    def traction(X, n):
        sxx, syy, sxy = kirsch_cart(X[:, 0], X[:, 1], A_HOLE)
        tx = sxx * n[:, 0] + sxy * n[:, 1]
        ty = sxy * n[:, 0] + syy * n[:, 1]
        return np.stack([tx, ty, np.zeros_like(tx)], axis=1)

    F = T.traction_load(m, faces, traction)
    K = T.stiffness(m, 1.0, 0.0)
    X = m.nodes
    fixed = np.concatenate([
        3 * np.flatnonzero(np.isclose(X[:, 0], 0.0, atol=1e-12)) + 0,   # u_x = 0 on x = 0
        3 * np.flatnonzero(np.isclose(X[:, 1], 0.0, atol=1e-12)) + 1,   # u_y = 0 on y = 0
        3 * np.flatnonzero(np.isclose(X[:, 2], 0.0, atol=1e-12)) + 2,   # u_z = 0 on z = 0
    ])
    u = T.solve(K, F, fixed, 0.0)
    S = T.nodal_stress(m, u, 1.0, 0.0)
    return m, S


@pytest.fixture(scope="module")
def kirsch_meshes():
    return {n: _kirsch_solve(n, n) for n in (8, 16, 32)}


# ---------------------------------------------------------------- T8.1 patch

def test_T8_1_patch_test_linear_field_exact():
    rng = np.random.default_rng(7)
    m = M.perturb_interior(M.box_mesh(1.0, 1.0, 1.0, 3, 3, 3), amp=0.15, seed=3)
    A = rng.normal(size=(3, 3))
    c = rng.normal(size=3)
    E, nu = 1.0, 0.3
    bnd = M.boundary_nodes(m)
    field = lambda X: X @ A.T + c  # noqa: E731
    vals = field(m.nodes[bnd])
    fixed = (3 * bnd[:, None] + np.arange(3)).ravel()
    K = T.stiffness(m, E, nu)
    u = T.solve(K, np.zeros(3 * m.n_nodes), fixed, vals.ravel())
    interior = np.setdiff1d(np.arange(m.n_nodes), bnd)
    u_nodes = u.reshape(-1, 3)
    assert np.max(np.abs(u_nodes[interior] - field(m.nodes[interior]))) <= 1e-10

    sig = T.gauss_stress(m, u, E, nu)               # (m, q, 6)
    sig_exact = _D_ref(E, nu) @ _voigt_strain(A)
    err = np.max(np.abs(sig - sig_exact))
    assert err <= 1e-8 * np.max(np.abs(sig_exact))


# ---------------------------------------------------------------- T8.2 rigid modes

def test_T8_2_single_tet_six_rigid_body_modes():
    m = M.tet_mesh([(0, 0, 0), (1.0, 0, 0), (0, 1.0, 0), (0.2, 0.3, 1.1)])
    K = T.stiffness(m, 1.0, 0.3).toarray()
    lam = np.linalg.eigvalsh(0.5 * (K + K.T))
    lam_abs = np.sort(np.abs(lam))
    assert np.all(lam_abs[:6] < 1e-10 * lam_abs.max())
    assert np.all(lam_abs[6:] > 1e-6 * lam_abs.max())


# ---------------------------------------------------------------- T8.3 cantilever

def _cantilever_tip(nx, ny, nz):
    L, h, w = 10.0, 1.0, 1.0
    E, nu = 1000.0, 0.0
    m = M.box_mesh(L, h, w, nx, ny, nz)
    X = m.nodes
    fixed = (3 * np.flatnonzero(np.isclose(X[:, 0], 0.0))[:, None] + np.arange(3)).ravel()
    tip_faces = M.boundary_faces(m, lambda p: np.isclose(p[:, 0], L))
    area = h * w
    F = T.traction_load(m, tip_faces, lambda P, n: np.tile([0.0, 1.0 / area, 0.0], (P.shape[0], 1)))
    K = T.stiffness(m, E, nu)
    u = T.solve(K, F, fixed, 0.0).reshape(-1, 3)
    tip_nodes = np.flatnonzero(np.isclose(X[:, 0], L))
    return u[tip_nodes, 1].mean(), F


def test_T8_3_cantilever_tip_deflection_and_convergence():
    # P = 1, L = 10, E = 1000, nu = 0 (G = 500), h = w = 1 (A = 1, I = 1/12), kappa = 5/6
    target = 1.0 * 10.0**3 / (3 * 1000.0 * (1.0 / 12.0)) + 1.0 * 10.0 / ((5.0 / 6.0) * 500.0)
    assert target == pytest.approx(4.024, abs=1e-12)
    errs = []
    for nx, ny, nz in [(10, 2, 2), (20, 2, 2), (40, 4, 4)]:
        tip, F = _cantilever_tip(nx, ny, nz)
        # total applied force must equal P = 1
        assert F[1::3].sum() == pytest.approx(1.0, abs=1e-12)
        errs.append(abs(tip - target) / target)
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), errs
    assert errs[-1] < 2e-2, errs


# ---------------------------------------------------------------- T8.4 Kirsch SCF

def test_T8_4_kirsch_scf_convergence(kirsch_meshes):
    errs = []
    for n in (8, 16, 32):
        m, S = kirsch_meshes[n]
        k = np.argmin(np.linalg.norm(m.nodes - np.array([0.0, A_HOLE, 0.0]), axis=1))
        sxx = S[k, 0]
        errs.append(abs(sxx - 3.0) / 3.0)
    print("Kirsch relative SCF errors (8, 16, 32):", errs)
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), errs
    assert errs[-1] < 0.03, errs


# ---------------------------------------------------------------- T8.5 Kirsch hoop stress

def test_T8_5_kirsch_hoop_stress_on_hole(kirsch_meshes):
    m, S = kirsch_meshes[32]
    angles = np.deg2rad([0.0, 22.5, 45.0, 67.5, 90.0])
    for th in angles:
        p = np.array([A_HOLE * np.cos(th), A_HOLE * np.sin(th), 0.0])
        k = np.argmin(np.linalg.norm(m.nodes - p, axis=1))
        assert np.linalg.norm(m.nodes[k] - p) < 1e-9
        s, c = np.sin(th), np.cos(th)
        stt = S[k, 0] * s * s + S[k, 1] * c * c - 2 * S[k, 5] * s * c
        exact = 1.0 - 2.0 * np.cos(2 * th)                  # sigma_thth at r = a, sigma_inf = 1
        assert stt == pytest.approx(exact, rel=5e-2, abs=1e-12)


# ---------------------------------------------------------------- T8.6 element volume

def test_T8_6_element_volume_exact():
    verts = [(0, 0, 0), (2.0, 0, 0), (0, 1.0, 0), (0, 0, 3.0)]   # V = 1
    m = M.tet_mesh(verts)
    vol = T.element_volumes(m)
    assert vol[0] == pytest.approx(1.0, abs=1e-14)
