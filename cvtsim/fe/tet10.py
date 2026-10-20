"""Linear static FE with 10-node tetrahedra (docs/theory/08_fe_tet10.md, sections 2 and 6).

DOF numbering: node n has DOFs 3n, 3n+1, 3n+2 for (u, v, w).
Voigt order for stress and strain: (xx, yy, zz, yz, xz, xy); shear strains are
engineering strains.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from .mesh import TET_EDGES, Mesh, _face_table

# ---- quadrature -------------------------------------------------------------
_a = 0.5854101966249685
_b = 0.1381966011250105
# 4-point tet rule (exact for degree 2). Barycentric (L1..L4) permutations of (a,b,b,b).
GAUSS4_NAT = np.array([[_b, _b, _b], [_a, _b, _b], [_b, _a, _b], [_b, _b, _a]])
GAUSS4_W = np.full(4, 1.0 / 24.0)  # sums to the reference volume 1/6

# 6-point triangle rule (Dunavant degree 4). Weights sum to the reference area 1/2.
_t1, _t2 = 0.445948490915965, 0.091576213509771
_w1, _w2 = 0.223381589678011 / 2.0, 0.109951743655322 / 2.0
_c1, _c2 = 1 - 2 * _t1, 1 - 2 * _t2
TRI6_NAT = np.array([
    [_t1, _t1], [_t1, _c1], [_c1, _t1],
    [_t2, _t2], [_t2, _c2], [_c2, _t2],
])
TRI6_W = np.array([_w1, _w1, _w1, _w2, _w2, _w2])

_DL_TET = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
_DL_TRI = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
_TRI_EDGES = ((0, 1), (1, 2), (0, 2))
_CHUNK = 4096


def _quad_shape(L, dL, edges):
    """Quadratic Lagrange shape functions in volume coordinates L (q, nv)."""
    q, nv = L.shape
    dim = dL.shape[1]
    N = np.empty((q, nv + len(edges)))
    dN = np.empty((q, nv + len(edges), dim))
    for i in range(nv):
        N[:, i] = L[:, i] * (2.0 * L[:, i] - 1.0)
        dN[:, i, :] = (4.0 * L[:, i] - 1.0)[:, None] * dL[i]
    for k, (i, j) in enumerate(edges):
        N[:, nv + k] = 4.0 * L[:, i] * L[:, j]
        dN[:, nv + k, :] = 4.0 * (L[:, i, None] * dL[j] + L[:, j, None] * dL[i])
    return N, dN


def tet_shape(xi):
    """Shape functions and natural derivatives at points xi (q, 3)."""
    xi = np.atleast_2d(np.asarray(xi, dtype=float))
    L = np.concatenate([1.0 - xi.sum(axis=1, keepdims=True), xi], axis=1)
    return _quad_shape(L, _DL_TET, [tuple(e) for e in TET_EDGES])


def tri_shape(st):
    """6-node triangle shape functions at natural points (s, t) of shape (q, 2)."""
    st = np.atleast_2d(np.asarray(st, dtype=float))
    L = np.concatenate([1.0 - st.sum(axis=1, keepdims=True), st], axis=1)
    return _quad_shape(L, _DL_TRI, _TRI_EDGES)


def _D(E: float, nu: float) -> np.ndarray:
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] = lam + 2.0 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def _grads(Xe, dN):
    """Jacobian dets and physical shape-function gradients.

    Xe: (m, 10, 3); dN: (q, 10, 3). Returns det (m, q), G (m, q, 10, 3) with
    G[..., n, i] = dN_n / dx_i.
    """
    J = np.einsum("eni,pnj->epij", Xe, dN)
    det = np.linalg.det(J)
    if np.any(det <= 0.0):
        raise ValueError("non-positive Jacobian: element is inverted or degenerate")
    invJ = np.linalg.inv(J)
    G = np.einsum("pnj,epji->epni", dN, invJ)
    return det, G


def _B(G):
    """Strain-displacement matrices (..., 6, 30) from gradients (..., 10, 3)."""
    sh = G.shape[:-2]
    B = np.zeros(sh + (6, 10, 3))
    B[..., 0, :, 0] = G[..., :, 0]
    B[..., 1, :, 1] = G[..., :, 1]
    B[..., 2, :, 2] = G[..., :, 2]
    B[..., 3, :, 1] = G[..., :, 2]
    B[..., 3, :, 2] = G[..., :, 1]
    B[..., 4, :, 0] = G[..., :, 2]
    B[..., 4, :, 2] = G[..., :, 0]
    B[..., 5, :, 0] = G[..., :, 1]
    B[..., 5, :, 1] = G[..., :, 0]
    return B.reshape(sh + (6, 30))


def element_volumes(mesh: Mesh, rule=(GAUSS4_NAT, GAUSS4_W)) -> np.ndarray:
    """Sum_g w_g det J_g for each element (volume for straight-edged elements)."""
    pts, w = rule
    _, dN = tet_shape(pts)
    X = mesh.nodes[mesh.elems]
    det, _ = _grads(X, dN)
    return det @ w


def stiffness(mesh: Mesh, E: float, nu: float, rule=(GAUSS4_NAT, GAUSS4_W)):
    """Assembled stiffness (3N x 3N, CSR)."""
    pts, w = rule
    _, dN = tet_shape(pts)
    D = _D(E, nu)
    ndof = 3 * mesh.n_nodes
    blocks = []
    for s0 in range(0, mesh.n_elems, _CHUNK):
        el = mesh.elems[s0:s0 + _CHUNK]
        Xe = mesh.nodes[el]
        det, G = _grads(Xe, dN)
        B = _B(G)                                              # (m, q, 6, 30)
        DB = np.einsum("ab,epbj->epaj", D, B)
        Ke = np.einsum("epai,epaj,ep->eij", B, DB, det * w[None, :])
        dofs = (3 * el[:, :, None] + np.arange(3)).reshape(len(el), 30)
        R = np.broadcast_to(dofs[:, :, None], Ke.shape).ravel()
        C = np.broadcast_to(dofs[:, None, :], Ke.shape).ravel()
        blocks.append(sp.coo_matrix((Ke.ravel(), (R, C)), shape=(ndof, ndof)).tocsr())
    K = blocks[0]
    for bk in blocks[1:]:
        K = K + bk
    return K.tocsr()


def traction_load(mesh: Mesh, faces, traction_fn) -> np.ndarray:
    """Consistent nodal loads from surface tractions on 6-node triangle faces.

    faces: (f, 6) node ids (corners then mid-edge nodes), e.g. from boundary_faces.
    traction_fn(points (k,3), outward_normals (k,3)) -> (k,3) traction per area.
    Normals are oriented outward using the element that owns each face.
    """
    faces =np.asarray(faces, dtype=np.int64).reshape(-1, 6)
    F = np.zeros(3 * mesh.n_nodes)
    if faces.shape[0] == 0:
        return F
    keys_all, _, opp_all = _face_table(mesh)
    order = np.argsort(keys_all)
    sk = keys_all[order]
    N_ = mesh.n_nodes
    c = np.sort(faces[:, :3], axis=1)
    fkeys = (c[:, 0] * N_ + c[:, 1]) * N_ + c[:, 2]
    pos = np.searchsorted(sk, fkeys)
    if np.any(pos >= sk.size) or np.any(sk[np.minimum(pos, sk.size - 1)] != fkeys):
        raise ValueError("face not found in mesh")
    opp = mesh.nodes[opp_all[order][pos]]

    Ns, dNs = tri_shape(TRI6_NAT)
    X = mesh.nodes[faces]                                    # (f, 6, 3)
    xs = np.einsum("fni,pn->fpi", X, dNs[:, :, 0])
    xt = np.einsum("fni,pn->fpi", X, dNs[:, :, 1])
    cr = np.cross(xs, xt)
    jac = np.linalg.norm(cr, axis=2)                         # (f, q)
    nrm = cr / jac[..., None]
    xq = np.einsum("fni,pn->fpi", X, Ns)
    sgn = np.sign(np.einsum("fpi,fpi->fp", nrm, xq - opp[:, None, :]))
    nrm = nrm * sgn[..., None]
    t = np.asarray(traction_fn(xq.reshape(-1, 3), nrm.reshape(-1, 3)), dtype=float)
    t = t.reshape(xq.shape)
    vals = np.einsum("pn,fpc,fp->fnc", Ns, t, jac * TRI6_W[None, :])
    idx = (3 * faces[:, :, None] + np.arange(3)).reshape(-1)
    np.add.at(F, idx, vals.reshape(-1))
    return F


def solve(K, F, fixed_dofs, fixed_vals=0.0) -> np.ndarray:
    """Solve K u = F with Dirichlet DOFs eliminated by row/column partitioning."""
    K = sp.csr_matrix(K)
    n = K.shape[0]
    fixed = np.asarray(fixed_dofs, dtype=np.int64).ravel()
    vals = np.broadcast_to(np.asarray(fixed_vals, dtype=float), fixed.shape)
    free = np.ones(n, dtype=bool)
    free[fixed] = False
    f = np.flatnonzero(free)
    u = np.zeros(n)
    u[fixed] = vals
    rhs = F[f] - K[f][:, fixed] @ vals
    u[f] = spsolve(K[f][:, f].tocsc(), rhs)
    return u


def _stress_at(mesh: Mesh, u, E, nu, pts):
    """Stress at natural points pts (q, 3) for every element: (m, q, 6)."""
    _, dN = tet_shape(pts)
    D = _D(E, nu)
    X = mesh.nodes[mesh.elems]
    ue = np.asarray(u, dtype=float).reshape(-1, 3)[mesh.elems]
    _, G = _grads(X, dN)
    gu = np.einsum("eni,epnj->epij", ue, G)                  # gu[i, j] = du_i / dx_j
    eps = np.stack([
        gu[..., 0, 0], gu[..., 1, 1], gu[..., 2, 2],
        gu[..., 1, 2] + gu[..., 2, 1],
        gu[..., 0, 2] + gu[..., 2, 0],
        gu[..., 0, 1] + gu[..., 1, 0],
    ], axis=-1)
    return np.einsum("ab,epb->epa", D, eps)


def gauss_stress(mesh: Mesh, u, E: float, nu: float, rule=(GAUSS4_NAT, GAUSS4_W)):
    """Stress at Gauss points: (m, q, 6)."""
    return _stress_at(mesh, u, E, nu, rule[0])


_NAT_CORNER = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
# natural coordinates of the 10 nodes (mid-edge nodes at edge midpoints)
_NAT10 = np.vstack([_NAT_CORNER, [0.5 * (_NAT_CORNER[a] + _NAT_CORNER[b]) for a, b in TET_EDGES]])


def nodal_stress(mesh: Mesh, u, E: float, nu: float) -> np.ndarray:
    """Nodal stress: each element's stress at its own nodes, then averaged (unweighted)
    over the elements sharing each node. Returns (N, 6)."""
    s = _stress_at(mesh, u, E, nu, _NAT10)                  # (m, 10, 6)
    out = np.zeros((mesh.n_nodes, 6))
    np.add.at(out, mesh.elems.ravel(), s.reshape(-1, 6))
    cnt = np.bincount(mesh.elems.ravel(), minlength=mesh.n_nodes).astype(float)
    return out / cnt[:, None]
