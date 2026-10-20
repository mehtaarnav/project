"""10-node tetrahedral meshes and generators (docs/theory/08_fe_tet10.md, section 3).

Local tet10 numbering (0-based): corners 0..3; mid-edge nodes 4..9 sit on
TET_EDGES = (0,1),(1,2),(0,2),(0,3),(1,3),(2,3).
Faces (TET_FACES) list the three corners followed by the mid-edge nodes on
edges (ab), (bc), (ca), matching the 6-node triangle ordering.
"""
from __future__ import annotations

import numpy as np

TET_EDGES = np.array([(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)], dtype=np.int64)
TET_FACES = np.array([
    [0, 1, 2, 4, 5, 6],
    [0, 1, 3, 4, 8, 7],
    [0, 2, 3, 6, 9, 7],
    [1, 2, 3, 5, 9, 8],
], dtype=np.int64)
TET_FACE_OPP = np.array([3, 2, 1, 0], dtype=np.int64)  # vertex opposite each face

# Kuhn (Freudenthal) split of a hex whose corners are numbered
# v0=(0,0,0) v1=(1,0,0) v2=(0,1,0) v3=(1,1,0) v4=(0,0,1) v5=(1,0,1) v6=(0,1,1) v7=(1,1,1).
_HEX_OFFSETS = np.array([(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
                         (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)], dtype=np.int64)
_KUHN = np.array([[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7],
                  [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]], dtype=np.int64)


class Mesh:
    """Quadratic tetrahedral mesh.

    nodes : (N, 3) coordinates; the first `n_corner` rows are corner nodes.
    elems : (M, 10) connectivity in the local numbering above.
    mid_edges : (N - n_corner, 2) corner pairs of each mid-edge node, if known.
    hole_radius : if set, mid-edge nodes whose two corners lie on the circle
        r = hole_radius (about the z axis) are projected onto that circle.
    """

    def __init__(self, nodes, elems, n_corner=None, mid_edges=None, hole_radius=None):
        self.nodes = np.asarray(nodes, dtype=float).reshape(-1, 3)
        self.elems = np.asarray(elems, dtype=np.int64).reshape(-1, 10)
        self.n_corner = self.nodes.shape[0] if n_corner is None else int(n_corner)
        self.mid_edges = None if mid_edges is None else np.asarray(mid_edges, dtype=np.int64)
        self.hole_radius = hole_radius

    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]

    @property
    def n_elems(self) -> int:
        return self.elems.shape[0]


def _place_mids(X, mid_edges, n_corner, hole_radius):
    a = X[mid_edges[:, 0]]
    b = X[mid_edges[:, 1]]
    mid = 0.5 * (a + b)
    if hole_radius is not None:
        tol = 1e-9 * max(1.0, hole_radius)
        on_a = np.abs(np.hypot(a[:, 0], a[:, 1]) - hole_radius) < tol
        on_b = np.abs(np.hypot(b[:, 0], b[:, 1]) - hole_radius) < tol
        both = on_a & on_b
        rxy = np.hypot(mid[both, 0], mid[both, 1])
        mid[both, 0:2] *= (hole_radius / rxy)[:, None]
    X[n_corner:] = mid
    return X


def _face_table(mesh: Mesh):
    """All element faces: sorted-corner keys, 6-node face connectivity, opposite vertex."""
    el = mesh.elems
    faces = el[:, TET_FACES].reshape(-1, 6)
    opp = el[:, TET_FACE_OPP].reshape(-1)
    c = np.sort(faces[:, :3], axis=1)
    N = mesh.n_nodes
    keys = (c[:, 0] * N + c[:, 1]) * N + c[:, 2]
    return keys, faces, opp


def _boundary(mesh: Mesh):
    keys, faces, opp = _face_table(mesh)
    uniq, inv, counts = np.unique(keys, return_inverse=True, return_counts=True)
    bmask = counts[inv] == 1
    return faces[bmask], opp[bmask], keys[bmask]


def boundary_faces(mesh: Mesh, predicate):
    """6-node faces on the outer boundary whose six nodes all satisfy predicate.

    predicate: callable taking an (k, 3) array of coordinates, returns (k,) bool.
    """
    faces, _, _ = _boundary(mesh)
    if faces.size == 0:
        return faces.reshape(0, 6)
    ok = np.asarray(predicate(mesh.nodes[faces.ravel()]), dtype=bool).reshape(-1, 6)
    return faces[ok.all(axis=1)]


def boundary_nodes(mesh: Mesh) -> np.ndarray:
    faces, _, _ = _boundary(mesh)
    return np.unique(faces)


def _build(xyz, tets, hole_radius=None) -> Mesh:
    """Build a tet10 mesh from corner coordinates and corner tetrahedra."""
    xyz = np.asarray(xyz, dtype=float)
    tets = np.asarray(tets, dtype=np.int64).copy()
    Nc = xyz.shape[0]
    p = xyz[tets]
    vol6 = np.einsum("mi,mi->m", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    neg = vol6 < 0
    tets[neg] = tets[neg][:, [0, 2, 1, 3]]          # make every tet positively oriented
    pairs = np.sort(tets[:, TET_EDGES], axis=2)       # (m, 6, 2)
    keys = pairs[..., 0] * Nc + pairs[..., 1]
    uniq, inv = np.unique(keys.ravel(), return_inverse=True)
    mid_ids = Nc + inv.reshape(-1, 6)
    mid_edges = np.stack([uniq // Nc, uniq % Nc], axis=1)
    X = np.vstack([xyz, np.zeros((len(uniq), 3))])
    X = _place_mids(X, mid_edges, Nc, hole_radius)
    elems = np.hstack([tets, mid_ids])
    return Mesh(X, elems, n_corner=Nc, mid_edges=mid_edges, hole_radius=hole_radius)


def tet_mesh(vertices) -> Mesh:
    """Single tet10 with straight edges."""
    v = np.asarray(vertices, dtype=float).reshape(4, 3)
    return _build(v, [[0, 1, 2, 3]])


def box_mesh(Lx, Ly, Lz, nx, ny, nz) -> Mesh:
    """Structured box [0,Lx]x[0,Ly]x[0,Lz]; each hex split into 6 Kuhn tets."""
    xs = np.linspace(0.0, Lx, nx + 1)
    ys = np.linspace(0.0, Ly, ny + 1)
    zs = np.linspace(0.0, Lz, nz + 1)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    xyz = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
    gid = lambda i, j, k: (i * (ny + 1) + j) * (nz + 1) + k  # noqa: E731
    I, J, K = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    I, J, K = I.ravel(), J.ravel(), K.ravel()
    hexes = np.stack([gid(I + di, J + dj, K + dk) for di, dj, dk in _HEX_OFFSETS], axis=1)
    tets = hexes[:, _KUHN].reshape(-1, 4)
    return _build(xyz, tets)


def plate_hole_quarter_mesh(a, W, t, n_r, n_theta, n_z, grading=1.0) -> Mesh:
    """Quarter plate [0,W]^2 minus a hole of radius a, thickness t (z in [0,t]).

    Structured polar grid: theta in [0, 90 deg] split into two blocks of n_theta
    angular elements each (block 1 reaches edge x=W, block 2 reaches edge y=W).
    Radial element sizes grow geometrically with ratio `grading` (>1 concentrates
    elements at the hole). Hexes are split into Kuhn tets; mid-edge nodes on the
    hole are projected onto the circle.
    """
    if not (0 < a < W and t > 0 and n_r >= 1 and n_theta >= 1 and n_z >= 1):
        raise ValueError("invalid plate-with-hole dimensions")
    th = np.linspace(0.0, np.pi / 2, 2 * n_theta + 1)
    Rout = np.empty_like(th)
    blk1 = th <= np.pi / 4 + 1e-14
    Rout[blk1] = W / np.cos(th[blk1])
    Rout[~blk1] = W / np.sin(th[~blk1])
    sizes = float(grading) ** np.arange(n_r)
    frac = np.concatenate([[0.0], np.cumsum(sizes)]) / sizes.sum()   # (n_r+1,)
    r = a + (Rout[None, :] - a) * frac[:, None]                       # (n_r+1, 2nθ+1)
    xy = np.stack([r * np.cos(th)[None, :], r * np.sin(th)[None, :]], axis=-1)  # (i,j,2)
    zc = t * np.arange(n_z + 1) / n_z
    xyz = np.empty((n_r + 1, 2 * n_theta + 1, n_z + 1, 3))
    xyz[..., 0:2] = xy[:, :, None, :]
    xyz[..., 2] = zc[None, None, :]
    xyz = xyz.reshape(-1, 3)
    NJ = 2 * n_theta + 1
    gid = lambda i, j, k: (i * NJ + j) * (n_z + 1) + k  # noqa: E731
    I, Jj, Kk = np.meshgrid(np.arange(n_r), np.arange(2 * n_theta), np.arange(n_z), indexing="ij")
    I, Jj, Kk = I.ravel(), Jj.ravel(), Kk.ravel()
    hexes = np.stack([gid(I + di, Jj + dj, Kk + dk) for di, dj, dk in _HEX_OFFSETS], axis=1)
    tets = hexes[:, _KUHN].reshape(-1, 4)
    return _build(xyz, tets, hole_radius=a)


def perturb_interior(mesh: Mesh, amp: float, seed: int) -> Mesh:
    """Random displacement of interior corner nodes (amp times the smallest edge).

    Mid-edge nodes are re-placed at the midpoints of their (moved) corners, or
    projected onto the hole circle if the mesh has one. Boundary nodes stay fixed.
    """
    rng = np.random.default_rng(seed)
    X = mesh.nodes.copy()
    bnd = set(boundary_nodes(mesh).tolist())
    interior = np.array([i for i in range(mesh.n_corner) if i not in bnd], dtype=np.int64)
    corners = mesh.elems[:, :4]
    p = X[corners]
    edge_len = np.stack([np.linalg.norm(p[:, a] - p[:, b], axis=1) for a, b in TET_EDGES], axis=1)
    h = edge_len.min()
    X[interior] += rng.uniform(-1.0, 1.0, size=(interior.size, 3)) * (amp * h)
    if mesh.mid_edges is not None:
        X = _place_mids(X, mesh.mid_edges, mesh.n_corner, mesh.hole_radius)
    return Mesh(X, mesh.elems, n_corner=mesh.n_corner, mid_edges=mesh.mid_edges,
                hole_radius=mesh.hole_radius)
