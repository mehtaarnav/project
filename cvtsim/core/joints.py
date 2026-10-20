"""Holonomic joints (docs/theory/01_rigid_body.md, sec. 4).

Conventions
-----------
* ``bi`` / ``bj`` are body indices into ``System.bodies``; ``None`` is ground.
  Ground has R = I, r = v = w = 0, so an attachment ``s`` on ground is a world point.
* ``si``, ``sj`` are attachment points in the body frame, relative to the CoM.
* ``d = r_j + rho_j - r_i - rho_i`` (rho = R s). Spherical uses
  ``Phi = r_i + rho_i - r_j - rho_j`` (= -d), matching its Jacobian block.
* ``evaluate(sys)`` returns ``(Phi (m,), Jl (m,12), gamma (m,))`` where ``Jl``
  holds the body-i block ``[v_i, w_i]`` and the body-j block ``[v_j, w_j]``.

The row algebra is written with plain Python floats: for 3-vectors this is
several times faster than small numpy calls, and the joints are evaluated
several times per time step.
"""
from __future__ import annotations

import math

import numpy as np

_ZERO = (0.0, 0.0, 0.0)
_GROUND_R = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)


# ---- small-vector algebra on 3-tuples ---------------------------------------

def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _neg(a):
    return (-a[0], -a[1], -a[2])


def _scale(s, a):
    return (s * a[0], s * a[1], s * a[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _mv(R, x):
    """R is a row-major 3x3 as a 9-sequence."""
    return (R[0] * x[0] + R[1] * x[1] + R[2] * x[2],
            R[3] * x[0] + R[4] * x[1] + R[5] * x[2],
            R[6] * x[0] + R[7] * x[1] + R[8] * x[2])


def _unit(a):
    a = np.asarray(a, dtype=float).reshape(3)
    n = float(np.linalg.norm(a))
    return (a[0] / n, a[1] / n, a[2] / n)


def _tup(a):
    return (float(a[0]), float(a[1]), float(a[2]))


def perp_basis(a):
    """Return (u, w) orthonormal and perpendicular to unit a, with u x w = a."""
    a = _unit(a)
    t = (1.0, 0.0, 0.0) if abs(a[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _cross(a, t)
    n = math.sqrt(_dot(u, u))
    u = _scale(1.0 / n, u)
    w = _cross(a, u)
    return u, w


class _Frame:
    """Kinematic data of one side of a joint (ground or a body), tuples of floats."""

    __slots__ = ('r', 'R', 'v', 'w', 'rho')

    def __init__(self, r, R, v, w, rho):
        self.r = r
        self.R = R
        self.v = v
        self.w = w
        self.rho = rho


# ---- row builders: each returns (phi, Ji (6), Jj (6), g) for one scalar row -

def _spherical_rows(A, B):
    """Three rows; returns list of (phi, Ji, Jj, g) with Ji, Jj 3-row lists of 6-tuples."""
    rho_i, rho_j = A.rho, B.rho
    phi = _sub(_add(A.r, rho_i), _add(B.r, rho_j))
    wi, wj = A.w, B.w
    g = _neg(_sub(_cross(wi, _cross(wi, rho_i)), _cross(wj, _cross(wj, rho_j))))
    # row k of -[rho_i]x is (..) ; written out explicitly
    ri, rj = rho_i, rho_j
    Ji = [(1.0, 0.0, 0.0, 0.0, ri[2], -ri[1]),
          (0.0, 1.0, 0.0, -ri[2], 0.0, ri[0]),
          (0.0, 0.0, 1.0, ri[1], -ri[0], 0.0)]
    Jj = [(-1.0, 0.0, 0.0, 0.0, -rj[2], rj[1]),
          (0.0, -1.0, 0.0, rj[2], 0.0, -rj[0]),
          (0.0, 0.0, -1.0, -rj[1], rj[0], 0.0)]
    return [(phi[k], Ji[k], Jj[k], g[k]) for k in range(3)]


def _dot1(A, B, a, b):
    """Phi = a . b with a on body i and b on body j (world vectors)."""
    wi, wj = A.w, B.w
    c = _cross(a, b)
    phi = _dot(a, b)
    Ji = (0.0, 0.0, 0.0) + c
    Jj = (0.0, 0.0, 0.0) + _neg(c)
    g = -(_dot(_cross(wi, _cross(wi, a)), b)
          + 2.0 * _dot(_cross(wi, a), _cross(wj, b))
          + _dot(a, _cross(wj, _cross(wj, b))))
    return (phi, Ji, Jj, g)


def _dot2(A, B, d, dd, a):
    """Phi = a . d with a on body i (world vector), d = r_j+rho_j-r_i-rho_i."""
    rho_i, rho_j = A.rho, B.rho
    wi, wj = A.w, B.w
    phi = _dot(a, d)
    Ji = _neg(a) + _sub(_cross(a, d), _cross(rho_i, a))
    Jj = a + _cross(rho_j, a)
    g = -(_dot(a, _sub(_cross(wj, _cross(wj, rho_j)), _cross(wi, _cross(wi, rho_i))))
          + 2.0 * _dot(_cross(wi, a), dd)
          + _dot(d, _cross(wi, _cross(wi, a))))
    return (phi, Ji, Jj, g)


def _pack(rows):
    """rows: list of (phi, Ji, Jj, g) with scalar phi/g and 6-tuple Ji/Jj."""
    phi = np.array([r[0] for r in rows], dtype=float)
    Jl = np.array([tuple(r[1]) + tuple(r[2]) for r in rows], dtype=float).reshape(-1, 12)
    g = np.array([r[3] for r in rows], dtype=float)
    return phi, Jl, g


class Joint:
    """Base class. Subclasses set ``nrows`` and implement ``evaluate``."""

    nrows = 0

    def __init__(self, bi, bj, si=(0.0, 0.0, 0.0), sj=(0.0, 0.0, 0.0)):
        self.bi = None if bi is None else int(bi)
        self.bj = None if bj is None else int(bj)
        self.si = _tup(np.asarray(si, dtype=float).reshape(3))
        self.sj = _tup(np.asarray(sj, dtype=float).reshape(3))

    def initialize(self, sys):
        """Hook called once by System.__init__."""

    def commit(self, sys):
        """Hook called after each accepted step (screw angle bookkeeping)."""

    def _frame(self, sys, idx, s):
        if idx is None:
            return _Frame(_ZERO, _GROUND_R, _ZERO, _ZERO, s)
        b = sys.bodies[idx]
        Rt = b.Rt
        return _Frame(_tup(b.r), Rt, _tup(b.v), _tup(b.w), _mv(Rt, s))

    def _kin(self, sys):
        A = self._frame(sys, self.bi, self.si)
        B = self._frame(sys, self.bj, self.sj)
        d = _sub(_add(B.r, B.rho), _add(A.r, A.rho))
        dd = _sub(_add(B.v, _cross(B.w, B.rho)), _add(A.v, _cross(A.w, A.rho)))
        return A, B, d, dd

    def evaluate(self, sys):
        raise NotImplementedError

    def phi(self, sys):
        return self.evaluate(sys)[0]

    def gamma(self, sys):
        return self.evaluate(sys)[2]

    def jac(self, sys):
        """Dense (m, 6*nbodies) Jacobian in the [v, w] per-body layout."""
        _, Jl, _ = self.evaluate(sys)
        J = np.zeros((Jl.shape[0], 6 * len(sys.bodies)))
        if self.bi is not None:
            J[:, 6 * self.bi:6 * self.bi + 6] += Jl[:, :6]
        if self.bj is not None:
            J[:, 6 * self.bj:6 * self.bj + 6] += Jl[:, 6:]
        return J


class Spherical(Joint):
    """Ball joint: 3 positional rows."""

    nrows = 3

    def evaluate(self, sys):
        A, B, _, _ = self._kin(sys)
        return _pack(_spherical_rows(A, B))


class Revolute(Joint):
    """Revolute joint about body-i axis ai (body-j axis aj aligned). 5 rows."""

    nrows = 5

    def __init__(self, bi, bj, si, sj, ai, aj):
        super().__init__(bi, bj, si, sj)
        self.ai = _unit(ai)
        self.aj = _unit(aj)
        self.ui, self.wi = perp_basis(self.ai)

    def evaluate(self, sys):
        A, B, _, _ = self._kin(sys)
        b = _mv(B.R, self.aj)
        rows = _spherical_rows(A, B)
        rows.append(_dot1(A, B, _mv(A.R, self.ui), b))
        rows.append(_dot1(A, B, _mv(A.R, self.wi), b))
        return _pack(rows)


class Prismatic(Joint):
    """Prismatic joint along body-i axis ai. 5 rows.

    Rotation is locked by dot-1 rows on the pairs (u_i, w_j), (u_i, a_j),
    (w_i, a_j); transverse translation by dot-2 rows with u_i, w_i.
    The body triads must be matched at the initial configuration (checked in
    ``initialize``); by default they are generated from ``ai`` and ``aj``.
    """

    nrows = 5

    def __init__(self, bi, bj, si, sj, ai, aj):
        super().__init__(bi, bj, si, sj)
        self.ai = _unit(ai)
        self.aj = _unit(aj)
        self.ui, self.wi = perp_basis(self.ai)
        self.wj = perp_basis(self.aj)[1]

    def _rot_rows(self, A, B):
        ui = _mv(A.R, self.ui)
        wi = _mv(A.R, self.wi)
        aj = _mv(B.R, self.aj)
        wj = _mv(B.R, self.wj)
        rows = [_dot1(A, B, ui, wj), _dot1(A, B, ui, aj), _dot1(A, B, wi, aj)]
        return rows, ui, wi

    def initialize(self, sys):
        A, B, _, _ = self._kin(sys)
        rows, _, _ = self._rot_rows(A, B)
        if max(abs(r[0]) for r in rows) > 1e-9:
            raise ValueError("Prismatic body triads are not matched at the initial configuration")

    def evaluate(self, sys):
        A, B, d, dd = self._kin(sys)
        rows, ui, wi = self._rot_rows(A, B)
        rows.append(_dot2(A, B, d, dd, ui))
        rows.append(_dot2(A, B, d, dd, wi))
        return _pack(rows)


class Screw(Joint):
    """Screw (helix) joint along body-i axis ai with lead ``pitch`` [m/rad]. 5 rows.

    Cylindrical part: dot-1 x2 (axis alignment) and dot-2 x2 (transverse offset).
    Screw row: Phi_s = a.d - p*theta - h0, where theta is the unwrapped
    relative rotation about a. theta = 0 at the configuration given at
    construction (``initialize``), and it is advanced only by ``commit``.
    """

    nrows = 5

    def __init__(self, bi, bj, si, sj, ai, aj, pitch):
        super().__init__(bi, bj, si, sj)
        self.ai = _unit(ai)
        self.aj = _unit(aj)
        self.p = float(pitch)
        self.ui, self.wi = perp_basis(self.ai)
        self.uj = perp_basis(self.aj)[0]
        self._raw_prev = 0.0
        self._acc = 0.0
        self._h0 = 0.0

    def _raw_angle_of(self, A, B):
        a = _mv(A.R, self.ai)
        ei = _mv(A.R, self.ui)
        ej = _mv(B.R, self.uj)
        ej = _sub(ej, _scale(_dot(ej, a), a))
        return math.atan2(_dot(a, _cross(ei, ej)), _dot(ei, ej))

    def _theta_of(self, A, B):
        raw = self._raw_angle_of(A, B)
        dth = (raw - self._raw_prev + math.pi) % (2.0 * math.pi) - math.pi
        return self._acc + dth, raw

    def theta(self, sys):
        """Unwrapped relative rotation about the screw axis (no state change)."""
        A, B, _, _ = self._kin(sys)
        return self._theta_of(A, B)[0]

    def initialize(self, sys):
        A, B, d, _ = self._kin(sys)
        self._raw_prev = self._raw_angle_of(A, B)
        self._acc = 0.0
        self._h0 = _dot(_mv(A.R, self.ai), d)

    def commit(self, sys):
        A, B, _, _ = self._kin(sys)
        theta, raw = self._theta_of(A, B)
        self._acc = theta
        self._raw_prev = raw

    def evaluate(self, sys):
        A, B, d, dd = self._kin(sys)
        a = _mv(A.R, self.ai)
        b = _mv(B.R, self.aj)
        ui = _mv(A.R, self.ui)
        wi_loc = _mv(A.R, self.wi)
        theta, _ = self._theta_of(A, B)
        wi, wj = A.w, B.w

        phi_d, Ji_d, Jj_d, g_d = _dot2(A, B, d, dd, a)
        p = self.p
        phi_s = phi_d - p * theta - self._h0
        # rotational part of the screw row: +p*a on body i, -p*a on body j
        Ji_s = (Ji_d[0], Ji_d[1], Ji_d[2], Ji_d[3] + p * a[0], Ji_d[4] + p * a[1], Ji_d[5] + p * a[2])
        Jj_s = (Jj_d[0], Jj_d[1], Jj_d[2], Jj_d[3] - p * a[0], Jj_d[4] - p * a[1], Jj_d[5] - p * a[2])
        g_s = g_d + p * _dot(_cross(wi, a), _sub(wj, wi))

        rows = [
            _dot1(A, B, ui, b),
            _dot1(A, B, wi_loc, b),
            _dot2(A, B, d, dd, ui),
            _dot2(A, B, d, dd, wi_loc),
            (phi_s, Ji_s, Jj_s, g_s),
        ]
        return _pack(rows)
