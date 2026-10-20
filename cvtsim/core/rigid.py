"""Rigid bodies and quaternion helpers (docs/theory/01_rigid_body.md, sec. 2-3).

Quaternions are stored as (w, x, y, z) and map body to world: R = quat_to_R(q).
"""
from __future__ import annotations

import math

import numpy as np


def skew(a):
    """Skew-symmetric matrix [a]x such that [a]x b = a x b."""
    a = np.asarray(a, dtype=float)
    return np.array([[0.0, -a[2], a[1]],
                     [a[2], 0.0, -a[0]],
                     [-a[1], a[0], 0.0]])


def cross3(a, b):
    """Cross product of two 3-vectors (faster than np.cross for small arrays)."""
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])


def quat_normalize(q):
    q = np.asarray(q, dtype=float)
    return q / np.linalg.norm(q)


def quat_mul(p, q):
    """Hamilton product p (x) q."""
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def quat_exp(rotvec):
    """Unit quaternion of the rotation by |rotvec| about rotvec/|rotvec|."""
    v = np.asarray(rotvec, dtype=float).reshape(3)
    th = np.linalg.norm(v)
    if th < 1e-12:
        return quat_normalize(np.array([1.0, 0.5 * v[0], 0.5 * v[1], 0.5 * v[2]]))
    s = np.sin(0.5 * th) / th
    return np.array([np.cos(0.5 * th), s * v[0], s * v[1], s * v[2]])


def quat_to_R(q):
    """Body-to-world rotation matrix of (possibly non-unit) quaternion q."""
    w, x, y, z = (float(c) for c in q)
    n = math.sqrt(w * w + x * x + y * y + z * z)
    w, x, y, z = w / n, x / n, y / n, z / n
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


class RigidBody:
    """Rigid body with CoM state in the world frame.

    Parameters: mass ``m`` [kg], body-frame CoM inertia ``Jb`` [kg m^2],
    CoM position ``r`` [m], quaternion ``q`` (w,x,y,z), CoM velocity ``v`` [m/s]
    and world-frame angular velocity ``w`` [rad/s].

    ``force`` [N] and ``torque`` [N m, about CoM, world frame] are applied
    loads; they default to zero and are added to the body's generalized force.
    """

    def __init__(self, m, Jb, r, q, v, w):
        self.m = float(m)
        self.Jb = np.asarray(Jb, dtype=float).reshape(3, 3)
        self.r = np.asarray(r, dtype=float).reshape(3)
        self.q = np.asarray(q, dtype=float).reshape(4)   # normalized by the setter
        self.v = np.asarray(v, dtype=float).reshape(3)
        self.w = np.asarray(w, dtype=float).reshape(3)
        self.force = np.zeros(3)
        self.torque = np.zeros(3)

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, value):
        self._q = quat_normalize(np.asarray(value, dtype=float).reshape(4))
        self._R = None
        self._Rt = None

    @property
    def R(self):
        """Body-to-world rotation (cached until q is set again)."""
        if self._R is None:
            self._R = quat_to_R(self._q)
        return self._R

    @property
    def Rt(self):
        """Row-major 9-tuple of R (cached)."""
        if self._Rt is None:
            self._Rt = tuple(self.R.ravel().tolist())
        return self._Rt

    @property
    def Jw(self):
        """World-frame inertia tensor R Jb R^T."""
        R = self.R
        return R @ self.Jb @ R.T
