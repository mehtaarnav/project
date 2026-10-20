"""Multibody system: constrained Newton-Euler dynamics (docs/theory/01_rigid_body.md, sec. 3).

Generalized velocity u = [v_1, w_1, v_2, w_2, ...] (world frame), 6 per body.
Equations: M udot - Jc^T lam = F, Jc udot = gamma - 2a Jc u - b^2 Phi (Baumgarte).
"""
from __future__ import annotations

import numpy as np

from .rigid import RigidBody, cross3, quat_exp, quat_mul, quat_normalize

_I3 = np.eye(3)


class System:
    """Rigid bodies connected by joints.

    Parameters
    ----------
    bodies : sequence of RigidBody
    joints : sequence of Joint (bodies referenced by index, None = ground)
    gravity : (3,) acceleration [m/s^2]
    baumgarte : (a, b) stabilization gains; (0, 0) gives pure index-1 dynamics
    project : if True, every step is followed by position and velocity projection
    """

    def __init__(self, bodies, joints=(), gravity=(0.0, 0.0, -9.81),
                 baumgarte=(0.0, 0.0), project=False):
        self.bodies = list(bodies)
        self.joints = list(joints)
        self.nb = len(self.bodies)
        self.gravity = np.asarray(gravity, dtype=float).reshape(3)
        self.baumgarte = (float(baumgarte[0]), float(baumgarte[1]))
        self.project = bool(project)
        self.last_lambda = np.zeros(0)
        for jt in self.joints:
            jt.initialize(self)

    # ---- state access ------------------------------------------------------
    def get_u(self):
        return np.concatenate([np.concatenate([b.v, b.w]) for b in self.bodies]) \
            if self.bodies else np.zeros(0)

    def set_u(self, u):
        for k, b in enumerate(self.bodies):
            b.v = u[6 * k:6 * k + 3].copy()
            b.w = u[6 * k + 3:6 * k + 6].copy()

    def _get_state(self):
        return np.concatenate([np.concatenate([b.r, b.q, b.v, b.w]) for b in self.bodies])

    def _set_state(self, y):
        for k, b in enumerate(self.bodies):
            s = y[13 * k:13 * k + 13]
            b.r = s[0:3].copy()
            b.q = quat_normalize(s[3:7])
            b.v = s[7:10].copy()
            b.w = s[10:13].copy()

    # ---- dynamics ----------------------------------------------------------
    def _mass_and_force(self):
        N = 6 * self.nb
        M = np.zeros((N, N))
        F = np.zeros(N)
        for k, b in enumerate(self.bodies):
            Jw = b.Jw
            M[6 * k:6 * k + 3, 6 * k:6 * k + 3] = b.m * _I3
            M[6 * k + 3:6 * k + 6, 6 * k + 3:6 * k + 6] = Jw
            F[6 * k:6 * k + 3] = b.m * self.gravity + b.force
            F[6 * k + 3:6 * k + 6] = b.torque - cross3(b.w, Jw @ b.w)
        return M, F

    def constraints(self):
        """Stacked (Phi, Jc, gamma) and per-joint row slices."""
        N = 6 * self.nb
        evals = [(jt, jt.evaluate(self)) for jt in self.joints]
        m_tot = sum(ev[1][0].size for ev in evals)
        J = np.zeros((m_tot, N))
        if m_tot == 0:
            return np.zeros(0), J, np.zeros(0)
        Phi = np.concatenate([ev[1][0] for ev in evals])
        gam = np.concatenate([ev[1][2] for ev in evals])
        r0 = 0
        for jt, (phi, Jl, _) in evals:
            r1 = r0 + phi.size
            if jt.bi is not None:
                J[r0:r1, 6 * jt.bi:6 * jt.bi + 6] += Jl[:, :6]
            if jt.bj is not None:
                J[r0:r1, 6 * jt.bj:6 * jt.bj + 6] += Jl[:, 6:]
            r0 = r1
        return Phi, J, gam

    def accelerations(self):
        """Solve the KKT system; returns udot. Lagrange multipliers in ``last_lambda``."""
        M, F = self._mass_and_force()
        u = self.get_u()
        Phi, J, gam = self.constraints()
        m = J.shape[0]
        N = M.shape[0]
        if m == 0:
            udot = np.linalg.solve(M, F)
            self.last_lambda = np.zeros(0)
            return udot
        a, b = self.baumgarte
        rhs = gam - 2.0 * a * (J @ u) - b * b * Phi
        K = np.zeros((N + m, N + m))
        K[:N, :N] = M
        K[:N, N:] = -J.T
        K[N:, :N] = J
        sol = np.linalg.solve(K, np.concatenate([F, rhs]))
        self.last_lambda = sol[N:]
        return sol[:N]

    def _step_symplectic_euler(self, h):
        udot = self.accelerations()
        self.set_u(self.get_u() + h * udot)
        for b in self.bodies:
            b.r = b.r + h * b.v
            b.q = quat_normalize(quat_mul(quat_exp(h * b.w), b.q))

    def _step_rk4(self, h):
        def rhs(y):
            self._set_state(y)
            udot = self.accelerations()
            out = np.empty_like(y)
            for k, b in enumerate(self.bodies):
                s = 13 * k
                out[s:s + 3] = b.v
                out[s + 3:s + 7] = 0.5 * quat_mul(np.concatenate([[0.0], b.w]), b.q)
                out[s + 7:s + 10] = udot[6 * k:6 * k + 3]
                out[s + 10:s + 13] = udot[6 * k + 3:6 * k + 6]
            return out

        y0 = self._get_state()
        k1 = rhs(y0)
        k2 = rhs(y0 + 0.5 * h * k1)
        k3 = rhs(y0 + 0.5 * h * k2)
        k4 = rhs(y0 + h * k3)
        y1 = y0 + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        self._set_state(y1)

    def step(self, h, method='symplectic_euler'):
        """Advance the system by one step of size h."""
        if method == 'symplectic_euler':
            self._step_symplectic_euler(h)
        elif method == 'rk4':
            self._step_rk4(h)
        else:
            raise ValueError(f"unknown method {method!r}")
        if self.project:
            self._project()
        for jt in self.joints:
            jt.commit(self)

    # ---- projection --------------------------------------------------------
    def _project(self, max_iter=5, tol=1e-12):
        # Position Gauss-Newton; the Jacobian of the last evaluation is at the final
        # positions, so it is reused for the velocity projection.
        for _ in range(max_iter):
            Phi, J, _ = self.constraints()
            if Phi.size == 0 or np.max(np.abs(Phi)) < tol:
                break
            dq = -J.T @ np.linalg.solve(J @ J.T, Phi)
            for k, b in enumerate(self.bodies):
                b.r = b.r + dq[6 * k:6 * k + 3]
                b.q = quat_normalize(quat_mul(quat_exp(dq[6 * k + 3:6 * k + 6]), b.q))
        else:
            # corrections were applied after the last evaluation: refresh the Jacobian
            _, J, _ = self.constraints()
        if J.shape[0] == 0:
            return
        M, _ = self._mass_and_force()
        u = self.get_u()
        Y = np.linalg.solve(M, J.T)
        lam = np.linalg.solve(J @ Y, J @ u)
        self.set_u(u - Y @ lam)

    # ---- diagnostics -------------------------------------------------------
    def energy(self):
        """Kinetic plus gravitational potential energy, U = -sum m g.r."""
        E = 0.0
        for b in self.bodies:
            E += 0.5 * b.m * float(b.v @ b.v) + 0.5 * float(b.w @ (b.Jw @ b.w))
            E -= b.m * float(self.gravity @ b.r)
        return E

    def constraint_violation(self):
        """Max over joints of the Euclidean norm of each joint's Phi."""
        worst = 0.0
        for jt in self.joints:
            worst = max(worst, float(np.linalg.norm(jt.phi(self))))
        return worst
