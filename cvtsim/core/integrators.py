"""Time integrators (docs/theory/02_integrators.md).

* Generalized-alpha (Chung & Hulbert 1993) for M x'' + f_int(x, x', t) = f_ext(t).
* Velocity Verlet and symplectic Euler for x'' = accel(x, v, t).
"""
from __future__ import annotations

import numpy as np


def genalpha_params(rho_inf):
    """Return (alpha_m, alpha_f, gamma, beta) for spectral radius rho_inf in [0, 1]."""
    rho = float(rho_inf)
    if not 0.0 <= rho <= 1.0:
        raise ValueError("rho_inf must lie in [0, 1]")
    am = (2.0 * rho - 1.0) / (rho + 1.0)
    af = rho / (rho + 1.0)
    gamma = 0.5 - am + af
    beta = 0.25 * (1.0 - am + af) ** 2
    return am, af, gamma, beta


def genalpha_amplification(omega_h, rho_inf, xi=0.0):
    """3x3 amplification matrix of the linear SDOF x'' + 2 xi w x' + w^2 x = 0.

    State s = (x_n, h v_n, h^2 a_n); s_{n+1} = G s_n.
    """
    am, af, gam, bet = genalpha_params(rho_inf)
    z = float(omega_h)
    zeta = 2.0 * float(xi)
    pX = np.array([1.0, 1.0, 0.5 - bet])   # X' = pX.s + bet*A'
    pV = np.array([0.0, 1.0, 1.0 - gam])   # V' = pV.s + gam*A'
    eA = np.array([0.0, 0.0, 1.0])
    eV = np.array([0.0, 1.0, 0.0])
    eX = np.array([1.0, 0.0, 0.0])
    # balance: A_{am} + zeta z V_{af} + z^2 X_{af} = 0 with (.)_{af} = (1-af)(.)' + af(.)_n
    c = (1.0 - am) + zeta * z * (1.0 - af) * gam + z * z * (1.0 - af) * bet
    r = am * eA + zeta * z * ((1.0 - af) * pV + af * eV) + z * z * ((1.0 - af) * pX + af * eX)
    Arow = -r / c
    G = np.vstack([pX + bet * Arow, pV + gam * Arow, Arow])
    return G


class GeneralizedAlpha:
    """Generalized-alpha integrator with Newton iterations on the acceleration.

    Parameters
    ----------
    M : array (n, n) or callable M(x) -> array
    f_int : callable f_int(x, v, t) -> (n,)
    f_ext : callable f_ext(t) -> (n,) or None (zero)
    rho_inf : high-frequency spectral radius in [0, 1]
    K, Cd : optional callables K(x, v, t), Cd(x, v, t) -> (n, n) tangents, or
        constant arrays. If omitted they are formed by forward differences.
    """

    def __init__(self, M, f_int, f_ext, rho_inf, K=None, Cd=None,
                 fd_eps=1e-7, tol=1e-10, max_iter=20):
        self.M = M
        self.f_int = f_int
        self.f_ext = f_ext
        self.K = K
        self.Cd = Cd
        self.fd_eps = float(fd_eps)
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.am, self.af, self.gamma, self.beta = genalpha_params(rho_inf)

    def _M(self, x):
        M = self.M(x) if callable(self.M) else self.M
        return np.atleast_2d(np.asarray(M, dtype=float))

    def _fint(self, x, v, t):
        return np.atleast_1d(np.asarray(self.f_int(x, v, t), dtype=float)).ravel()

    def _fext(self, t, n):
        if self.f_ext is None:
            return np.zeros(n)
        return np.atleast_1d(np.asarray(self.f_ext(t), dtype=float)).ravel()

    def _tangent(self, given, which, x, v, t, f0):
        n = x.size
        if given is not None:
            Mat = given(x, v, t) if callable(given) else given
            return np.atleast_2d(np.asarray(Mat, dtype=float))
        T = np.zeros((n, n))
        for j in range(n):
            eps = self.fd_eps * max(1.0, abs(x[j] if which == 'x' else v[j]))
            xp = x.copy(); vp = v.copy()
            if which == 'x':
                xp[j] += eps
            else:
                vp[j] += eps
            T[:, j] = (self._fint(xp, vp, t) - f0) / eps
        return T

    def _initial_accel(self, x0, v0, t0):
        n = x0.size
        return np.linalg.solve(self._M(x0),
                               self._fext(t0, n) - self._fint(x0, v0, t0))

    def _step(self, xn, vn, an, tn, h):
        am, af, gam, bet = self.am, self.af, self.gamma, self.beta
        n = xn.size
        t1 = tn + h
        t_af = (1.0 - af) * t1 + af * tn
        a = an.copy()

        def kin(a_):
            x1 = xn + h * vn + h * h * ((0.5 - bet) * an + bet * a_)
            v1 = vn + h * ((1.0 - gam) * an + gam * a_)
            return x1, v1

        def residual(a_):
            x1, v1 = kin(a_)
            a_am = (1.0 - am) * a_ + am * an
            x_af = (1.0 - af) * x1 + af * xn
            v_af = (1.0 - af) * v1 + af * vn
            R = self._M(x_af) @ a_am + self._fint(x_af, v_af, t_af) - self._fext(t_af, n)
            return R, x_af, v_af

        scale = max(1.0, float(np.max(np.abs(self._fext(t_af, n)))))
        for _ in range(self.max_iter):
            R, x_af, v_af = residual(a)
            if np.max(np.abs(R)) <= self.tol * scale:
                x1, v1 = kin(a)
                return x1, v1, a
            Kt = self._tangent(self.K, 'x', x_af, v_af, t_af, self._fint(x_af, v_af, t_af))
            Ct = self._tangent(self.Cd, 'v', x_af, v_af, t_af, self._fint(x_af, v_af, t_af))
            Jac = (1.0 - am) * self._M(x_af) + (1.0 - af) * (gam * h * Ct + bet * h * h * Kt)
            a = a + np.linalg.solve(Jac, -R)
        R, _, _ = residual(a)
        if np.max(np.abs(R)) <= self.tol * scale:
            x1, v1 = kin(a)
            return x1, v1, a
        raise RuntimeError("generalized-alpha Newton iteration did not converge")

    def run(self, x0, v0, t0, h, nsteps):
        """Integrate nsteps steps. Returns (t, X, V, A), each with nsteps+1 rows."""
        x = np.atleast_1d(np.asarray(x0, dtype=float)).ravel()
        v = np.atleast_1d(np.asarray(v0, dtype=float)).ravel()
        n = x.size
        T = t0 + h * np.arange(nsteps + 1)
        X = np.empty((nsteps + 1, n))
        V = np.empty((nsteps + 1, n))
        A = np.empty((nsteps + 1, n))
        X[0], V[0] = x, v
        A[0] = self._initial_accel(x, v, t0)
        for k in range(nsteps):
            x, v, a = self._step(X[k], V[k], A[k], T[k], h)
            X[k + 1], V[k + 1], A[k + 1] = x, v, a
        return T, X, V, A


def velocity_verlet(accel, x0, v0, t0, h, nsteps):
    """Velocity Verlet with lagged velocity in the second acceleration.

    accel(x, v, t) -> (n,). Returns (t, X, V), each with nsteps+1 rows.
    """
    x = np.atleast_1d(np.asarray(x0, dtype=float)).ravel()
    v = np.atleast_1d(np.asarray(v0, dtype=float)).ravel()
    T = t0 + h * np.arange(nsteps + 1)
    X = np.empty((nsteps + 1, x.size))
    V = np.empty((nsteps + 1, x.size))
    X[0], V[0] = x, v
    for k in range(nsteps):
        t = T[k]
        vh = v + 0.5 * h * np.asarray(accel(x, v, t), dtype=float)
        x = x + h * vh
        vhalf = vh + 0.5 * h * np.asarray(accel(x, vh, T[k + 1]), dtype=float)
        v = vhalf
        X[k + 1], V[k + 1] = x, v
    return T, X, V


def symplectic_euler(accel, x0, v0, t0, h, nsteps):
    """Semi-implicit Euler: v_{n+1} = v_n + h a(x_n, v_n, t_n); x_{n+1} = x_n + h v_{n+1}."""
    x = np.atleast_1d(np.asarray(x0, dtype=float)).ravel()
    v = np.atleast_1d(np.asarray(v0, dtype=float)).ravel()
    T = t0 + h * np.arange(nsteps + 1)
    X = np.empty((nsteps + 1, x.size))
    V = np.empty((nsteps + 1, x.size))
    X[0], V[0] = x, v
    for k in range(nsteps):
        v = v + h * np.asarray(accel(x, v, T[k]), dtype=float)
        x = x + h * v
        X[k + 1], V[k + 1] = x, v
    return T, X, V
