"""Single-cylinder firing pulses (docs/theory/11_engine_friction_uq.md, section 1).

Crank angle theta [rad], 0 = TDC firing, 4 pi cycle: expansion [0, pi], exhaust [pi, 2 pi],
intake [2 pi, 3 pi], compression [3 pi, 4 pi]. On the closed part (psi = theta, or theta - 4 pi,
in [-pi, pi]) the first law with polytropic exponent n, d(p V^n) = (n - 1) V^(n-1) dQ, gives

    p(psi) = [p_amb V1^n + (n - 1) Q_in I(psi)] / V(psi)^n,   I = int_{-pi}^{psi} V^(n-1) x_b' dpsi

which is linear in Q_in. Closing Q_in by the map (indicated work = 4 pi T_m) gives

    T_gas(theta, T_m) = g_mot(theta) + T_m g_Q(theta)

with g_mot = A s' (p_mot - p_amb) (motored, zero mean) and g_Q = 4 pi A s' p_Q / W_Q (unit mean),
so the cycle mean is T_m exactly (11 section 1.1). Both are tabulated on a uniform periodic grid
of `n` nodes and interpolated linearly in the jitted loop; the table means are made exact on the
interpolant (trapezoid = exact integral of a piecewise-linear function).
"""
import numpy as np

from cvtsim.core.jit import maybe_njit

FOUR_PI = 4.0 * np.pi


@maybe_njit
def crank_kinematics(theta, r, lam):
    """Slider-crank piston displacement from TDC and its theta-derivatives (s, s', s'').

    s = r(1 - cos) + l(1 - D), D = sqrt(1 - lam^2 sin^2), l = r/lam."""
    sn = np.sin(theta); cs = np.cos(theta)
    D = np.sqrt(1.0 - lam * lam * sn * sn)
    s = r * (1.0 - cs) + (r / lam) * (1.0 - D)
    sp = r * sn * (1.0 + lam * cs / D)
    spp = r * cs + r * lam * ((cs * cs - sn * sn) / D + lam * lam * sn * sn * cs * cs / D ** 3)
    return s, sp, spp


@maybe_njit
def rec_torque(theta, w, m_rec, r, lam):
    """Reciprocating velocity torque T_rec = -m_rec s' s'' w^2 and the inertia m_rec s'^2.

    With J(theta) = J_f + m_rec s'^2 the crank equation J w' = T + T_rec conserves
    1/2 J w^2 exactly (11 section 1.2)."""
    _, sp, spp = crank_kinematics(theta, r, lam)
    return -m_rec * sp * spp * w * w, m_rec * sp * sp


@maybe_njit
def table_eval(tab, theta):
    """Periodic linear interpolation of tab (n nodes over [0, 4 pi)) at theta."""
    n = tab.shape[0]
    q = theta / FOUR_PI
    q = (q - np.floor(q)) * n
    i = int(q)
    if i >= n:
        i = 0
        q = 0.0
    f = q - i
    j = i + 1 if i + 1 < n else 0
    return tab[i] + f * (tab[j] - tab[i])


@maybe_njit
def gas_torque(theta, T_m, gm, gq):
    """Gas torque T_gas = g_mot(theta) + T_m g_Q(theta) from the tables."""
    return table_eval(gm, theta) + T_m * table_eval(gq, theta)


def wiebe(psi, theta_s, dtheta, a, m):
    """Burned fraction x_b and dx_b/dpsi (zero before theta_s)."""
    t = np.clip((psi - theta_s) / dtheta, 0.0, None)
    xb = 1.0 - np.exp(-a * t ** (m + 1))
    dxb = a * (m + 1) * t ** m * np.exp(-a * t ** (m + 1)) / dtheta
    return xb, dxb


def engine_geometry(e):
    """(r, lam, A, V_c, V1) from the engine params dict (bore, stroke, rod_length, r_c)."""
    r = 0.5 * e["stroke"]
    A = 0.25 * np.pi * e["bore"] ** 2
    V_d = A * e["stroke"]
    V_c = V_d / (e["compression_ratio"] - 1.0)
    return r, r / e["rod_length"], A, V_c, V_c + V_d


def _wiebe_args(e):
    return (np.radians(e["wiebe_theta_s_deg"]), np.radians(e["wiebe_duration_deg"]),
            e["wiebe_a"], e["wiebe_m"])


def work_per_heat(e):
    """Continuum indicated work per unit heat W_Q = int x_b' [1 - (V/V1)^(n-1)] dpsi over
    [theta_s, pi] (first law over the closed part; adaptive quadrature)."""
    from scipy.integrate import quad
    r, lam, A, V_c, V1 = engine_geometry(e)
    n = e["n_poly"]
    ts, dt, a, m = _wiebe_args(e)

    def f(psi):
        V = V_c + A * crank_kinematics(psi, r, lam)[0]
        return wiebe(psi, ts, dt, a, m)[1] * (1.0 - (V / V1) ** (n - 1.0))
    pts = [ts + k * dt for k in (0.5, 1.0, 1.5, 2.0) if ts + k * dt < np.pi]
    return quad(f, ts, np.pi, points=pts, limit=400, epsabs=0.0, epsrel=1e-13)[0]


def pulse_table(e, n=None):
    """Tables (gm, gq) on n nodes theta_j = 4 pi j/n, plus diagnostics.

    n must be a multiple of 48 so that nodes fall on BDC (+-pi) and on theta_s (multiples of
    15 deg); the cell integrals of I use 8-point Gauss-Legendre. Returns dict(gm, gq, theta,
    mean_raw (trapezoid mean of A s' p_Q per unit heat), W_Q, n)."""
    n = int(n or e["pulse_table_n"])
    if n % 48:
        raise ValueError("pulse_table_n must be a multiple of 48")
    r, lam, A, V_c, V1 = engine_geometry(e)
    npoly, p_amb = e["n_poly"], e["p_amb"]
    ts, dt, a, m = _wiebe_args(e)
    theta = FOUR_PI * np.arange(n) / n
    # closed part psi in [-pi, pi]: nodes j = 3n/4 .. n-1 (psi = theta - 4 pi) and 0 .. n/4
    psi = np.concatenate([theta[3 * n // 4:] - FOUR_PI, theta[:n // 4 + 1]])
    xg, wg = np.polynomial.legendre.leggauss(8)
    h = psi[1] - psi[0]
    q = (psi[:-1, None] + 0.5 * h * (xg[None, :] + 1.0))
    Vq = V_c + A * crank_kinematics(q, r, lam)[0]
    cell = 0.5 * h * np.sum(wg[None, :] * Vq ** (npoly - 1.0) * wiebe(q, ts, dt, a, m)[1], axis=1)
    I = np.concatenate([[0.0], np.cumsum(cell)])
    s, sp, _ = crank_kinematics(psi, r, lam)
    V = V_c + A * s
    pQ = (npoly - 1.0) * I / V ** npoly
    pm = p_amb * ((V1 / V) ** npoly - 1.0)
    gQ_c, gm_c = A * sp * pQ, A * sp * pm
    gQ = np.zeros(n); gm = np.zeros(n)
    k = n // 4
    gQ[3 * n // 4:] = gQ_c[:k]; gQ[:k + 1] = gQ_c[k:]
    gm[3 * n // 4:] = gm_c[:k]; gm[:k + 1] = gm_c[k:]
    mean_raw = gQ.mean()          # trapezoid mean over the periodic grid
    return dict(gm=gm - gm.mean(), gq=gQ / mean_raw, theta=theta, mean_raw=mean_raw,
                W_Q=work_per_heat(e), n=n, A=A, V1=V1, pQ=pQ, pm=pm, psi=psi)
