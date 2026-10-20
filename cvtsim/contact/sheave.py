"""Analytic cone-sheave V-belt contact (docs/theory/04_contact.md).

Penalty normal force (compressive only) on each flank, regularized Coulomb
friction in the X-Y plane, no mesh contact search. Pure NumPy kernels,
numba-compatible via maybe_njit.

Frame: positive rotation about +Z. Forces are those ON THE BELT NODES unless
stated; sheave torques are the torques ON the sheave half about +Z.
"""
import numpy as np

from cvtsim.core.jit import maybe_njit


@maybe_njit
def capstan_ratio(mu, theta, beta):
    """Tight/slack tension ratio at full slip in a V-groove: exp(mu*theta/sin(beta)) (C4)."""
    return np.exp(mu * theta / np.sin(beta))


@maybe_njit
def axial_force_incipient(Tq, R, beta, mu):
    """Axial thrust per flank at full slip: F_a = cos(beta)*Tq/(2*mu*R) (C5)."""
    return np.cos(beta) * Tq / (2.0 * mu * R)


@maybe_njit
def sheave_contact(x, v, cx, cy, r_star, r_star_dot, beta, k_n, c_n, mu, v_eps,
                   r_outer, w_fixed, w_mov):
    """Contact of belt nodes with one V-sheave pair (fixed + movable flank).

    Returns (f (N,2) force on nodes, Fa thrust on movable half (N),
             tq_fixed, tq_mov (N m, torque on each half about +Z),
             P_diss (W, friction only, >= 0), U_pen (J),
             delta (N,) penetration, N_flank (N,) per-flank normal force).
    Contact only where r < r_star and r < r_outer (C2); zero force elsewhere.
    """
    N = x.shape[0]
    f = np.zeros((N, 2))
    delta = np.zeros(N)
    Nf = np.zeros(N)
    Fa = 0.0
    tq_f = 0.0
    tq_m = 0.0
    Pd = 0.0
    Up = 0.0
    s_b = np.sin(beta)
    c_b = np.cos(beta)
    for k in range(N):
        dx = x[k, 0] - cx
        dy = x[k, 1] - cy
        r = np.sqrt(dx * dx + dy * dy)
        if r < r_star and r < r_outer and r > 0.0:
            ex = dx / r
            ey = dy / r
            d = (r_star - r) * s_b
            rdot = ex * v[k, 0] + ey * v[k, 1]
            ddot = (r_star_dot - rdot) * s_b
            Nn = k_n * d + c_n * ddot
            if Nn < 0.0:
                Nn = 0.0
            delta[k] = d
            Nf[k] = Nn
            Up += k_n * d * d  # two flanks, each contributing 1/2 k_n d^2
            Fa += Nn * c_b
            fr = 2.0 * Nn * s_b
            fx = fr * ex
            fy = fr * ey
            for flank in range(2):
                if flank == 0:
                    w = w_fixed
                else:
                    w = w_mov
                ux = -w * dy  # w zhat x (x - c)
                uy = w * dx
                sx = v[k, 0] - ux
                sy = v[k, 1] - uy
                den = np.sqrt(sx * sx + sy * sy + v_eps * v_eps)
                gx = -mu * Nn * sx / den  # friction ON the node
                gy = -mu * Nn * sy / den
                fx += gx
                fy += gy
                Pd += -(gx * sx + gy * sy)
                tau = -(dx * gy - dy * gx)  # torque on sheave half: (x-c) x (-g)
                if flank == 0:
                    tq_f += tau
                else:
                    tq_m += tau
            f[k, 0] = fx
            f[k, 1] = fy
    return f, Fa, tq_f, tq_m, Pd, Up, delta, Nf
