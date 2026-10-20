"""Segmented belt kernels (docs/theory/03_belt.md).

Plain-NumPy, array-in/array-out, numba-compatible (maybe_njit). Frames: belt
nodes move in the X-Y plane; positive rotation about +Z.

Conventions:
  * x, v: (N, 2) node positions / velocities.
  * L0: scalar rest segment length (lumped mass m = rhoA*L0).
  * closed=True: N elements (k, k+1 mod N), bending at every node.
    closed=False: N-1 elements, bending at interior nodes 1..N-2 only.
  * Returned "U_strain" is the total elastic energy (axial + bending).
"""
import numpy as np

from cvtsim.core.jit import maybe_njit


@maybe_njit
def belt_forces(x, v, L0, EA, eta, EI, cb, closed):
    """Internal belt forces on nodes.

    Returns (f (N,2) force on each node, U_strain (J), P_diss (W)).
    Axial: T = EA*(eps + eta*eps_dot). Bending: M = (EI/L0)*phi + (cb/L0)*phi_dot,
    force F = -M * grad(phi).
    """
    N = x.shape[0]
    f = np.zeros((N, 2))
    U = 0.0
    P = 0.0
    if closed:
        ne = N
    else:
        ne = N - 1
    # axial elements
    for e in range(ne):
        k = e
        k1 = (e + 1) % N
        dx = x[k1, 0] - x[k, 0]
        dy = x[k1, 1] - x[k, 1]
        l = np.sqrt(dx * dx + dy * dy)
        tx = dx / l
        ty = dy / l
        eps = (l - L0) / L0
        epsd = (tx * (v[k1, 0] - v[k, 0]) + ty * (v[k1, 1] - v[k, 1])) / L0
        T = EA * (eps + eta * epsd)
        U += 0.5 * EA * L0 * eps * eps
        P += EA * eta * L0 * epsd * epsd
        f[k1, 0] -= T * tx
        f[k1, 1] -= T * ty
        f[k, 0] += T * tx
        f[k, 1] += T * ty
    # bending at nodes
    if closed:
        start = 0
        stop = N
    else:
        start = 1
        stop = N - 1
    for k in range(start, stop):
        km = (k - 1) % N
        kp = (k + 1) % N
        ex = x[k, 0] - x[km, 0]
        ey = x[k, 1] - x[km, 1]
        fx = x[kp, 0] - x[k, 0]
        fy = x[kp, 1] - x[k, 1]
        e2 = ex * ex + ey * ey
        f2 = fx * fx + fy * fy
        cross = ex * fy - ey * fx
        dot = ex * fx + ey * fy
        phi = np.arctan2(cross, dot)
        # gradients (zhat x u = (-u_y, u_x))
        gp_x = -fy / f2  # d phi / d x_{k+1}
        gp_y = fx / f2
        gm_x = -ey / e2  # d phi / d x_{k-1}  (zhat x e)/|e|^2 = (-e_y, e_x)/|e|^2
        gm_y = ex / e2
        g_x = -(gp_x + gm_x)  # d phi / d x_k
        g_y = -(gp_y + gm_y)
        # phi_dot = grad(phi).v. Since grad_{k-1}+grad_k+grad_{k+1} = 0, this equals
        # the relative-velocity form below, which avoids cancellation of large terms.
        phid = (gm_x * (v[km, 0] - v[k, 0]) + gm_y * (v[km, 1] - v[k, 1])
                + gp_x * (v[kp, 0] - v[k, 0]) + gp_y * (v[kp, 1] - v[k, 1]))
        M = (EI / L0) * phi + (cb / L0) * phid
        U += 0.5 * (EI / L0) * phi * phi
        P += (cb / L0) * phid * phid
        f[km, 0] -= M * gm_x
        f[km, 1] -= M * gm_y
        f[k, 0] -= M * g_x
        f[k, 1] -= M * g_y
        f[kp, 0] -= M * gp_x
        f[kp, 1] -= M * gp_y
    return f, U, P


@maybe_njit
def belt_tensions(x, v, L0, EA, eta, closed):
    """Element tensions T = EA*(eps + eta*eps_dot), shape (Ne,), N."""
    N = x.shape[0]
    if closed:
        ne = N
    else:
        ne = N - 1
    T = np.zeros(ne)
    for e in range(ne):
        k = e
        k1 = (e + 1) % N
        dx = x[k1, 0] - x[k, 0]
        dy = x[k1, 1] - x[k, 1]
        l = np.sqrt(dx * dx + dy * dy)
        tx = dx / l
        ty = dy / l
        eps = (l - L0) / L0
        epsd = (tx * (v[k1, 0] - v[k, 0]) + ty * (v[k1, 1] - v[k, 1])) / L0
        T[e] = EA * (eps + eta * epsd)
    return T


@maybe_njit
def belt_length(r1, r2, C):
    """Exact open-belt (external tangent) path length around circles
    (0,0,r1) and (C,0,r2): 2 C cos(phi) + pi (r1+r2) + 2 phi (r2-r1),
    with sin(phi) = (r2-r1)/C (docs/theory/05_cvt.md section 2)."""
    phi = np.arcsin((r2 - r1) / C)
    return 2.0 * C * np.cos(phi) + np.pi * (r1 + r2) + 2.0 * phi * (r2 - r1)


@maybe_njit
def init_belt_on_path(rp, rs, C, N):
    """N nodes equally spaced by arc length along the open belt path around
    circles (0,0,rp) and (C,0,rs), counter-clockwise (belt moves +Z sense).
    Path order: bottom tangent (+X), arc on secondary, top tangent (-X), arc on
    primary. Returns (x (N,2), L0) with L0 = Lpath/N (closed loop)."""
    phi = np.arcsin((rs - rp) / C)
    s = np.sin(phi)
    c = np.cos(phi)
    a = C * c  # tangent length
    Lp = belt_length(rp, rs, C)
    # segment lengths and boundaries
    l0 = a
    l1 = rs * (np.pi + 2.0 * phi)
    l2 = a
    l3 = rp * (np.pi - 2.0 * phi)
    b0 = l0
    b1 = b0 + l1
    b2 = b1 + l2
    x = np.zeros((N, 2))
    L0 = Lp / N
    for k in range(N):
        sk = k * L0
        if sk > Lp:
            sk = Lp
        if sk <= b0:
            t = sk / a
            px1 = -rp * s
            py1 = -rp * c
            px2 = C - rs * s
            py2 = -rs * c
            x[k, 0] = px1 + t * (px2 - px1)
            x[k, 1] = py1 + t * (py2 - py1)
        elif sk <= b1:
            th = (-np.pi / 2.0 - phi) + (sk - b0) / rs
            x[k, 0] = C + rs * np.cos(th)
            x[k, 1] = rs * np.sin(th)
        elif sk <= b2:
            t = (sk - b1) / a
            px2 = C - rs * s
            py2 = rs * c
            px1 = -rp * s
            py1 = rp * c
            x[k, 0] = px2 + t * (px1 - px2)
            x[k, 1] = py2 + t * (py1 - py2)
        else:
            th = (np.pi / 2.0 + phi) + (sk - b2) / rp
            x[k, 0] = rp * np.cos(th)
            x[k, 1] = rp * np.sin(th)
    return x, L0


@maybe_njit
def node_mass(rhoA, L0):
    """Lumped node mass m = rhoA * L0 (kg)."""
    return rhoA * L0


@maybe_njit
def flank_stiffness(G, h_c, L0, b):
    """Per-flank per-node normal stiffness k_n = E_t h_c L0 / (b/2), E_t = 3G (B1)."""
    Et = 3.0 * G
    return Et * h_c * L0 / (0.5 * b)
