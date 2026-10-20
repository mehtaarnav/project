"""Timoshenko shaft on spring bearings (docs/theory/07_shafts.md, sections 2-3, 6).

Each node has six DOFs ordered as
    0: u_x   1: u_y   2: u_z   3: psi_x   4: psi_y   5: phi_z
with psi_x ~ du_x/dz (rotation about +y, conjugate to My) and
psi_y ~ du_y/dz (rotation about -x, conjugate to -Mx).

Loads are tuples (z, Fx, Fy, Fz, Mx, My, Tz) applied at a node.
Bearing reactions are the forces the bearing applies on the shaft (-k u).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

_UX, _UY, _UZ, _PX, _PY, _PZ = range(6)
_TOL = 1e-12


def cowper_kappa(nu: float) -> float:
    """Shear coefficient of a solid circular section (Cowper 1966)."""
    return 6.0 * (1.0 + nu) / (7.0 + 6.0 * nu)


@dataclass
class ShaftResult:
    z: np.ndarray                 # node stations (n,)
    diam: np.ndarray              # element diameters (n-1,)
    u: np.ndarray                 # nodal DOFs (n, 6)
    bearing_z: np.ndarray         # (nb,)
    reactions: np.ndarray         # (nb, 4): Fx, Fy, Fz, Tz applied by bearing ON shaft
    loads: np.ndarray             # (nl, 7): z, Fx, Fy, Fz, Mx, My, Tz (applied)
    _all: np.ndarray = field(repr=False, default=None)  # (nl+nb, 7) for equilibrium

    def _elem_index(self, z: float) -> int:
        k = int(np.searchsorted(self.z, z, side="left")) - 1
        return min(max(k, 0), len(self.diam) - 1)

    def internal(self, z: float):
        """Stress resultants just to the left of station z (equilibrium of z' < z).

        Returns (Vx, Vy, Mx, My, T, N); N is tension positive.
        """
        A = self._all
        m = A[:, 0] < z - _TOL
        d = A[m, 0] - z
        Fx, Fy, Fz = A[m, 1].sum(), A[m, 2].sum(), A[m, 3].sum()
        Mx = np.sum(A[m, 4] - d * A[m, 2])
        My = np.sum(A[m, 5] + d * A[m, 1])
        T = A[m, 6].sum()
        return (-Fx, -Fy, -Mx, -My, -T, -Fz)

    def stress(self, z: float, Ktb: float = 1.0, Kts: float = 1.0, Kta: float = 1.0):
        """Nominal and peak stresses at station z.

        Returns (sig_b, tau, sig_a, sig_vm): nominal bending and torsional
        shear magnitudes, signed nominal axial stress, and the peak von Mises
        stress sqrt((Ktb*sig_b + Kta*sig_a)^2 + 3 (Kts*tau)^2) (A7.3).
        """
        d = self.diam[self._elem_index(z)]
        _, _, Mx, My, T, N = self.internal(z)
        M = math.hypot(Mx, My)
        sig_b = 32.0 * M / (math.pi * d**3)
        tau = 16.0 * abs(T) / (math.pi * d**3)
        sig_a = 4.0 * N / (math.pi * d**2)
        sig_vm = math.sqrt((Ktb * sig_b + Kta * sig_a) ** 2 + 3.0 * (Kts * tau) ** 2)
        return sig_b, tau, sig_a, sig_vm


class Shaft:
    """Stepped Timoshenko shaft. `stations_z` (n) with `diameters` (n-1) per element."""

    torsion_ground_k: float = 1e12

    def __init__(self, stations_z, diameters, E: float, nu: float, kappa: float | None = None):
        z = np.asarray(stations_z, dtype=float).ravel()
        d = np.asarray(diameters, dtype=float).ravel()
        if z.size < 2 or np.any(np.diff(z) <= 0):
            raise ValueError("stations_z must be strictly increasing with >= 2 entries")
        if d.size != z.size - 1 or np.any(d <= 0):
            raise ValueError("need one positive diameter per element")
        self.z = z
        self.d = d
        self.E = float(E)
        self.nu = float(nu)
        self.G = self.E / (2.0 * (1.0 + self.nu))
        self.kappa = cowper_kappa(self.nu) if kappa is None else float(kappa)
        self.bearings: list[tuple[float, float, float, float]] = []

    # ----- model building -------------------------------------------------
    def add_bearing(self, z: float, k_r: float = 1e12, k_a: float = 0.0, k_t: float = 0.0) -> None:
        if not (self.z[0] - _TOL <= z <= self.z[-1] + _TOL):
            raise ValueError("bearing outside shaft")
        self.bearings.append((float(z), float(k_r), float(k_a), float(k_t)))

    def _mesh(self, extra_z) -> tuple[np.ndarray, np.ndarray]:
        zs = np.unique(np.concatenate([self.z, np.asarray(extra_z, float).ravel()]))
        keep = np.concatenate([[True], np.diff(zs) > _TOL])
        zs = zs[keep]
        mid = 0.5 * (zs[:-1] + zs[1:])
        idx = np.clip(np.searchsorted(self.z, mid) - 1, 0, len(self.d) - 1)
        return zs, self.d[idx]

    def _node(self, zs: np.ndarray, z: float) -> int:
        k = int(np.argmin(np.abs(zs - z)))
        if abs(zs[k] - z) > 1e-9 * max(1.0, abs(z)):
            raise ValueError(f"station {z} not in mesh")
        return k

    # ----- solve -----------------------------------------------------------
    def solve(self, loads=(), torque_ground_z: float | None = None) -> ShaftResult:
        if not self.bearings:
            raise ValueError("at least one bearing is required")
        loads = [tuple(float(v) for v in ld) + (0.0,) * (7 - len(ld)) for ld in loads]
        loads_arr = np.array(loads, dtype=float).reshape(-1, 7)
        brg = [list(b) for b in self.bearings]

        # Rigid axial and torsional grounds (A7.4 / solver convention): if no
        # thrust (k_a) or torsion (k_t) stiffness is given anywhere, a rigid
        # ground is placed at the torque-reaction bearing (default: first added).
        gz = self.bearings[0][0] if torque_ground_z is None else float(torque_ground_z)
        if not any(b[2] > 0 for b in brg):
            for b in brg:
                if abs(b[0] - gz) < _TOL:
                    b[2] = self.torsion_ground_k
                    break
            else:
                raise ValueError("torque_ground_z must match a bearing")
        if not any(b[3] > 0 for b in brg):
            for b in brg:
                if abs(b[0] - gz) < _TOL:
                    b[3] = self.torsion_ground_k
                    break
            else:
                raise ValueError("torque_ground_z must match a bearing")

        bz = np.array([b[0] for b in brg])
        zs, de = self._mesh(np.concatenate([bz, loads_arr[:, 0]]))
        n = zs.size
        ndof = 6 * n
        K = np.zeros((ndof, ndof))
        E, G, kap = self.E, self.G, self.kappa
        for e in range(n - 1):
            L = zs[e + 1] - zs[e]
            dd = de[e]
            A = math.pi * dd**2 / 4.0
            I = math.pi * dd**4 / 64.0
            J = 2.0 * I
            EI, EA, GJ, kGA = E * I, E * A, G * J, kap * G * A
            Phi = 12.0 * EI / (kGA * L * L)
            c = EI / ((1.0 + Phi) * L**3)
            kb = c * np.array([
                [12.0, 6 * L, -12.0, 6 * L],
                [6 * L, (4 + Phi) * L * L, -6 * L, (2 - Phi) * L * L],
                [-12.0, -6 * L, 12.0, -6 * L],
                [6 * L, (2 - Phi) * L * L, -6 * L, (4 + Phi) * L * L],
            ])
            a, b = e, e + 1
            dx = [6 * a + _UX, 6 * a + _PX, 6 * b + _UX, 6 * b + _PX]
            dy = [6 * a + _UY, 6 * a + _PY, 6 * b + _UY, 6 * b + _PY]
            K[np.ix_(dx, dx)] += kb
            K[np.ix_(dy, dy)] += kb
            ka = EA / L * np.array([[1.0, -1.0], [-1.0, 1.0]])
            dz = [6 * a + _UZ, 6 * b + _UZ]
            K[np.ix_(dz, dz)] += ka
            kt = GJ / L * np.array([[1.0, -1.0], [-1.0, 1.0]])
            dp = [6 * a + _PZ, 6 * b + _PZ]
            K[np.ix_(dp, dp)] += kt

        F = np.zeros(ndof)
        for ld in loads_arr:
            k = self._node(zs, ld[0])
            F[6 * k + _UX] += ld[1]
            F[6 * k + _UY] += ld[2]
            F[6 * k + _UZ] += ld[3]
            F[6 * k + _PX] += ld[5]   # My is conjugate to psi_x
            F[6 * k + _PY] -= ld[4]   # Mx is conjugate to -psi_y
            F[6 * k + _PZ] += ld[6]

        bnodes = []
        for (z, k_r, k_a, k_t) in brg:
            k = self._node(zs, z)
            bnodes.append(k)
            K[6 * k + _UX, 6 * k + _UX] += k_r
            K[6 * k + _UY, 6 * k + _UY] += k_r
            K[6 * k + _UZ, 6 * k + _UZ] += k_a
            K[6 * k + _PZ, 6 * k + _PZ] += k_t

        u = np.linalg.solve(K, F).reshape(n, 6)
        reactions = np.zeros((len(brg), 4))
        for i, (z, k_r, k_a, k_t) in enumerate(brg):
            k = bnodes[i]
            reactions[i] = (-k_r * u[k, _UX], -k_r * u[k, _UY], -k_a * u[k, _UZ], -k_t * u[k, _PZ])

        allrows = [loads_arr]
        br_arr = np.zeros((len(brg), 7))
        br_arr[:, 0] = bz
        br_arr[:, 1:4] = reactions[:, :3]
        br_arr[:, 6] = reactions[:, 3]
        allrows.append(br_arr)
        all_arr = np.vstack(allrows)

        return ShaftResult(z=zs, diam=de, u=u, bearing_z=bz, reactions=reactions,
                           loads=loads_arr, _all=all_arr)


def influence_matrix(shaft: Shaft, load_z, eval_z, unit_components=(1, 1, 1, 1, 1, 1)):
    """Linear response to unit loads (for superposition).

    Columns: each load station in `load_z` times each active component
    (Fx, Fy, Fz, Mx, My, Tz) with unit magnitude. Rows: bearing reactions
    (Fx, Fy, Fz, Tz per bearing) followed by internal (Vx, Vy, Mx, My, T, N)
    at each station of `eval_z`. Bearings must be added before calling.
    """
    cols = []
    names = []
    for zl in np.atleast_1d(load_z):
        for c in range(6):
            if not unit_components[c]:
                continue
            ld = [0.0] * 6
            ld[c] = 1.0
            names.append((zl, c))
            cols.append((zl, ld))
    out = []
    for zl, ld in cols:
        res = shaft.solve([(zl, *ld)])
        vec = [*res.reactions.ravel()]
        for ze in np.atleast_1d(eval_z):
            vec.extend(res.internal(ze))
        out.append(vec)
    return np.array(out).T


def first_bending_frequency(shaft: Shaft, mass_per_length: float, lumped=()) -> float:
    """Lowest lateral (x-plane) natural frequency in Hz, lumped-mass model.

    Translational masses are lumped at the nodes (half of each adjacent
    element's mass plus any `lumped` (z, m) point masses). Rotary inertia is
    neglected, so rotations are condensed out statically (A7.1 quasi-static
    check uses this value).
    """
    lumped = list(lumped)
    extra = [b[0] for b in shaft.bearings] + [m_[0] for m_ in lumped]
    zs, de = shaft._mesh(extra)
    n = zs.size
    # assemble bending stiffness of the x-plane and bearing springs
    K = np.zeros((2 * n, 2 * n))  # DOFs: (u_x, psi_x) per node
    for e in range(n - 1):
        L = zs[e + 1] - zs[e]
        dd = de[e]
        I = math.pi * dd**4 / 64.0
        EI = shaft.E * I
        kGA = shaft.kappa * shaft.G * math.pi * dd**2 / 4.0
        Phi = 12.0 * EI / (kGA * L * L)
        c = EI / ((1.0 + Phi) * L**3)
        kb = c * np.array([
            [12.0, 6 * L, -12.0, 6 * L],
            [6 * L, (4 + Phi) * L * L, -6 * L, (2 - Phi) * L * L],
            [-12.0, -6 * L, 12.0, -6 * L],
            [6 * L, (2 - Phi) * L * L, -6 * L, (4 + Phi) * L * L],
        ])
        idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(idx, idx)] += kb
    for (z, k_r, _, _) in shaft.bearings:
        k = int(np.argmin(np.abs(zs - z)))
        K[2 * k, 2 * k] += k_r
    mass = np.zeros(n)
    for e in range(n - 1):
        L = zs[e + 1] - zs[e]
        mass[e] += 0.5 * mass_per_length * L
        mass[e + 1] += 0.5 * mass_per_length * L
    for (z, m) in lumped:
        mass[int(np.argmin(np.abs(zs - z)))] += m
    # static condensation of psi (rotations) onto translational DOFs
    ut = np.arange(n) * 2
    pt = ut + 1
    Kuu = K[np.ix_(ut, ut)]
    Kup = K[np.ix_(ut, pt)]
    Kpp = K[np.ix_(pt, pt)]
    Kc = Kuu - Kup @ np.linalg.solve(Kpp, Kup.T)
    Mi = 1.0 / np.sqrt(mass)
    A = Mi[:, None] * Kc * Mi[None, :]
    w2 = np.linalg.eigvalsh(0.5 * (A + A.T))
    if w2[0] <= 0.0:
        raise ValueError("shaft is not restrained: add bearings")
    return math.sqrt(w2[0]) / (2.0 * math.pi)
