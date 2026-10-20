"""Primary ramp and secondary helix profiles (docs/theory/09_tuning_vehicle_belt.md, sections 1-2).

Both profiles are piecewise polynomials: knots `xk` (n+1,) and per-segment ascending
coefficients `c` (n, 7) in `s = x - xk[i]` (degree <= 6). Outside `[xk[0], xk[-1]]` the
function continues linearly with the end value and slope (zero curvature), so it stays C1.
`pp_eval` returns (f, f', f'') exactly inside the jitted loop: the profile is evaluated,
not resampled, so there is no table-resolution error to converge.

* Ramp `r_w(u)`: straight (`r_w0 + kappa u`), a polynomial of degree <= 6, or a table
  interpolated with a monotone cubic (PCHIP). Then `kappa = r_w'` and `kappa' = r_w''`.
* Helix `h(x)`: angle table `alpha(x)` at cam radius `r_c`, `h' = 1/(r_c tan alpha)`.
  `steps`: alpha piecewise constant; `h'` changes over a length `blend` centred on each
  break with the smootherstep `10t^3 - 15t^4 + 6t^5`, so `h''` and `h'''` are continuous
  (a sharp break is a Dirac in `h''`, an impulse the explicit integrator cannot carry).
  The blend is symmetric, so `h(x)` outside the blend zones equals the integral of the
  sharp step exactly. `linear`: `h'` linear between table points (`h''` jumps at the points).

Where the second derivative jumps (PCHIP knots, `linear` helix points, the end of a
profile) the equations of motion have a force discontinuity: an integrator crossing it
makes an O(h) energy error per crossing, negligible at the launch time step.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator

from cvtsim.core.jit import maybe_njit

NC = 7  # coefficients per segment (degree <= 6)


@maybe_njit
def pp_eval(xk, c, x):
    """(f, f', f'') of the piecewise polynomial (xk, c) at x, linear outside the knots."""
    n = c.shape[0]
    if x <= xk[0]:
        s = x - xk[0]
        return c[0, 0] + c[0, 1] * s, c[0, 1], 0.0
    i = 0
    if x >= xk[n]:
        i = n - 1
    else:
        while x >= xk[i + 1]:
            i += 1
    s = min(x, xk[n]) - xk[i]
    K = c.shape[1]
    f = c[i, K - 1]
    d = (K - 1) * c[i, K - 1]
    dd = (K - 1) * (K - 2) * c[i, K - 1]
    for j in range(K - 2, -1, -1):
        f = f * s + c[i, j]
        if j >= 1:
            d = d * s + j * c[i, j]
        if j >= 2:
            dd = dd * s + j * (j - 1) * c[i, j]
    if x >= xk[n]:
        return f + d * (x - xk[n]), d, 0.0
    return f, d, dd


class _Profile:
    def __init__(self, xk, c, kind):
        self.xk = np.ascontiguousarray(xk, dtype=np.float64)
        self.c = np.ascontiguousarray(c, dtype=np.float64)
        self.kind = kind
        if self.c.shape != (len(self.xk) - 1, NC) or np.any(np.diff(self.xk) <= 0.0):
            raise ValueError("profile knots must be strictly increasing with one coefficient row per segment")

    def __call__(self, x):
        """(f, f', f'') at scalar x."""
        return pp_eval(self.xk, self.c, float(x))

    def table(self, n=101):
        """[[x, f], ...] sampled over the knot range (viewer geometry)."""
        x = np.linspace(self.xk[0], self.xk[-1], n)
        return [[float(a), float(self(a)[0])] for a in x]


class Ramp(_Profile):
    """Flyweight radius r_w(u) over travel from full open u; call -> (r_w, kappa, kappa')."""

    @classmethod
    def straight(cls, r_w0, kappa, u_max):
        return cls([0.0, u_max], [[r_w0, kappa] + [0.0] * (NC - 2)], "straight")

    @classmethod
    def poly(cls, coeffs, u_max):
        if len(coeffs) > NC:
            raise ValueError(f"ramp polynomial degree must be <= {NC - 1}")
        a = np.zeros(NC)
        a[:len(coeffs)] = coeffs
        return cls([0.0, u_max], [a], "poly")

    @classmethod
    def pchip(cls, u, r_w):
        u, r_w = np.asarray(u, float), np.asarray(r_w, float)
        if np.any(np.diff(r_w) <= 0.0):
            raise ValueError("ramp table r_w(u) must be strictly increasing")
        c = np.zeros((len(u) - 1, NC))
        c[:, :4] = PchipInterpolator(u, r_w).c[::-1].T
        return cls(u, c, "pchip")

    @classmethod
    def from_spec(cls, spec, u_max):
        """spec: {type: straight, r_w0, kappa} | {type: poly, coeffs} | {type: table, table: [[u, r_w], ...]}."""
        t = spec["type"]
        if t == "straight":
            return cls.straight(spec["r_w0"], spec["kappa"], u_max)
        if t == "poly":
            return cls.poly(spec["coeffs"], u_max)
        if t == "table":
            a = np.asarray(spec["table"], float)
            return cls.pchip(a[:, 0], a[:, 1])
        raise ValueError(f"unknown ramp type {t!r}")


class Helix(_Profile):
    """Relative twist h(x) of the secondary movable half; call -> (h, h', h'')."""

    @classmethod
    def from_angles(cls, r_cam, x, alpha_deg, x_max, interp="steps", blend=2e-3):
        x = np.asarray(x, float)
        hp = 1.0 / (r_cam * np.tan(np.radians(np.asarray(alpha_deg, float))))
        if x[0] != 0.0 or np.any(np.diff(x) <= 0.0):
            raise ValueError("helix table must start at x = 0 and increase")
        # segments (x0, x1, h'(x0), h'(x1), shape) with shape 'lin' or 'smooth'
        seg = []
        if interp == "steps":
            if len(x) > 1 and (blend <= 0.0 or blend >= np.min(np.diff(x)) or x[1] - 0.5 * blend <= 0.0):
                raise ValueError("helix blend must be positive and shorter than the angle steps")
            a = 0.0
            for b, h0, h1 in zip(x[1:], hp[:-1], hp[1:]):
                seg += [(a, b - 0.5 * blend, h0, h0, "lin"), (b - 0.5 * blend, b + 0.5 * blend, h0, h1, "smooth")]
                a = b + 0.5 * blend
            if a < x_max or not seg:
                seg.append((a, max(x_max, a + blend), hp[-1], hp[-1], "lin"))
        elif interp == "linear":
            seg = [(x0, x1, h0, h1, "lin") for x0, x1, h0, h1 in zip(x[:-1], x[1:], hp[:-1], hp[1:])]
            if x[-1] < x_max or not seg:
                seg.append((x[-1], max(x_max, x[-1] + 1e-3), hp[-1], hp[-1], "lin"))
        else:
            raise ValueError(f"unknown helix interpolation {interp!r}")
        X = [seg[0][0]] + [q[1] for q in seg]
        c = np.zeros((len(seg), NC))
        H = 0.0
        for i, (x0, x1, h0, h1, shape) in enumerate(seg):
            d, dh = x1 - x0, h1 - h0
            if shape == "lin":   # h' = h0 + dh s/d
                c[i, :3] = H, h0, 0.5 * dh / d
                H += h0 * d + 0.5 * dh * d
            else:                # h' = h0 + dh (10 t^3 - 15 t^4 + 6 t^5), t = s/d
                c[i] = H, h0, 0.0, 0.0, 2.5 * dh / d ** 3, -3.0 * dh / d ** 4, dh / d ** 5
                H += h0 * d + 0.5 * dh * d
        return cls(X, c, interp)

    @classmethod
    def from_spec(cls, spec, x_max):
        """spec: {r_cam, angles: [[x, deg], ...], interp: steps|linear, blend}."""
        a = np.asarray(spec["angles"], float).reshape(-1, 2)
        return cls.from_angles(spec["r_cam"], a[:, 0], a[:, 1], x_max,
                               spec.get("interp", "steps"), spec.get("blend", 2e-3))
