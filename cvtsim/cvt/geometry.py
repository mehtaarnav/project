"""CVT belt geometry (docs/theory/05_cvt.md, section 2)."""
import numpy as np
from scipy.optimize import brentq


def _phi(r1, r2, C):
    return np.arcsin((r2 - r1) / C)


def belt_length(r1, r2, C):
    """Open-belt length (G1) for pitch radii r1 (primary), r2 (secondary), centre distance C."""
    phi = _phi(r1, r2, C)
    return 2.0 * C * np.cos(phi) + np.pi * (r1 + r2) + 2.0 * phi * (r2 - r1)


def derive_geometry(beta, s_p, s_s, i_low, i_high):
    """Derive belt-radius limits, centre distance C and belt length L.

    beta: sheave half-angle in radians (flank measured from the radial plane).
    s_p, s_s: primary closing / secondary opening strokes [m].
    i_low, i_high: ratio limits r_s/r_p at the low and high ends.
    """
    tb = np.tan(beta)
    dp = s_p / (2.0 * tb)
    ds = s_s / (2.0 * tb)
    rp_min = (ds + i_high * dp) / (i_low - i_high)
    rp_max = rp_min + dp
    rs_max = i_low * rp_min
    rs_min = i_high * rp_max

    def f(C):
        return belt_length(rp_min, rs_max, C) - belt_length(rp_max, rs_min, C)

    lo = max(abs(rs_max - rp_min), abs(rs_min - rp_max)) * (1.0 + 1e-9)
    C = brentq(f, lo, 1.0, xtol=1e-15, rtol=1e-15, maxiter=500)
    L = belt_length(rp_min, rs_max, C)
    return dict(rp_min=rp_min, rp_max=rp_max, rs_min=rs_min, rs_max=rs_max,
                C=C, L=L, beta=beta, tan_beta=tb, delta_p=dp, delta_s=ds)


def r_star_primary(x_p, g):
    """Free primary pitch radius for primary closing travel x_p (04 section 2)."""
    return g["rp_min"] + x_p / (2.0 * g["tan_beta"])


def r_star_secondary(x_s, g):
    """Free secondary pitch radius for secondary opening travel x_s (04 section 2)."""
    return g["rs_max"] - x_s / (2.0 * g["tan_beta"])


def secondary_radius_from_xp(x_p, g):
    """Secondary pitch radius r_s that keeps the belt length at L (G1)."""
    rp = r_star_primary(x_p, g)
    C, L = g["C"], g["L"]

    def f(rs):
        return belt_length(rp, rs, C) - L

    # belt_length is strictly increasing in r_s (dL/dr_s = pi + 2 phi > 0), so the root is unique.
    lo = max(1e-9, rp - 0.99 * C)
    hi = rp + 0.99 * C
    return brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=500)


def ratio_from_xp(x_p, g):
    """Speed ratio i = r_s / r_p at primary closing travel x_p."""
    rp = r_star_primary(x_p, g)
    return secondary_radius_from_xp(x_p, g) / rp
