"""CVT belt geometry (docs/theory/05_cvt.md, section 2)."""
import numpy as np
from scipy.optimize import brentq


def _phi(r1, r2, C):
    return np.arcsin((r2 - r1) / C)


def belt_length(r1, r2, C):
    """Open-belt length (G1) for pitch radii r1 (primary), r2 (secondary), centre distance C."""
    phi = _phi(r1, r2, C)
    return 2.0 * C * np.cos(phi) + np.pi * (r1 + r2) + 2.0 * phi * (r2 - r1)


def _rs_closure(rp, L, C):
    """Secondary pitch radius that closes the belt length L at primary radius rp (G1)."""
    # belt_length is strictly increasing in r_s (dL/dr_s = pi + 2 phi > 0): unique root.
    return brentq(lambda rs: belt_length(rp, rs, C) - L, 1e-6, rp + 0.999 * C,
                  xtol=1e-15, rtol=1e-15, maxiter=500)


def derive_geometry(beta, L, C, rp_min, rs_max, s_p, i_high):
    """Operating geometry from measured hardware (05_cvt.md section 2).

    beta: flank half-angle [rad]; L: belt pitch length; C: centre distance;
    rp_min: primary minimum pitch radius (belt on the hub, sheave open);
    rs_max: secondary maximum pitch radius (sheave fully closed);
    s_p: primary mechanical stroke; i_high: overdrive ratio stop (r_s/r_p).
    """
    tb = np.tan(beta)
    rs_low = _rs_closure(rp_min, L, C)
    if rs_low > rs_max:
        raise ValueError("belt too long for the secondary at low ratio")
    rp_max = brentq(lambda rp: _rs_closure(rp, L, C) / rp - i_high,
                    rp_min, rp_min + s_p / (2.0 * tb), xtol=1e-15, rtol=1e-15, maxiter=500)
    rs_min = _rs_closure(rp_max, L, C)
    return dict(rp_min=rp_min, rp_max=rp_max, rs_low=rs_low, rs_min=rs_min, rs_max=rs_max,
                ratio_low=rs_low / rp_min, x_s0=2.0 * tb * (rs_max - rs_low),
                x_p_stop=2.0 * tb * (rp_max - rp_min), C=C, L=L, beta=beta, tan_beta=tb)


def r_star_primary(x_p, g):
    """Free primary pitch radius for primary closing travel x_p (04 section 2)."""
    return g["rp_min"] + x_p / (2.0 * g["tan_beta"])


def r_star_secondary(x_s, g):
    """Free secondary pitch radius for secondary opening travel x_s (04 section 2)."""
    return g["rs_max"] - x_s / (2.0 * g["tan_beta"])


def secondary_radius_from_xp(x_p, g):
    """Secondary pitch radius r_s that keeps the belt length at L (G1)."""
    return _rs_closure(r_star_primary(x_p, g), g["L"], g["C"])


def ratio_from_xp(x_p, g):
    """Speed ratio i = r_s / r_p at primary closing travel x_p."""
    rp = r_star_primary(x_p, g)
    return secondary_radius_from_xp(x_p, g) / rp
