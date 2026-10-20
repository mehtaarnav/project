"""Peterson-type stress-concentration factors (curve-fit approximations).

The curve-fit coefficients are UNVERIFIED against Pilkey & Pilkey, Peterson's
Stress Concentration Factors, 3rd ed. (docs/theory/07_shafts.md, section 3).
Do not use for design until they are checked against the charts.
"""
from __future__ import annotations

import math


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def kt_shoulder_bending(D: float, d: float, r: float) -> float:
    """Shoulder fillet, bending. Valid for 0.1 <= t/r <= 2, t = (D - d)/2."""
    _require(D > d > 0 and r > 0, "need D > d > 0 and r > 0")
    t = 0.5 * (D - d)
    tr = t / r
    _require(0.1 - 1e-12 <= tr <= 2.0 + 1e-12, f"bending fit valid for 0.1 <= t/r <= 2, got {tr:g}")
    s = math.sqrt(tr)
    x = 2.0 * t / D
    c1 = 0.947 + 1.206 * s - 0.131 * tr
    c2 = 0.022 - 3.405 * s + 0.915 * tr
    c3 = 0.869 + 1.777 * s - 0.555 * tr
    c4 = -0.810 + 0.422 * s - 0.260 * tr
    return c1 + c2 * x + c3 * x**2 + c4 * x**3


def kt_shoulder_torsion(D: float, d: float, r: float) -> float:
    """Shoulder fillet, torsion. Valid for 0.25 <= t/r <= 4, t = (D - d)/2."""
    _require(D > d > 0 and r > 0, "need D > d > 0 and r > 0")
    t = 0.5 * (D - d)
    tr = t / r
    _require(0.25 - 1e-12 <= tr <= 4.0 + 1e-12, f"torsion fit valid for 0.25 <= t/r <= 4, got {tr:g}")
    s = math.sqrt(tr)
    x = 2.0 * t / D
    c1 = 0.905 + 0.783 * s - 0.075 * tr
    c2 = -0.437 - 1.969 * s + 0.553 * tr
    c3 = 1.557 + 1.073 * s - 0.578 * tr
    c4 = -1.061 + 0.171 * s + 0.086 * tr
    return c1 + c2 * x + c3 * x**2 + c4 * x**3


def kt_hole_plate_tension(d: float, H: float) -> float:
    """Net-section Kt of a finite-width plate with a central hole, tension."""
    _require(H > 0 and 0 <= d < H, "need 0 <= d < H")
    q = 1.0 - d / H
    return 2.0 + 0.284 * q - 0.600 * q**2 + 1.32 * q**3
