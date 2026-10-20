"""Rainflow counting, shaft fatigue (Goodman, Marin, Basquin/Miner). Spec: docs/theory/12_fatigue_bearings.md.

Units: stresses in Pa inside the functions; Su, Sy in MPa only where a
Marin/Peterson coefficient needs them (documented per function). Diameters in
m for stress, mm for kb and the Peterson constant.
UNVERIFIED: Peterson a-constant, Marin Table 6-2 coefficients, Basquin 0.9 Su
anchor (see spec section 2.4 and 2.5).
"""
from __future__ import annotations

import math
from typing import Sequence

# Marin surface factor coefficients ka = a * Su^b, Su in MPa. UNVERIFIED (Shigley Table 6-2, memory).
MARIN_KA = {
    "ground": (1.58, -0.085),
    "machined": (4.51, -0.265),
    "hot_rolled": (57.7, -0.718),
    "forged": (272.0, -0.995),
}


def turning_points(x: Sequence[float]) -> list[float]:
    """Endpoints and local extrema of a sampled history, consecutive duplicates removed."""
    pts: list[float] = []
    for v in x:
        if not pts or v != pts[-1]:
            pts.append(float(v))
    if len(pts) < 3:
        return pts
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        if (pts[i] - pts[i - 1]) * (pts[i + 1] - pts[i]) < 0:
            out.append(pts[i])
    out.append(pts[-1])
    return out


def rainflow(x: Sequence[float]) -> list[tuple[float, float, float]]:
    """ASTM E1049-85 rainflow count (three-point method).

    Returns (range, mean, count) tuples; count is 1.0 for a full cycle and
    0.5 for a half cycle. The starting point is discarded as a half cycle when
    its range is first closed; the residue is returned as half cycles.
    """
    pts = turning_points(x)
    stack: list[float] = []
    start_alive = True
    out: list[tuple[float, float, float]] = []
    for p in pts:
        stack.append(p)
        while len(stack) >= 3:
            X = abs(stack[-1] - stack[-2])
            Y = abs(stack[-2] - stack[-3])
            if X < Y:
                break
            if len(stack) == 3 and start_alive:
                # Y contains the starting point: half cycle, discard start only.
                out.append((Y, 0.5 * (stack[0] + stack[1]), 0.5))
                stack.pop(0)
                start_alive = False
            else:
                # Full cycle: drop the two points forming Y, keep the last.
                out.append((Y, 0.5 * (stack[-2] + stack[-3]), 1.0))
                last = stack.pop()
                stack.pop()
                stack.pop()
                stack.append(last)
    for a, b in zip(stack[:-1], stack[1:]):
        out.append((abs(b - a), 0.5 * (a + b), 0.5))
    return out


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def peterson_a_mm(Su_MPa: float) -> float:
    """Peterson material constant a [mm], steel bending. UNVERIFIED form (Shigley, Peterson)."""
    _require(Su_MPa > 0, "Su must be positive")
    return 0.0254 * (2070.0 / Su_MPa) ** 1.8


def notch_sensitivity(r_mm: float, Su_MPa: float) -> float:
    """Peterson notch sensitivity q = 1 / (1 + a / r)."""
    _require(r_mm > 0, "notch radius must be positive")
    return 1.0 / (1.0 + peterson_a_mm(Su_MPa) / r_mm)


def fatigue_kf(Kt: float, q: float) -> float:
    """Fatigue stress-concentration factor Kf = 1 + q (Kt - 1)."""
    _require(Kt >= 1.0 and 0.0 <= q <= 1.0, "need Kt >= 1 and 0 <= q <= 1")
    return 1.0 + q * (Kt - 1.0)


def marin_ka(Su_MPa: float, finish: str) -> float:
    """Surface factor ka = a Su^b (Su in MPa)."""
    _require(finish in MARIN_KA, f"finish must be one of {sorted(MARIN_KA)}")
    a, b = MARIN_KA[finish]
    return a * Su_MPa**b


def marin_kb(d_mm: float) -> float:
    """Size factor kb = 1.24 d^-0.107 (d in mm), valid 2.79 <= d <= 51 mm."""
    _require(2.79 <= d_mm <= 51.0, "kb fit valid for 2.79 <= d <= 51 mm")
    return 1.24 * d_mm**-0.107


def endurance_limit(Su_MPa: float, finish: str, d_mm: float,
                    kc: float = 1.0, kd: float = 1.0, ke: float = 1.0, kf: float = 1.0) -> float:
    """Se = ka kb kc kd ke kf Se', Se' = 0.5 Su for Su <= 1400 MPa. Returns Pa."""
    _require(Su_MPa <= 1400.0, "Se' = 0.5 Su only for Su <= 1400 MPa")
    se_prime = 0.5 * Su_MPa
    se = marin_ka(Su_MPa, finish) * marin_kb(d_mm) * kc * kd * ke * kf * se_prime
    return se * 1e6


def bending_stress_amplitude(M_amp: float, d: float, Kf: float = 1.0) -> float:
    """Fully reversed bending stress amplitude Kf 32 M / (pi d^3) [Pa], mean 0."""
    _require(d > 0, "diameter must be positive")
    return Kf * 32.0 * abs(M_amp) / (math.pi * d**3)


def shaft_cycle_stresses(M_amp: float, T_hist: Sequence[float], d: float,
                         Kf: float, Kfs: float) -> list[tuple[float, float, float]]:
    """Per-cycle von Mises amplitude and mean for a rotating shaft.

    M_amp: bending moment amplitude [N m] from a fixed-direction load (fully
    reversed each revolution, mean 0). T_hist: torque history [N m], rainflowed.
    d: diameter [m]. Returns (sa', sm', count) with
    sa' = sqrt((Kf sa)^2 + 3 (Kfs ta)^2), sm' = sqrt(3) Kfs |tm|.
    """
    _require(d > 0, "diameter must be positive")
    sa_b = bending_stress_amplitude(M_amp, d, Kf)
    out = []
    for rng, mean, n in rainflow(T_hist):
        ta = 16.0 * (0.5 * rng) / (math.pi * d**3)
        tm = 16.0 * mean / (math.pi * d**3)
        sa = math.sqrt((sa_b) ** 2 + 3.0 * (Kfs * ta) ** 2)
        sm = math.sqrt(3.0) * Kfs * abs(tm)
        out.append((sa, sm, n))
    return out


def goodman_n(sa: float, sm: float, Se: float, Su: float) -> float:
    """Modified Goodman factor of safety n = 1 / (sa/Se + sm/Su)."""
    _require(Se > 0 and Su > 0, "Se and Su must be positive")
    denom = sa / Se + sm / Su
    _require(denom > 0, "amplitude and mean cannot both be zero")
    return 1.0 / denom


def yield_n(Sy: float, sa: float, sm: float) -> float:
    """Yield check n_y = Sy / (sa + sm) (Shigley, conservative peak-stress form)."""
    _require(sa + sm > 0, "peak stress must be positive")
    return Sy / (sa + sm)


def basquin_constants(Su: float, Se: float) -> tuple[float, float]:
    """S = A N^b through (1e3, 0.9 Su) and (1e6, Se). Returns (A, b). UNVERIFIED anchor."""
    _require(0 < Se < 0.9 * Su, "need 0 < Se < 0.9 Su")
    b = math.log10(Se / (0.9 * Su)) / 3.0
    A = 0.9 * Su * 1e3 ** (-b)
    return A, b


def cycles_to_failure(S: float, Su: float, Se: float) -> float:
    """N = (S / A)^(1/b) for a fully reversed amplitude S [same units as Su]."""
    A, b = basquin_constants(Su, Se)
    return (S / A) ** (1.0 / b)


def miner_damage(cycles: Sequence[tuple[float, float, float]], Su: float, Se: float) -> float:
    """Miner damage D = sum n_i / N_i over (sa, sm, count) cycles.

    Each cycle is converted to a fully reversed Goodman-equivalent amplitude
    sa / (1 - sm / Su). Equivalent amplitudes below Se contribute no damage
    (endurance cutoff). Units: sa, sm, Su, Se in the same unit (Pa).
    """
    D = 0.0
    for sa, sm, n in cycles:
        sa_eq = sa / (1.0 - sm / Su)
        if sa_eq < Se:
            continue
        D += n / cycles_to_failure(sa_eq, Su, Se)
    return D
