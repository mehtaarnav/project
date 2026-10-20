"""Rolling bearing L10 life (ISO 281 basic rating life) and PTO side-load check.

Spec: docs/theory/12_fatigue_bearings.md section 3. Loads in N, C in N, speed in rpm.
C values come from cvtsim/io/materials.yaml (bearings section) and are PLACEHOLDER
until the part numbers are confirmed.
"""
from __future__ import annotations

from typing import Sequence

P_BALL = 3.0
P_ROLLER = 10.0 / 3.0


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def l10_revolutions(C: float, P: float, p: float = P_BALL) -> float:
    """Basic rating life L10 = (C/P)^p * 1e6 revolutions."""
    _require(C > 0 and P > 0, "C and P must be positive")
    return (C / P) ** p * 1e6


def equivalent_load(P_list: Sequence[float], revs: Sequence[float], p: float = P_BALL) -> float:
    """Equivalent load over a spectrum, weighted by revolutions at each load.

    P_eq = (sum P_i^p n_i / sum n_i)^(1/p), n_i revolutions at load P_i (ISO 281 averaging).
    """
    _require(len(P_list) == len(revs) and len(P_list) > 0, "P and revs must match and be non-empty")
    _require(all(n >= 0 for n in revs) and sum(revs) > 0, "revolutions must be >= 0 with positive sum")
    num = sum((P**p) * n for P, n in zip(P_list, revs))
    return (num / sum(revs)) ** (1.0 / p)


def mean_rpm(rpm_history: Sequence[float]) -> float:
    """Time-average speed of uniformly sampled rpm history."""
    _require(len(rpm_history) > 0, "empty rpm history")
    return sum(rpm_history) / len(rpm_history)


def l10_hours(C: float, P_eq: float, rpm_history: Sequence[float], p: float = P_BALL) -> float:
    """L10 life in hours: L10 / (60 * mean rpm). Needs a non-zero mean speed."""
    rpm = mean_rpm(rpm_history)
    _require(rpm > 0, "mean rpm must be positive")
    return l10_revolutions(C, P_eq, p) / (60.0 * rpm)


def pto_side_load_ratio(F_radial: Sequence[float], F_allow: float) -> float:
    """Max |F_radial| / F_allow. F_allow is the PTO allowable [N], PLACEHOLDER until Rehlko supplies it."""
    _require(F_allow > 0, "allowable side load must be positive")
    _require(len(F_radial) > 0, "empty load history")
    return max(abs(F) for F in F_radial) / F_allow
