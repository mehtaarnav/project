"""Belt bending fatigue and heat-fade diagnostic (docs/theory/13 sec 3-4). Post-processing only.

All constants are PLACEHOLDER or ASSUMED (cvtsim/io/belt_life.yaml). A life is a relative
index until N_ref, eps_ref and k are calibrated; do not read it as a prediction.
`p` comes from inputs(params): the YAML inputs plus C, h_n, mu0 from the simulation parameters.
"""
import pathlib

import numpy as np
import yaml

from cvtsim.post.efficiency import _col

YAML = pathlib.Path(__file__).resolve().parents[1] / "io" / "belt_life.yaml"


def inputs(params):
    """Plain input dict: belt_life.yaml values plus geometry and belt entries of `params`
    (cvtsim.io.params.load()). h_n = h_n_fraction * flank_height (13 sec 3.2)."""
    with open(YAML) as f:
        raw = yaml.safe_load(f)
    p = {}
    for sec in raw.values():
        for key, node in sec.items():
            v = node["value"]
            p[key] = float(v) if isinstance(v, str) else v   # YAML 1.1: 1.0e6 reads as str
    p["C"] = float(params["derived"]["geometry"]["C"])
    p["mu0"] = float(params["belt"]["mu"])
    p["h_n"] = p["h_n_fraction"] * float(params["belt"]["flank_height"])
    return p


def wrap_angle_and_length(r_self, r_other, C):
    """Wrap angle of an open belt on a sheave of effective radius r_self, the other sheave at
    r_other, centre distance C: theta = pi - 2 asin((r_other - r_self)/C) (argument clipped to
    [-1, 1], so theta stays in [0, 2 pi]). Returns (theta, L = r_self * theta)."""
    theta = np.pi - 2.0 * np.arcsin(np.clip((r_other - r_self) / C, -1.0, 1.0))
    return theta, r_self * theta


def bending_strain(h_n, r):
    """Outer-fibre bending strain over a sheave of radius r: eps = h_n / r."""
    return h_n / r


def temperature_factor(T, T_ref, dT_half):
    """Damage multiplier 2**((T - T_ref)/dT_half): life halves per dT_half above T_ref (ASSUMED)."""
    return 2.0 ** ((T - T_ref) / dT_half)


def _sheave_damage(res, CH, p):
    """Per-sheave damage rates [1/s] and the belt temperature factor, arrays over the run."""
    C = p["C"]
    T = _col(res, CH, "T_belt")
    fT = temperature_factor(T, p["T_ref"], p["dT_half"])
    out = {}
    for k, other, w in (("p", "s", "w_e"), ("s", "p", "w_s")):
        r = _col(res, CH, f"r_eff_{k}")
        _, L = wrap_angle_and_length(np.where(r > 0.0, r, 1.0), _col(res, CH, f"r_eff_{other}"), C)
        contact = (r > 0.0) & (L > 0.0)                     # no belt on this sheave: no damage
        n_dot = np.where(contact, r * np.abs(_col(res, CH, w)) / np.where(contact, L, 1.0), 0.0)
        eps = bending_strain(p["h_n"], np.where(r > 0.0, r, 1.0))
        per_cycle = (eps / p["eps_ref"]) ** p["k"] / p["N_ref"]   # 1 / N_f (Miner)
        out[k] = fT * n_dot * per_cycle
    return out


def damage_rate(res, CH, p):
    """Total damage rate Ddot(t) = sum over sheaves [1/s] (13 sec 3.3)."""
    d = _sheave_damage(res, CH, p)
    return d["p"] + d["s"]


def _trapz(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def belt_damage(res, CH, p):
    """Accumulated Miner damage D over the run and projected life (13 sec 3.4).

    life_s = T_run / (D * duty): repeated launches to D = 1, if `duty` of operating time is
    spent in launch-like loading. Returns dict with D, D_p, D_s, T_run, life_s, life_h.
    """
    t = _col(res, CH, "t")
    d = _sheave_damage(res, CH, p)
    D = _trapz(d["p"] + d["s"], t)
    D_p = _trapz(d["p"], t)
    D_s = _trapz(d["s"], t)
    T_run = float(t[-1] - t[0])
    life_s = T_run / (D * p["duty"]) if D > 0.0 else float("inf")
    return dict(D=D, D_p=D_p, D_s=D_s, T_run=T_run, life_s=life_s, life_h=life_s / 3600.0)


def friction_fade(res, CH, p):
    """Diagnostic only, not fed back (13 sec 4): mu(T_b) = mu0 (1 - c (T_b - T_ref)) for
    c = c_lo and c = c_hi. Returns dict with arrays mu_lo, mu_hi (same length as res)."""
    T = _col(res, CH, "T_belt")
    dT = T - p["T_ref"]
    return dict(mu_lo=p["mu0"] * (1.0 - p["c_lo"] * dT),
                mu_hi=p["mu0"] * (1.0 - p["c_hi"] * dT))
