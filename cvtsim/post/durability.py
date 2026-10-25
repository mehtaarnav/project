"""Launch durability: shaft fatigue, bearing life and PTO side load from one launch.

Spec: docs/theory/12_fatigue_bearings.md section 7. Inputs are the shaft histories of
cvtsim.fe.shaft_loads.shaft_histories and the constants in cvtsim/io/materials.yaml.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import yaml

from cvtsim.fe import bearings as brg
from cvtsim.fe.fatigue import basquin_constants, endurance_limit, fatigue_kf, notch_sensitivity, rainflow

MATERIALS = Path(__file__).resolve().parents[1] / "io" / "materials.yaml"
SQ3 = math.sqrt(3.0)


def load_materials(path=MATERIALS):
    with open(path) as f:
        return yaml.safe_load(f)


def revolutions(t, w):
    """Revolutions per sample interval, dN_i = (|w_i| + |w_i+1|)/2 dt / 2 pi (7.1)."""
    w = np.abs(np.asarray(w, float))
    return 0.5 * (w[1:] + w[:-1]) * np.diff(np.asarray(t, float)) / (2.0 * np.pi)


def _mid(x):
    x = np.asarray(x, float)
    return 0.5 * (x[1:] + x[:-1])


def _damage(sa, sm, n, Su, Se):
    """Miner sum of n cycles at (sa, sm) arrays; Goodman-equivalent amplitude, endurance cutoff (2.5)."""
    A, b = basquin_constants(Su, Se)
    seq = sa / np.maximum(1.0 - sm / Su, 1e-12)
    seq = np.where(sm >= Su, np.inf, seq)
    hit = seq >= Se
    return float(np.sum(n[hit] / (seq[hit] / A) ** (1.0 / b)))


def _window_stats(T, win):
    """Moving half-range and mean of T over a centred window of `win` samples (7.3)."""
    from scipy.ndimage import maximum_filter1d, minimum_filter1d, uniform_filter1d
    T = np.asarray(T, float)
    win = max(int(win), 1)
    return (0.5 * (maximum_filter1d(T, win, mode="nearest") - minimum_filter1d(T, win, mode="nearest")),
            uniform_filter1d(T, win, mode="nearest"))


def shaft_fatigue(t, w, M, T, crit, Kt, Su, Sy, finish, win=1):
    """Fatigue and yield of one critical section (7.1-7.3). Stresses in Pa, Su/Sy in Pa.

    crit: dict(d, r) [m]; Kt: (Kt_bending, Kt_torsion); win: samples in one firing period, the
    window for the local torque amplitude and mean. Returns dict of scalars.
    """
    d = crit["d"]
    q = notch_sensitivity(crit["r"] * 1e3, Su / 1e6)
    Kf, Kfs = fatigue_kf(Kt[0], q), Kt[1]
    Se = endurance_limit(Su / 1e6, finish, d * 1e3)
    zb, zt = 32.0 / (math.pi * d ** 3), 16.0 / (math.pi * d ** 3)
    M, T = np.asarray(M, float), np.asarray(T, float)
    dN = revolutions(t, w)
    # bending family: one fully reversed cycle per revolution, mean from the torque
    sa_b = Kf * zb * np.abs(_mid(M))
    sm_b = SQ3 * Kfs * zt * np.abs(_mid(T))
    D_b = _damage(sa_b, sm_b, dN, Su, Se)
    seq_b = sa_b / np.maximum(1.0 - sm_b / Su, 1e-12)
    rev_above = float(dN[seq_b >= Se].sum())
    # torsion family: rainflow of T(t)
    cyc = rainflow(T) if len(T) > 1 else []
    if cyc:
        rng, mean, cnt = (np.array(c, float) for c in zip(*cyc))
        D_t = _damage(SQ3 * Kfs * zt * 0.5 * rng, SQ3 * Kfs * zt * np.abs(mean), cnt, Su, Se)
    else:
        D_t = 0.0
    # Goodman on simultaneous stresses: bending at each sample, torque amplitude and mean over
    # one firing period around it (7.3). Its minimum is an infinite-life check.
    ta, tm = _window_stats(T, win)
    sa = np.sqrt((Kf * zb * M) ** 2 + 3.0 * (Kfs * zt * ta) ** 2)
    sm = SQ3 * Kfs * zt * np.abs(tm)
    n_t = 1.0 / np.maximum(sa / Se + sm / Su, 1e-300)
    i = int(np.argmin(n_t))
    # yield: peak nominal von Mises (no notch factors, ductile material)
    vm = np.sqrt((zb * M) ** 2 + 3.0 * (zt * T) ** 2)
    D = D_b + D_t
    return dict(Kf=Kf, Kfs=Kfs, q=q, Se_MPa=Se / 1e6, N_rev=float(dN.sum()),
                D_bend=D_b, D_tors=D_t, D=D,
                launches_to_failure=(1.0 / D if D > 0 else None),  # None: no cycle above Se, life not bounded by this model
                rev_above_Se=rev_above,
                n_goodman=float(n_t[i]), t_goodman_min=float(np.asarray(t)[i]),
                sa_MPa=float(sa[i]) / 1e6, sm_MPa=float(sm[i]) / 1e6,
                vm_nominal_peak_MPa=float(vm.max()) / 1e6,
                n_yield=(Sy / vm.max() if vm.max() > 0 else math.inf))


def bearing_life(t, w, R, C=None, p=brg.P_BALL, target_launches=None, target_hours=None):
    """Launch-severity bearing results (7.4). R: radial reaction history [N]."""
    dN = revolutions(t, w)
    N = float(dN.sum())
    P = brg.equivalent_load(list(_mid(R)), list(dN), p) if N > 0 else 0.0
    rpm = N / (t[-1] - t[0]) * 60.0
    out = dict(P_eq_N=P, R_peak_N=float(np.max(R)), N_rev=N, rpm_mean=rpm)
    if target_launches:
        out["C_req_N"] = P * (target_launches * N / 1e6) ** (1.0 / p)
    if target_hours:
        out["C_req_hours_N"] = P * (60.0 * rpm * target_hours / 1e6) ** (1.0 / p)
    if C:
        out["launches_to_L10"] = brg.l10_revolutions(C, P, p) / N
        out["L10_h"] = brg.l10_revolutions(C, P, p) / (60.0 * rpm)
    return out


def side_load(Fx, Fy, a, F_allow=None):
    """PTO side load at the belt line and moment at the mounting face (7.5)."""
    F = np.hypot(Fx, Fy)
    out = dict(F_hub_peak_N=float(F.max()), M_face_peak_Nm=float(F.max() * a))
    if F_allow:
        out["side_load_ratio"] = brg.pto_side_load_ratio(list(F), F_allow)
    return out


def sampling_ok(t, rpm_max, f_pto=None):
    """True if the output resolves the firing torque (and the PTO mode when given) (7.6)."""
    dt = float(np.median(np.diff(t)))
    ok = dt <= 120.0 / (20.0 * max(rpm_max, 1.0))
    if f_pto:
        ok = ok and dt <= 1.0 / (10.0 * f_pto)
    return ok


def _rating(mats, key):
    if key is None:
        return None, brg.P_BALL
    b = mats["bearings"][key]
    return b["C_N"]["value"], (brg.P_BALL if b["type"] == "ball" else brg.P_ROLLER)


def _material(mats, key):
    m = mats["materials"][key]
    Su, Sy = m["Su_MPa"]["value"], m["Sy_MPa"]["value"]
    if Su is None or Sy is None:
        raise ValueError(f"material {key!r} has no Su/Sy in materials.yaml")
    return Su * 1e6, Sy * 1e6


def launch_durability(res, CH, sh, p, mats=None):
    """All durability results of one launch, as a nested dict (7)."""
    mats = mats or load_materials()
    c = lambda k: res[:, CH.index(k)]
    t = c("t")
    du = p["durability"]
    out = dict(sampled_ok=sampling_ok(t, (c("w_e") * 30 / np.pi).max(),
                                      p["derived"].get("pto", {}).get("f_n")))
    for name, tag, w in (("primary", "prim", c("w_c")), ("secondary", "sec", c("w_s"))):
        cfg, sd = du[name], p["shafts"][name]
        Su, Sy = _material(mats, cfg["material"])
        w_run = np.median(np.abs(c("w_e"))[c("w_e") > 1.0]) if (c("w_e") > 1.0).any() else 1.0
        win = 4.0 * np.pi / w_run / float(np.median(np.diff(t)))  # one firing period, 4-stroke
        fat = shaft_fatigue(t, w, sh[f"M_{tag}"], sh[f"T_{tag}"], sd["critical"],
                            sh["kt"][name], Su, Sy, cfg["finish"], win)
        brs = []
        for j, key in enumerate(cfg["bearings"]):
            C, pe = _rating(mats, key)
            brs.append(dict(name=key, **bearing_life(t, w, sh[f"F_bearing_{tag}{j+1}"], C, pe,
                                                    du["target_launches"])))
        out[name] = dict(fatigue=fat, bearings=brs)
    F_allow = mats["pto"]["allowable_side_load_N"]["value"]
    out["pto"] = side_load(c("Fx_p"), c("Fy_p"), du["belt_line_offset"], F_allow)
    return out


def summary_metrics(d):
    """Flat scalar metrics for tune runs and the uncertainty bands."""
    m = dict(pto_side_peak_N=d["pto"]["F_hub_peak_N"], pto_face_moment_Nm=d["pto"]["M_face_peak_Nm"],
             durability_sampled_ok=float(d["sampled_ok"]))
    for name, tag in (("primary", "prim"), ("secondary", "sec")):
        f = d[name]["fatigue"]
        m[f"n_goodman_{tag}"] = f["n_goodman"]
        m[f"n_yield_{tag}"] = f["n_yield"]
        m[f"damage_{tag}"] = f["D"]
        m[f"rev_above_Se_{tag}"] = f["rev_above_Se"]
        for j, b in enumerate(d[name]["bearings"]):
            m[f"P_eq_{tag}{j+1}_N"] = b["P_eq_N"]
            m[f"C_req_{tag}{j+1}_N"] = b.get("C_req_N", math.nan)
            if "launches_to_L10" in b:
                m[f"L10_launches_{tag}{j+1}"] = b["launches_to_L10"]
    return m
