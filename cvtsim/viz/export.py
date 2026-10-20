"""Viewer data export (docs/viewer_contract.md, "Python -> viewer")."""
import json
import pathlib

import numpy as np

SIG = 5
MAX_SERIES = 2000


def _sig(a, n=SIG):
    """Round to n significant figures (zeros and non-finite values pass through)."""
    a = np.asarray(a, dtype=float)
    out = a.copy()
    nz = np.isfinite(a) & (a != 0)
    mag = np.floor(np.log10(np.abs(a[nz])))
    out[nz] = np.round(a[nz] * 10.0 ** (n - 1 - mag)) / 10.0 ** (n - 1 - mag)
    return out


def _lst(a):
    """Nested list rounded to 5 significant figures; NaN/inf become null."""
    r = _sig(a)
    if r.ndim == 0:
        return float(r) if np.isfinite(r) else None
    if r.ndim == 1 and np.isfinite(r).all():
        return r.tolist()
    return [_lst(x) for x in r]


def _decimate(n, m=MAX_SERIES):
    return np.arange(n) if n <= m else np.unique(np.round(np.linspace(0, n - 1, m)).astype(int))


def _table(d):
    """Dict of equal-length arrays, decimated to MAX_SERIES points, as lists."""
    idx = _decimate(len(np.asarray(next(iter(d.values())))))
    return {k: _lst(np.asarray(v, dtype=float)[idx]) for k, v in d.items()}


def _round(x):
    if isinstance(x, dict):
        return {k: _round(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_round(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return _lst(float(x))
    return x


def _geometry(p):
    g = p["derived"]["geometry"]
    pr, sc, bt, gm = p["primary"], p["secondary"], p["belt"], p["geometry"]
    gap = pr["engagement_gap"]
    if "ramp" in pr:
        ramp = pr["ramp"]
    else:
        u = np.linspace(0.0, gap + g["x_p_stop"], 21)
        ramp = np.column_stack([u, pr["r_w0"] + pr["kappa"] * u]).tolist()
    if "helix" in sc:
        helix = sc["helix"]
    else:
        x = np.linspace(0.0, sc["stroke"], 21)
        helix = np.column_stack([x, x / p["derived"]["p_h"]]).tolist()
    return _round({
        "C": g["C"], "beta_deg": p["sheave"]["beta_deg"], "rp_min": g["rp_min"],
        "rp_mech_max": gm["rp_mech_max"], "rs_mech_min": gm["rs_mech_min"],
        "rs_mech_max": gm["rs_mech_max"], "belt_width": bt["top_width"],
        "belt_thickness": bt["flank_height"], "gap": gap, "x_p_stop": g["x_p_stop"],
        "s_s": sc["stroke"], "n_weights": pr["n_weights"], "r_w0": pr["r_w0"],
        "ramp": ramp, "helix": helix, "r_cam": sc["r_cam"],
        "spring_p": {"k": pr["k_spring"], "F0": pr["F_preload"]},
        "spring_s": {"k": sc["k_linear"], "F0": p["derived"]["F_l0"], "k_t": sc["k_torsion"],
                     "theta_t0_deg": sc["torsion_preload_deg"]},
    })


def _meta(run_dir, p):
    commit = ""
    mp = run_dir / "meta.json"
    if mp.exists():
        commit = json.loads(mp.read_text()).get("commit", "")
    return _round({"run": run_dir.name, "commit": commit, "dt": p["sim"]["dt"],
                   "tune": p.get("tune", "default")})


def write_viewer(run_dir, res, CH, snaps, p, sh=None):
    """Write run_dir/viewer.json and return its path.

    res: (n_out, len(CH)) channel array; snaps: simulate() frame dict; p: params.load() dict;
    sh: optional dict of equal-length arrays including "t" (bearing and stress channels).
    """
    run_dir = pathlib.Path(run_dir)
    res = np.asarray(res, dtype=float)
    out = {
        "meta": _meta(run_dir, p),
        "geometry": _geometry(p),
        "frames": {k: _lst(snaps[k]) for k in
                   ("t", "belt", "tension", "x_p", "x_s", "theta_e", "theta_s", "r_w")},
        "series": _table({c: res[:, i] for i, c in enumerate(CH)}),
        "shafts": _table(sh) if sh else {},
    }
    path = run_dir / "viewer.json"
    path.write_text(json.dumps(out, separators=(",", ":")))
    return path
