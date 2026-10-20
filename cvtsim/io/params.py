"""Parameter loading. YAML entries are {value:..., source:...}; `load()` returns a
flat-ish dict of plain values plus derived quantities, keeping sources separately."""
import pathlib
import numpy as np
import yaml

DEFAULT = pathlib.Path(__file__).with_name("dr26_default.yaml")


def _num(x):
    """YAML 1.1 reads exponents without a sign (e.g. 30.15e6) as strings."""
    if isinstance(x, str):
        try:
            return float(x)
        except ValueError:
            return x
    if isinstance(x, list):
        return [_num(v) for v in x]
    return x


def _strip(node):
    if isinstance(node, dict):
        if "value" in node:
            return _num(node["value"])
        return {k: _strip(v) for k, v in node.items() if k not in ("source", "unit", "note")}
    return _num(node)


def load_raw(path=DEFAULT):
    with open(path) as f:
        return yaml.safe_load(f)


def load(path=DEFAULT, overrides=None):
    """Plain-valued parameter dict with derived entries under p['derived']."""
    raw = load_raw(path)
    p = _strip(raw)
    if overrides:
        for key, val in overrides.items():
            d = p
            parts = key.split(".")
            for part in parts[:-1]:
                d = d[part]
            d[parts[-1]] = val
    from cvtsim.cvt.geometry import derive_geometry
    beta = np.radians(p["sheave"]["beta_deg"])
    g = derive_geometry(beta, p["primary"]["stroke"], p["secondary"]["stroke"],
                        p["ratio"]["low"], p["ratio"]["high"])
    s = p["secondary"]
    p_h = s["r_cam"] * np.tan(np.radians(s["cam_angle_deg"]))
    th0 = np.radians(s["torsion_preload_deg"])
    F_t0 = s["k_torsion"] * th0 / p_h
    p["derived"] = dict(geometry=g, beta=beta, p_h=p_h, theta_t0=th0,
                        F_torsion0=F_t0, F_l0=s["F_preload_equiv"] - F_t0,
                        m_w=p["primary"]["m_weight"] + p["primary"]["m_cover"])
    return p
