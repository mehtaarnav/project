"""Parameter loading. YAML entries are {value:..., source:...}; `load()` returns a
flat-ish dict of plain values plus derived quantities, keeping sources separately.

A tune (docs/theory/09 section 5) selects named parts from `tuning/components.yaml` and
adds scalar adjustments; `load(tune=...)` merges it over the base file before overrides."""
import pathlib
import numpy as np
import yaml

DEFAULT = pathlib.Path(__file__).with_name("dr26_default.yaml")
TUNING = pathlib.Path(__file__).with_name("tuning")
COMPONENTS = TUNING / "components.yaml"

# tune slot -> (components.yaml section, parameter section it fills)
SLOTS = {"primary_spring": ("primary_springs", "primary"),
         "flyweights": ("flyweights", "primary"),
         "ramp": ("ramps", "primary"),
         "secondary_spring": ("secondary_springs", "secondary"),
         "torsion_spring": ("torsion_springs", "secondary"),
         "helix": ("helixes", "secondary"),
         "surface": ("vehicle_surfaces", "vehicle")}
_SECONDARY_PRELOAD = ("F_preload_equiv", "F_preload_linear", "equiv_reference")


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
        return {k: _strip(v) for k, v in node.items() if k not in ("source", "unit", "note", "description")}
    return _num(node)


def load_raw(path=DEFAULT):
    with open(path) as f:
        return yaml.safe_load(f)


def tune_path(tune):
    """Tune name (file in tuning/tunes/) or path to a tune YAML."""
    p = pathlib.Path(tune)
    return p if p.suffix in (".yaml", ".yml") and p.exists() else TUNING / "tunes" / f"{tune}.yaml"


SCENARIOS = pathlib.Path(__file__).with_name("scenarios")


def scenario_names():
    return sorted(q.stem for q in SCENARIOS.glob("*.yaml"))


def load_scenario(name):
    """Scenario file (cvtsim/io/scenarios/<name>.yaml or a path) as a plain dict."""
    q = pathlib.Path(name)
    q = q if q.suffix in (".yaml", ".yml") and q.exists() else SCENARIOS / f"{name}.yaml"
    d = load_raw(q)
    d["name"] = q.stem
    return d


def tune_names():
    return sorted(q.stem for q in (TUNING / "tunes").glob("*.yaml"))


def component(section, name):
    """Plain values of one part in components.yaml, e.g. component("vehicle_surfaces", "dirt")."""
    return _strip(load_raw(COMPONENTS)[section][name])


def load_tune_raw(tune):
    """The tune file plus the raw (sourced) entries of the parts it selects, for run configs."""
    t = load_raw(tune_path(tune))
    comps = load_raw(COMPONENTS)
    parts = {slot: {name: comps[SLOTS[slot][0]][name]} for slot, name in t["components"].items()}
    return dict(tune=t, parts=parts)


def _apply_tune(p, tune):
    tr = load_tune_raw(tune)
    t = tr["tune"]
    for slot, part in tr["parts"].items():
        if slot not in SLOTS:
            raise KeyError(f"unknown tune slot {slot!r}")
        (name, entry), = part.items()
        vals = _strip(entry)
        sec = p[SLOTS[slot][1]]
        if slot == "ramp":
            sec["ramp"] = vals
        elif slot == "helix":
            sec["helix"] = vals
            sec["r_cam"] = vals["r_cam"]
            sec["cam_angle_deg"] = vals["angles"][0][1]
        else:
            if slot == "secondary_spring":
                for key in _SECONDARY_PRELOAD:
                    sec.pop(key, None)
            sec.update(vals)
    adj = _strip(t.get("adjust", {}))
    if "primary_shim" in adj:
        p["primary"]["shim"] = adj.pop("primary_shim")
    if "torsion_hole_deg" in adj:
        p["secondary"]["torsion_preload_deg"] = adj.pop("torsion_hole_deg")
    if adj:
        raise KeyError(f"unknown tune adjustments {sorted(adj)}")
    for sec in ("engine", "vehicle", "belt", "sim", "ratio"):
        p[sec].update(_strip(t.get(sec, {})))
    p["tune"] = dict(name=pathlib.Path(tune).stem, components=dict(t["components"]))


def check_high_ratio(g, gm, stroke, i_high):
    """Reject an overdrive stop the hardware cannot reach (primary radius, secondary radius, travel)."""
    bad = []
    if g["rp_max"] > gm["rp_mech_max"] + 1e-9:
        bad.append(f"primary radius {g['rp_max'] * 1e3:.2f} mm > limit {gm['rp_mech_max'] * 1e3:.2f} mm")
    if g["rs_min"] < gm["rs_mech_min"] - 1e-9:
        bad.append(f"secondary radius {g['rs_min'] * 1e3:.2f} mm < limit {gm['rs_mech_min'] * 1e3:.2f} mm")
    if g["x_p_stop"] > stroke + 1e-9:
        bad.append(f"primary travel {g['x_p_stop'] * 1e3:.2f} mm > stroke {stroke * 1e3:.2f} mm")
    if bad:
        raise ValueError(f"ratio.high = {i_high} is not reachable: " + "; ".join(bad))


def load(path=DEFAULT, overrides=None, tune=None):
    """Plain-valued parameter dict with derived entries under p['derived'].

    tune: name in tuning/tunes/ (or a path); its parts are merged before `overrides`.
    """
    raw = load_raw(path)
    p = _strip(raw)
    if tune is not None:
        _apply_tune(p, tune)
    if overrides:
        for key, val in overrides.items():
            d = p
            parts = key.split(".")
            for part in parts[:-1]:
                d = d[part]
            d[parts[-1]] = val
    from cvtsim.cvt.geometry import derive_geometry
    from cvtsim.cvt.profiles import Ramp, Helix
    beta = np.radians(p["sheave"]["beta_deg"])
    gm = p["geometry"]
    g = derive_geometry(beta, gm["belt_pitch_length"], gm["center_distance"],
                        gm["rp_mech_min"], gm["rs_mech_max"], p["primary"]["stroke"],
                        p["ratio"]["high"])
    check_high_ratio(g, gm, p["primary"]["stroke"], p["ratio"]["high"])
    pr, s = p["primary"], p["secondary"]
    # primary: shim offset adds k_p * shim to the preload (09 section 1)
    shim = pr.get("shim", 0.0)
    if shim:
        pr["F_preload_spring"] = pr["F_preload"]
        pr["F_preload"] = pr["F_preload"] + pr["k_spring"] * shim
    u_max = pr["stroke"]  # ramp coordinate u = x_p from belt contact; u < 0 extrapolates linearly
    ramp = Ramp.from_spec(pr["ramp"], u_max) if "ramp" in pr else Ramp.straight(pr["r_w0"], pr["kappa"], u_max)
    if "ramp" in pr:  # static 05 functions read the ramp at belt contact (u = 0)
        pr["r_w0"], pr["kappa"] = float(ramp(0.0)[0]), float(ramp(0.0)[1])
    # secondary
    if "torsion_holes_deg" in s and s["torsion_preload_deg"] not in s["torsion_holes_deg"]:
        raise ValueError(f"torsion pretwist {s['torsion_preload_deg']} deg is not a hole position {s['torsion_holes_deg']}")
    if "helix" in s:
        helix = Helix.from_spec(s["helix"], s["stroke"])
    else:
        helix = Helix.from_angles(s["r_cam"], [0.0], [s["cam_angle_deg"]], s["stroke"])
    p_h = s["r_cam"] * np.tan(np.radians(s["cam_angle_deg"]))
    th0 = np.radians(s["torsion_preload_deg"])
    F_t0 = s["k_torsion"] * th0 / p_h
    if "F_preload_linear" in s:
        F_l0 = s["F_preload_linear"]
    else:
        # brief convention: the equivalent preload includes the torsion share at a reference set-up
        ref = s.get("equiv_reference", s)
        F_l0 = s["F_preload_equiv"] - ref["k_torsion"] * np.radians(ref["torsion_preload_deg"]) / (
            ref["r_cam"] * np.tan(np.radians(ref["cam_angle_deg"])))
    b = p["belt"]
    m_belt = b["rho"] * b["area"] * gm["belt_pitch_length"]
    hA = b["h_conv"] * gm["belt_pitch_length"] * (b["top_width"] + 2.0 * b["flank_height"])
    thermal = dict(m_belt=m_belt, C_b=m_belt * b["c_p"], hA=hA, tau=m_belt * b["c_p"] / hA,
                   T_amb=b["T_amb_C"] + 273.15)
    # engine rotor split and PTO torsion (11 section 2); J_e is the one-rotor inertia of 05/09
    e = p["engine"]
    if e["pto_model"] not in ("torsional", "rigid"):
        raise ValueError(f"engine.pto_model must be 'torsional' or 'rigid', not {e['pto_model']!r}")
    e["J_e"] = e["J_f"] + e["J_c"]
    mat = p["shafts"]["material"]
    G_sh = mat["E"] / (2.0 * (1.0 + mat["nu"]))
    d_eff = e["pto_d"] - 0.5 * e["key_width"]
    L_eff = e["pto_length"] + e["crank_web_length"]
    k_pto = G_sh * np.pi * d_eff ** 4 / (32.0 * L_eff)
    J_red = e["J_f"] * e["J_c"] / (e["J_f"] + e["J_c"])
    pto = dict(G=G_sh, d_eff=d_eff, L_eff=L_eff, k=k_pto,
               c=2.0 * e["pto_damping_ratio"] * np.sqrt(k_pto * J_red),
               f_n=np.sqrt(k_pto / J_red) / (2.0 * np.pi))
    p["derived"] = dict(geometry=g, beta=beta, p_h=p_h, theta_t0=th0,
                        F_torsion0=F_t0, F_l0=F_l0,
                        m_w=pr["m_weight"] + pr["m_cover"], ramp=ramp, helix=helix,
                        belt_thermal=thermal, pto=pto)
    return p
