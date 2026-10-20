"""Reader for the team's CVT tuning workbook (CVT_Inputs template used by CVT_model_v4).

Sheets: Configurations, Primary Spring, Secondary Spring, Flyweights, Secondary Cam,
Constants, Dyno Data. Every value is converted to SI using the units row of its sheet.
`configuration(path, name)` returns the parameter overrides for one configuration in the
dotted-key form accepted by `cvtsim.io.params.load(overrides=...)`.
"""
import numpy as np
import openpyxl

IN, LBF, G = 0.0254, 4.4482216152605, 1e-3
TO_SI = {"in": IN, "m": 1.0, "mm": 1e-3, "g": G, "kg": 1.0, "deg": 1.0, "dim": 1.0, "rpm": 1.0,
         "lbf/in": LBF / IN, "N/m": 1.0, "ft-lbf/deg": LBF * 0.3048 * 180 / np.pi,
         "ft-lbf": LBF * 0.3048, "N*m": 1.0}


def _rows(ws):
    return [r for r in ws.iter_rows(values_only=True) if any(c is not None for c in r)]


def _table(ws):
    """Sheet with a header row and a units row -> {key: {column: SI value}}."""
    rows = _rows(ws)
    head, units = rows[0], rows[1]
    out = {}
    for r in rows[2:]:
        out[str(r[0]).strip()] = {str(h).strip(): float(v) * TO_SI[str(u).strip()]
                                  for h, u, v in zip(head[1:], units[1:], r[1:]) if v is not None}
    return out


def read(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    const = {}
    for r in _rows(wb["Constants"])[1:]:
        if r[1] and r[2] is not None:
            const[str(r[1])] = float(r[2]) * TO_SI.get(str(r[3]).strip(), 1.0)
    dyno = _rows(wb["Dyno Data"])
    du = [TO_SI[str(u)] for u in dyno[1]]
    dyno = np.array([[float(a) * du[0], float(b) * du[1]] for a, b in dyno[2:]])
    head, *configs = _rows(wb["Configurations"])
    configs = {str(r[0]): dict(zip(head, r)) for r in configs}
    return dict(primary=_table(wb["Primary Spring"]), secondary=_table(wb["Secondary Spring"]),
                flyweights=_table(wb["Flyweights"]), cam=_table(wb["Secondary Cam"]),
                constants=const, dyno=dyno, configurations=configs)


def configuration(path, name):
    """Overrides for one named configuration. Spring preload = rate x (minimum preload
    + primary initial travel xPi), as in CVT_model_v4 (rStateFix/wStateFix)."""
    w = read(path)
    c, k = w["configurations"][name], w["constants"]
    ps = w["primary"][str(c["Primary Spring"])]
    ss = w["secondary"][str(c["Secondary Spring"])]
    fw = w["flyweights"][str(c["Flyweights"]).replace(".0", "")]
    cam = w["cam"][str(c["Secondary Cam"])]
    k_p = ps["SPRING RATE"]
    return {
        "primary.k_spring": k_p,
        "primary.F_preload": k_p * (ps["Minimum Preload"] + k["xPi"]),
        "primary.m_weight": fw["Weight"],
        "primary.m_cover": k["mC"],
        "primary.n_weights": int(k["nFly"]),
        "primary.r_w0": k["rFly"],
        "primary.kappa": 1.0 / (1 / np.tan(np.radians(k["thetaP"])) + 1 / np.tan(np.radians(k["alpha"]))),
        "secondary.k_linear": ss["LINEAR SPRING RATE"],
        "secondary.k_torsion": ss["TORSIONAL SPRING RATE"],
        "secondary.torsion_preload_deg": float(c["Secondary Pretension"]),
        "secondary.F_preload_equiv": ss["LINEAR SPRING RATE"] * ss["Minimum Linear Preload"]
        + ss["TORSIONAL SPRING RATE"] * np.radians(float(c["Secondary Pretension"]))
        / (cam["Ramp Radius"] * np.tan(np.radians(cam["Ramp Aangle"]))),
        "secondary.r_cam": cam["Ramp Radius"],
        "secondary.cam_angle_deg": cam["Ramp Aangle"],
        "sheave.beta_deg": k["theta"],
        "belt.mu": k["muT"],
        "vehicle.final_drive": k["GR"],
        "vehicle.tire_radius": k["dTire"] / 2,
        "geometry.belt_pitch_length": k["L"],
        "geometry.center_distance": k["xc"],
        "geometry.rp_mech_min": k["rPmin"], "geometry.rp_mech_max": k["rPmax"],
        "geometry.rs_mech_min": k["rSmin"], "geometry.rs_mech_max": k["rSmax"],
        "engine.rpm": [float(v) for v in w["dyno"][:, 0]],
        "engine.torque": [float(v) for v in w["dyno"][:, 1]],
    }
