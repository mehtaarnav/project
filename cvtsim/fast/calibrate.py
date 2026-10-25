"""Calibrate the fast CVT model against a full-model run (docs/theory/15 section 5).

    python -m cvtsim.fast.calibrate <full runs to fit ...> --validate <full runs to check ...>

Fits the constants in FIT by least squares on engine rpm (after engagement), primary travel and
vehicle speed of every fit run, writes cvtsim/fast/calibration.yaml and reports the section 5
acceptance checks on the fit runs and on the validation runs (not used in the fit).
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import least_squares

from cvtsim.cvt.system import CH
from cvtsim.fast.model import FCH, build_fast, simulate_fast
from cvtsim.io.params import load

CAL_FILE = Path(__file__).with_name("calibration.yaml")
FIT = ("K", "g1", "g2", "c_s", "c_s2")
X0 = np.array([np.log(1e-5), 0.5, 0.0, 0.87, 1.0])
LO = np.array([np.log(1e-8), -2.0, -2.0, -0.5, -2.0])
HI = np.array([np.log(1e-2), 5.0, 5.0, 4.0, 4.0])


def load_full(run):
    """(t, rpm, x_p [mm], v [m/s], scenario) of a full-model run directory."""
    run = Path(run)
    with gzip.open(run / "series.csv.gz", "rt") as f:
        cols = f.readline().strip().split(",")
        res = np.loadtxt(f, delimiter=",")
    c = lambda k: res[:, cols.index(k)]
    cfg = yaml.safe_load((run / "config.yaml").read_text())
    sc = (cfg.get("parameters") or {}).get("scenario")
    tune = json.loads((run / "meta.json").read_text())["tags"][0]
    full = dict(t=c("t"), rpm=c("w_e") * 30 / np.pi, x_p=c("x_p") * 1e3, v=c("v"),
                engaged=c("Fa_p") > 1.0, scenario=sc, tune=tune, name=run.name)
    full["p"] = load(tune=tune)
    sel = np.zeros(len(full["t"]), bool)
    sel[::40] = True                       # 10 ms spacing from 4 kHz output
    full["mask"] = sel & (full["t"] > full["t"][np.argmax(full["engaged"])])
    return full


def cal_from(x):
    return dict(K=float(np.exp(x[0])), g1=float(x[1]), g2=float(x[2]), c_s=float(x[3]), c_s2=float(x[4]))


def run_on(full, cal, p=None):
    p = p or full["p"]
    fm = build_fast(p, scenario=full["scenario"], cal=cal)
    r = simulate_fast(fm, float(full["t"][-1]), nout=int(full["t"][-1] / 1e-2) + 1)
    c = lambda k: np.interp(full["t"], r[:, 0], r[:, FCH.index(k)])
    return dict(rpm=c("w_e") * 30 / np.pi, x_p=c("x_p") * 1e3, v=c("v"), F_b=c("F_b"))


def residuals(x, fulls):
    out = []
    for full in fulls:
        f, m = run_on(full, cal_from(x)), full["mask"]
        out += [(f["rpm"] - full["rpm"])[m] / np.mean(full["rpm"][m]), (f["x_p"] - full["x_p"])[m] / 1.5,
                (f["v"] - full["v"])[m] / max(full["v"][-1], 1.0)]
    return np.concatenate(out)


def _clean(d):
    return {k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in d.items()}


def checks(full, f):
    """Section 5 acceptance checks; returns dict of values and pass flags."""
    eng = np.argmax(full["engaged"])
    m = np.arange(len(full["t"])) > eng
    rpm_err = np.sqrt(np.mean((f["rpm"] - full["rpm"])[m] ** 2)) / np.mean(full["rpm"][m])
    xp_err = np.sqrt(np.mean((f["x_p"] - full["x_p"])[m] ** 2))
    v_err = abs(f["v"][-1] - full["v"][-1]) / full["v"][-1]
    i_f = np.argmax(f["x_p"] >= 0.0)
    eng_err = abs(f["rpm"][i_f] - full["rpm"][eng]) / full["rpm"][eng]
    return dict(rpm_rms=rpm_err, rpm_ok=rpm_err < 0.05, xp_rms_mm=xp_err, xp_ok=xp_err < 1.5,
                v_end=v_err, v_ok=v_err < 0.03, engage=eng_err, engage_ok=eng_err < 0.03)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.fast.calibrate")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--validate", nargs="*", default=[])
    a = ap.parse_args(argv)
    fulls = [load_full(r) for r in a.runs]
    fit = least_squares(residuals, X0, bounds=(LO, HI), args=(fulls,), x_scale=0.5)
    cal = cal_from(fit.x)
    out = dict(calibration=cal, cost=float(fit.cost),
               fit={f["name"]: _clean(checks(f, run_on(f, cal))) for f in fulls},
               validation={v["name"]: _clean(checks(v, run_on(v, cal))) for v in map(load_full, a.validate)})
    CAL_FILE.write_text(yaml.safe_dump(out, sort_keys=False))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
