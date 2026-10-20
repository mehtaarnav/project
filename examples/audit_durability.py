"""Independent audit of one launch's durability numbers (docs/theory/12 section 7).

Recomputes, without cvtsim.post.durability or cvtsim.fe.fatigue/bearings, from a run's series:
revolutions, Kf, Se, per-revolution bending damage, the simultaneous (time-local) Goodman factor
next to the peak-envelope factor, and one bearing's equivalent load and required rating.
Inputs M(t), T(t), R(t) come from cvtsim.fe.shaft_loads (verified by spec 07 tests).

    python examples/audit_durability.py runs/<tune-default run>
"""
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

from cvtsim.cvt.system import CH
from cvtsim.fe.shaft_loads import shaft_histories
from cvtsim.io.params import load

run = Path(sys.argv[1])
with gzip.open(run / "series.csv.gz", "rt") as f:
    cols = f.readline().strip().split(",")
    res = np.loadtxt(f, delimiter=",")
assert cols == CH, "series columns differ from this checkout's CH"
meta = json.loads((run / "meta.json").read_text())
p = load(tune="default")
sh = shaft_histories(res, p, CH)
mats = yaml.safe_load(open("cvtsim/io/materials.yaml"))
c = lambda k: res[:, CH.index(k)]
t = c("t")
out = dict(run=run.name, run_commit=meta.get("commit"), dt_out_s=float(np.median(np.diff(t))))

for name, tag, w in (("primary", "prim", c("w_c")), ("secondary", "sec", c("w_s"))):
    sd, cfg = p["shafts"][name], p["durability"][name]
    m = mats["materials"][cfg["material"]]
    Su, Sy = m["Su_MPa"]["value"], m["Sy_MPa"]["value"]          # MPa
    d, r = sd["critical"]["d"], sd["critical"]["r"]
    Ktb, Kts = sh["kt"][name]
    a = 0.0254 * (2070.0 / Su) ** 1.8                              # mm, 12 sec 2.3
    q = 1.0 / (1.0 + a / (r * 1e3))
    Kf, Kfs = 1.0 + q * (Ktb - 1.0), Kts
    ka = {"ground": (1.58, -0.085), "machined": (4.51, -0.265)}[cfg["finish"]]
    Se = ka[0] * Su ** ka[1] * 1.24 * (d * 1e3) ** -0.107 * 0.5 * Su
    bexp = math.log10(Se / (0.9 * Su)) / 3.0
    A = 0.9 * Su * 1e3 ** (-bexp)
    zb, zt = 32.0 / (math.pi * d ** 3) / 1e6, 16.0 / (math.pi * d ** 3) / 1e6   # MPa per N m
    M, T = sh[f"M_{tag}"], sh[f"T_{tag}"]
    nrev, D = 0.0, 0.0
    worst = (math.inf, None)
    for i in range(len(t) - 1):
        dn = 0.5 * (abs(w[i]) + abs(w[i + 1])) * (t[i + 1] - t[i]) / (2 * math.pi)
        nrev += dn
        sa = Kf * zb * 0.5 * (abs(M[i]) + abs(M[i + 1]))
        sm = math.sqrt(3) * Kfs * zt * abs(0.5 * (T[i] + T[i + 1]))
        seq = sa / (1.0 - sm / Su)
        if seq >= Se:
            D += dn / (seq / A) ** (1.0 / bexp)
        ng = 1.0 / (sa / Se + sm / Su) if sa + sm > 0 else math.inf
        if ng < worst[0]:
            worst = (ng, float(t[i]))
    env_sa = Kf * zb * np.abs(M).max()
    env_sm = math.sqrt(3) * Kfs * zt * np.abs(T).max()
    out[name] = dict(Su=Su, Kf=Kf, Kfs=Kfs, Se_MPa=Se, N_rev=nrev, D_bend_per_launch=D,
                     n_goodman_simultaneous_min=worst[0], at_t_s=worst[1],
                     n_goodman_envelope_no_torque_range=1.0 / (env_sa / Se + env_sm / Su),
                     t_of_peak_M=float(t[np.argmax(np.abs(M))]), t_of_peak_T=float(t[np.argmax(np.abs(T))]),
                     peak_M_Nm=float(np.abs(M).max()), peak_T_Nm=float(np.abs(T).max()))

# one bearing: PTO-side main (primary bearing 2)
R, w = sh["F_bearing_prim2"], c("w_c")
dn = 0.5 * (np.abs(w[1:]) + np.abs(w[:-1])) * np.diff(t) / (2 * math.pi)
Rm = 0.5 * (R[1:] + R[:-1])
Peq = (np.sum(Rm ** 3 * dn) / dn.sum()) ** (1 / 3)
nL = p["durability"]["target_launches"]
out["bearing_prim2"] = dict(P_eq_N=float(Peq), N_rev=float(dn.sum()),
                            C_req_N_for_target_launches=float(Peq * (nL * dn.sum() / 1e6) ** (1 / 3)),
                            C_listed_UNVERIFIED_N=mats["bearings"]["ch440_main_pto_side"]["C_N"]["value"])
print(json.dumps(out, indent=2))
