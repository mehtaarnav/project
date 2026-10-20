import json

import numpy as np

from cvtsim.io import params
from cvtsim.viz.export import write_viewer, _sig

CH = ["t", "w_e", "x_p"]


def _data(n=5000, F=7, N=12):
    rng = np.random.default_rng(0)
    t = np.linspace(0, 8, n)
    res = np.column_stack([t, 300 + 100 * rng.random(n), 0.01 * rng.random(n)])
    ang = np.linspace(0, 2 * np.pi, N, endpoint=False)
    snaps = dict(t=np.linspace(0, 8, F),
                 belt=np.stack([np.column_stack([np.cos(ang), np.sin(ang)]) * 0.123456789] * F),
                 tension=rng.random((F, N)) * 1234.56789)
    for k in ("x_p", "x_s", "theta_e", "theta_s", "r_w"):
        snaps[k] = rng.random(F) * 0.0123456789
    return res, snaps


def test_sig_rounding():
    r = _sig(np.array([123456.789, 0.000123456789, 0.0, -9.87654321]))
    assert np.allclose(r, [123460.0, 0.00012346, 0.0, -9.8765], rtol=0, atol=1e-12)


def test_write_viewer_accepts_shaft_histories(tmp_path):
    """Contract with the real producer: shaft_histories output goes straight into write_viewer."""
    from cvtsim.cvt.system import CH as SIM_CH, build, relax, simulate
    from cvtsim.fe.shaft_loads import shaft_histories
    p = params.load()
    m = build(p)
    relax(m)
    res, snaps = simulate(m, 0.01, nout=50, n_frames=5)
    sh = shaft_histories(res, p, SIM_CH)
    d = json.loads(write_viewer(tmp_path, res, SIM_CH, snaps, p, sh).read_text())
    assert set(d["shafts"]) == set(sh["cols"])
    assert len(d["shafts"]["t"]) == len(d["series"]["t"])
    mech = d["mechanism"]
    assert all(len(v) == len(d["frames"]["t"]) for v in mech.values())
    a0 = mech["helix_angle_deg"][0]
    assert abs(a0 - 45.0) < 1e-3                                   # straight 45 deg helix
    assert abs(mech["torsion_twist"][0] - (np.pi / 6 + mech["helix_twist"][0])) < 1e-4
    assert abs(mech["spring_p_compression"][0] - d["series"]["F_spr_p"][0] / 13208.0) < 1e-6
    assert "peak_belt_tension" in d["events"]


def test_write_viewer(tmp_path):
    p = params.load()
    res, snaps = _data()
    sh = {"t": np.linspace(0, 8, 300), "F_brg": np.arange(300) * 1.234567}
    d = json.loads(write_viewer(tmp_path, res, CH, snaps, p, sh).read_text())
    assert set(d) == {"meta", "geometry", "frames", "series", "shafts"}
    assert d["meta"]["dt"] == 1.25e-7 and d["meta"]["tune"] == "default"
    g = d["geometry"]
    assert g["C"] == 0.232 and g["beta_deg"] == 12 and g["n_weights"] == 3
    assert g["spring_s"]["F0"] == 381.62 and g["spring_p"]["k"] == 13208
    assert abs(g["x_p_stop"] - 0.0175241) < 1e-6           # 05 T5.1b, secondary-limited
    gap = 0.0089662                                          # travel from full open to belt contact
    assert abs(g["ramp"][0][0] + gap) < 1e-9 and abs(g["ramp"][0][1] - (0.04445 - 0.5 * gap)) < 1e-6  # 5 s.f. export rounding
    assert abs(g["ramp"][-1][1] - (0.04445 + 0.5 * 0.0175241)) < 1e-5
    assert abs(g["helix"][-1][1] - 0.02448 / 0.0381) < 1e-3
    s = d["series"]
    assert len(s["t"]) <= 2000 and len(s["t"]) == len(s["w_e"]) == len(s["x_p"])
    assert s["t"][0] == 0 and s["t"][-1] == 8
    fr = d["frames"]
    assert len(fr["belt"]) == 7 and len(fr["belt"][0]) == 12 and len(fr["belt"][0][0]) == 2
    assert fr["belt"][0][0][0] == 0.12346
    assert all(float(f"{v:.5g}") == v for v in fr["tension"][0])
    assert len(d["shafts"]["t"]) == 300 and d["shafts"]["F_brg"][1] == 1.2346


def test_short_series_not_decimated(tmp_path):
    p = params.load()
    res, snaps = _data(n=50)
    d = json.loads(write_viewer(tmp_path, res, CH, snaps, p).read_text())
    assert len(d["series"]["t"]) == 50 and d["shafts"] == {}
