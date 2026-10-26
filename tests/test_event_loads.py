"""Event loads and replays: docs/theory/17_event_loads.md section 4."""
import math

import numpy as np
import pytest

from cvtsim import event_loads as el


def test_t17_1_2_hub_forces():
    Fx, Fy = el.hub_forces(400.0, 300.0, 0.06, 0.06, 0.232)
    assert Fx == pytest.approx(800.0) and Fy == pytest.approx(0.0)
    Fx, Fy = el.hub_forces(400.0, 300.0, 0.03, 0.03 + 0.116, 0.232)
    assert Fx == pytest.approx(math.sqrt(3) * 400.0) and Fy == pytest.approx(150.0)


def test_t17_3_pick_windows():
    t = np.arange(0.0, 20.0, 0.1)
    y = np.zeros_like(t)
    y[[10, 12, 60, 150]] = [5.0, 9.0, 7.0, 6.0]
    assert el.pick_windows(t, y, 3) == pytest.approx([1.2, 6.0, 15.0])


def test_t17_5_sweep_final_drive():
    rows = el.sweep(["default"], "acceleration", [8.1, 9.0])
    assert [r["final_drive"] for r in rows] == [8.1, 9.0]
    assert rows[1]["rpm_mean"] > rows[0]["rpm_mean"]


def test_event_durability_fields():
    from cvtsim.events import run_event
    from cvtsim.io.params import load
    p = load(tune="default")
    _, res, _ = run_event("default", "acceleration")
    d, sh = el.event_durability(res, p, hours=1)
    assert d["pto_side_peak_N"] > 0 and "durability_sampled_ok" not in d
    assert d["damage_sec_per_duty"] == pytest.approx(d["damage_sec"] * 3600 / d["event_time_s"])


@pytest.mark.slow
def test_t17_4_warm_start_matches_uninterrupted_run():
    import gzip
    from pathlib import Path
    from cvtsim.cvt.system import CH, build, simulate, warm_start
    from cvtsim.io.params import load
    run = Path("runs/20261025-201508-tune-cvtech_0651_1100/series.csv.gz")
    if not run.exists():
        pytest.skip("full-model series not present")
    with gzip.open(run, "rt") as f:
        f.readline()
        r = np.loadtxt(f, delimiter=",")
    c = lambda k: r[:, CH.index(k)]
    t = c("t")
    i0 = int(np.argmin(abs(t - 4.0)))
    m = build(load(tune="cvtech_0651_1100"))
    warm_start(m, c("w_e")[i0], c("x_p")[i0], c("x_s")[i0], c("w_s")[i0], c("v")[i0], c("T_belt")[i0],
               t0=t[i0], theta_e=c("theta_e")[i0])
    res, _ = simulate(m, 0.8, nout=801, n_frames=2)
    g = lambda k: res[:, CH.index(k)]
    w, o = g("t") > t[i0] + 0.2, (t > t[i0] + 0.2) & (t <= t[i0] + 0.8)
    for a, b in ((np.abs(g("T_pto"))[w], np.abs(c("T_pto"))[o]), (g("T_tight")[w], c("T_tight")[o]),
                 (np.hypot(g("Fx_p"), g("Fy_p"))[w], np.hypot(c("Fx_p"), c("Fy_p"))[o])):
        assert a.max() == pytest.approx(b.max(), rel=0.05)
