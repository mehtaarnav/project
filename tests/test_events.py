"""Events on the fast model: docs/theory/16_events.md section 5 (T16.1-T16.4)."""
import math

import numpy as np
import pytest

from cvtsim import events
from cvtsim.drivetrain.vehicle import G_ACC
from cvtsim.fast.model import FCH, build_fast, simulate_fast
from cvtsim.io.params import load


def test_t16_1_course_sampling():
    s, g, tr, R = events.course(dict(segments=[[10.0, 0.0, 0.0, 0.0, 1.0], [10.0, 0.0, 10.0, 0.0, 0.8]]))
    assert len(s) == 21 and s[-1] == pytest.approx(20.0)
    assert math.degrees(np.interp(15.0, s, g)) == pytest.approx(5.0)
    assert tr[-1] == 0.8 and tr[5] == 1.0


def test_t16_2_target_speed():
    s = np.arange(0.0, 51.0)
    R = np.where(s >= 40, 10.0, 0.0)
    v = events.target_speed(s, np.zeros_like(s), np.ones_like(s), R, 0.7, 0.5)
    vc = math.sqrt(0.7 * G_ACC * 10.0)
    assert v[45] == pytest.approx(vc)
    assert v[20] == pytest.approx(math.sqrt(vc ** 2 + 2 * 0.5 * 20.0))


def test_t16_3_acceleration_matches_time_scenario():
    m, res, _ = events.run_event("default", "acceleration")
    r = simulate_fast(build_fast(load(tune="default")), 8.0, nout=80001)
    s = np.cumsum(np.r_[0.0, 0.5 * (r[1:, FCH.index("v")] + r[:-1, FCH.index("v")]) * np.diff(r[:, 0])])
    t_cross = np.interp(45.72, s, r[:, 0])
    assert m["finished"] and m["time_s"] == pytest.approx(t_cross, abs=2e-3)


def test_t16_4_endurance_energy_and_finish():
    m, res, _ = events.run_event("default", "endurance")
    c = lambda k: res[:, FCH.index(k)]
    E, Wi, Wd = c("E_kin")[:-1], c("W_in")[:-1], c("W_diss")[:-1]
    assert m["finished"] and m["ratio_min"] == pytest.approx(0.8, abs=1e-6)
    assert np.max(np.abs((E - E[0]) - (Wi - Wd))) / Wi[-1] < 1e-3


def test_every_course_loads():
    for name in events.track_names():
        s, g, tr, R = events.course(events.load_track(name))
        assert s[-1] > 0 and np.all(np.diff(s) > 0)
