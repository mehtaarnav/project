"""Coupled launch model (cvtsim.cvt.system): energy balance, idle behaviour, time-step convergence."""
import numpy as np
import pytest

from cvtsim.io.params import load
from cvtsim.cvt.system import build, relax, simulate, CH


def _run(t_end, **kw):
    m = build(load(), **kw)
    relax(m)
    return simulate(m, t_end, nout=200)[0]


@pytest.fixture(scope="module")
def short_run():
    return _run(0.15)


def col(r, k):
    return r[:, CH.index(k)]


def test_energy_balance(short_run):
    E, Win, Wd = col(short_run, "E_mech"), col(short_run, "W_in"), col(short_run, "W_diss")
    resid = (E - E[0]) - (Win - Wd)
    assert np.max(np.abs(resid)) < 1e-4 * Win[-1]


def test_no_drive_below_engagement(short_run):
    rpm = col(short_run, "w_e") * 30 / np.pi
    below = rpm < 1700.0  # derived engagement 1725.7 rpm
    assert below.sum() > 5
    assert np.all(col(short_run, "Fa_p")[below] == 0.0)
    assert np.max(np.abs(col(short_run, "w_s")[below])) < 0.05


def test_time_step_convergence(short_run):
    half = _run(0.15, dt=load()["sim"]["dt"] / 2)
    for k in ("w_e", "w_s", "x_p", "Fa_p"):
        a, b = col(short_run, k)[-1], col(half, k)[-1]
        assert abs(a - b) <= 1e-3 * max(abs(b), 1e-3), (k, a, b)
