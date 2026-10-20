"""Coupled launch model (cvtsim.cvt.system): energy balance, idle behaviour, time-step convergence."""
import numpy as np
import pytest

from cvtsim.io.params import load
from cvtsim.cvt.system import build, relax, simulate, CH


def _run(t_end, **kw):
    m = build(load(), **kw)
    relax(m)
    return simulate(m, t_end, nout=400)[0]


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


@pytest.mark.slow
def test_time_step_convergence():
    """dt vs dt/2 after engagement: means over 0.6-0.8 s agree to 1e-3.

    Engagement is a sheave impact on the slack belt followed by chatter; single samples during
    it are not dt-convergent (a sequence of impacts), and the settled sheave position x_p is
    path-dependent to about 0.2 % (7 um) until clutch friction is modelled (REPORT, known
    limitations), so x_p is not part of this check. Speeds and forces are.
    """
    dt = load()["sim"]["dt"]
    full, half = _run(0.8, dt=dt), _run(0.8, dt=dt / 2)
    w = full[:, 0] > 0.6
    for k in ("w_e", "w_s", "v", "Fa_p", "Fa_s", "T_tight"):
        a, b = col(full, k)[w].mean(), col(half, k)[w].mean()
        assert abs(a - b) <= 1e-3 * abs(b), (k, a, b)
