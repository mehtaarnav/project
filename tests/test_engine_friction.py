"""Spec 11 targets T11.1-T11.6 (docs/theory/11_engine_friction_uq.md, sections 1-3, 5).

Numbers come from docs/theory/targets.py (closed forms, no cvtsim import). The mean-torque
equivalence compares against a launch stored before spec 11 (commit 3eea3d6, tests/data)."""
import gzip
import math
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from cvtsim.io.params import load
from cvtsim.cvt.system import build, relax, simulate, CH, ST, engine_torque, helix_friction
from cvtsim.drivetrain.pulses import crank_kinematics, gas_torque, rec_torque, pulse_table

ROOT = pathlib.Path(__file__).resolve().parents[1]
RPM = 2.0 * math.pi / 60.0
# spec 11 off: mean engine torque, one rigid rotor, frictionless clutches
MEAN = {"engine.pulses": False, "engine.pto_model": "rigid", "primary.mu_roller": 0.0,
        "secondary.mu_helix": 0.0, "sheave.mu_bushing": 0.0}


@pytest.fixture(scope="module")
def tg():
    out = subprocess.run([sys.executable, str(ROOT / "docs/theory/targets.py")],
                         capture_output=True, text=True, check=True).stdout
    return {l.split()[0]: float(l.split()[1]) for l in out.splitlines() if l.strip()}


@pytest.fixture(scope="module")
def eng():
    return load()["engine"]


def _geom(e):
    r = 0.5 * e["stroke"]
    return r, r / e["rod_length"]


def test_T11_1_slider_crank(tg, eng):
    r, lam = _geom(eng)
    assert crank_kinematics(math.pi, r, lam)[0] == pytest.approx(tg["t11_1_s_180"], rel=1e-15)
    assert crank_kinematics(0.0, r, lam)[0] == 0.0
    th = np.linspace(0.0, 2 * np.pi, 721)
    d = 1e-5
    s, sp, spp = crank_kinematics(th, r, lam)
    fd1 = (crank_kinematics(th + d, r, lam)[0] - crank_kinematics(th - d, r, lam)[0]) / (2 * d)
    fd2 = (crank_kinematics(th + d, r, lam)[1] - crank_kinematics(th - d, r, lam)[1]) / (2 * d)
    assert np.max(np.abs(fd1 - sp)) <= 1e-9 * np.max(np.abs(sp))
    assert np.max(np.abs(fd2 - spp)) <= 1e-9 * np.max(np.abs(spp))


def test_T11_2_pulse_mean_equals_map(tg, eng):
    """Cycle mean of T_gas + T_rec at constant speed = governed map torque. Sampling at the table
    nodes and midpoints integrates the piecewise-linear table exactly."""
    p = load()
    tab = pulse_table(eng)
    gm, gq = tab["gm"], tab["gq"]
    r, lam = _geom(eng)
    th = 4 * np.pi * np.arange(2 * len(gm)) / (2 * len(gm))
    w_rpm, T_rpm = np.asarray(eng["rpm"], float), np.asarray(eng["torque"], float)
    for rpm in (2400, 3000, 3600):
        w = rpm * RPM
        T_m, _ = engine_torque(w, w_rpm, T_rpm, eng["rpm_cut"], eng["governor_band"])
        T = np.array([gas_torque(a, T_m, gm, gq) + rec_torque(a, w, eng["m_rec"], r, lam)[0] for a in th])
        assert abs(T.mean() - tg[f"t11_2_Tmap_{rpm}"]) <= 1e-9 * tg[f"t11_2_Tmap_{rpm}"], rpm
        assert T.max() > 5 * T_m  # a real pulse, not the mean
    assert p["engine"]["pulses"]


def test_T11_2_table_mean_converges(tg, eng):
    """Table resolution: the trapezoid mean of the unnormalized gas torque per joule converges as
    h^2 to the continuum first-law value and is within 1e-9 at pulse_table_n."""
    ref = tg["t11_gas_mean_per_J"]
    n = int(eng["pulse_table_n"])
    err = [abs(pulse_table(eng, n // 2 ** k)["mean_raw"] / ref - 1) for k in (2, 1, 0)]
    assert err[-1] < 1e-9
    assert 3.5 < err[0] / err[1] < 4.5 and 3.5 < err[1] / err[2] < 4.5


def test_T11_3_reciprocating_torque_zero_mean(eng):
    r, lam = _geom(eng)
    th = 2 * np.pi * np.arange(3600) / 3600
    w = 3600 * RPM
    T = rec_torque(th, w, eng["m_rec"], r, lam)[0]
    assert abs(T.mean()) <= 1e-12 * np.abs(T).max()


def _launch(p, t_end, nout=400):
    m = build(p)
    relax(m)
    return simulate(m, t_end, nout=nout)[0]


def test_T11_4_energy_closure_pulses_pto_friction():
    """Default model (pulses, PTO spring and damper, helix, roller and bushing friction), 0.5 s."""
    res = _launch(load(), 0.5)
    c = lambda k: res[:, CH.index(k)]
    E, Win, Wd = c("E_mech"), c("W_in"), c("W_diss")
    assert np.max(np.abs((E - E[0]) - (Win - Wd))) < 1e-4 * Win[-1]
    # the new terms are active
    assert np.ptp(c("T_pulse")) > 100.0
    assert np.max(np.abs(c("w_e") - c("w_c"))) > 0.1
    assert np.max(np.abs(c("F_fric_p"))) > 10.0 and np.max(np.abs(c("F_fric_s"))) > 10.0
    assert np.max(c("T_pto")) > np.max(c("T_eng"))


def test_T11_5_pto_torsional_frequency(tg):
    """Free vibration of crank and clutch on the PTO spring: no engine torque, no damping, the
    primary open at its stop (flyweights at r_w(-g) join the clutch inertia), 1 mrad twist."""
    p = load(overrides={"engine.torque": [0.0] * 7, "engine.pulses": False,
                        "engine.pto_damping_ratio": 0.0, "engine.idle_rpm": 0.0})
    assert p["derived"]["pto"]["k"] == pytest.approx(tg["t11_5_k_pto"], rel=1e-12)
    assert p["derived"]["pto"]["f_n"] == pytest.approx(tg["t11_5_f_spec"], rel=1e-12)
    m = build(p)
    relax(m)
    m["st"][ST.index("theta_e")] = 1e-3
    res = simulate(m, 0.05, nout=5000)[0]
    t, T = res[:, 0], res[:, CH.index("T_pto")]
    i = np.where((T[:-1] < 0.0) & (T[1:] >= 0.0))[0]
    tz = t[i] - T[i] * (t[i + 1] - t[i]) / (T[i + 1] - T[i])
    f = (len(tz) - 1) / (tz[-1] - tz[0])
    assert len(tz) >= 8
    assert f == pytest.approx(tg["t11_5_f_open_clutch"], rel=1e-3)
    # >= 20 steps per period at the launch step
    assert 1.0 / (f * p["sim"]["dt"]) >= 20


def test_T11_6_helix_friction_quasistatic(tg):
    """45 deg helix at x = 0, tau_mov = 20 N m: closing force h' tau_c - Q_f while opening and closing."""
    p = load()
    s, d = p["secondary"], p["derived"]
    h, hp, _ = d["helix"](0.0)
    tau_c = 20.0 + s["k_torsion"] * (d["theta_t0"] + h)
    mu, v_x = s["mu_helix"], p["sheave"]["v_fric"]
    F = {sgn: hp * tau_c - helix_friction(tau_c, hp, s["r_cam"], mu, sgn * 100 * v_x, v_x) for sgn in (1, -1)}
    assert F[1] == pytest.approx(tg["t11_6_F_close_opening"], rel=1e-12)
    assert F[-1] == pytest.approx(tg["t11_6_F_close_closing"], rel=1e-12)
    assert F[1] / F[-1] == pytest.approx(tg["t11_6_ratio"], rel=1e-12)
    assert helix_friction(tau_c, hp, s["r_cam"], mu, 0.0, v_x) == 0.0


def test_mean_torque_model_reproduced():
    """pulses off, rigid PTO, frictionless clutches: the 0.15 s launch from before spec 11
    (commit 3eea3d6), every stored channel, to 1e-9 of its range."""
    with gzip.open(ROOT / "tests/data/prespec11_launch_0p15s.csv.gz", "rt") as f:
        f.readline()
        names = f.readline().strip().split(",")
        ref = np.loadtxt(f, delimiter=",")
    res = _launch(load(overrides=MEAN), 0.15, nout=200)
    assert res.shape[0] == ref.shape[0]
    for j, k in enumerate(names):
        a, b = res[:, CH.index(k)], ref[:, j]
        assert np.abs(a - b).max() <= 1e-9 * max(np.abs(b).max(), 1e-12), k
    c = lambda k: res[:, CH.index(k)]
    assert np.array_equal(c("w_c"), c("w_e")) and np.array_equal(c("theta_c"), c("theta_e"))
    assert np.all(c("F_fric_p") == 0.0) and np.all(c("F_fric_s") == 0.0)


@pytest.mark.slow
def test_pto_rigid_limit():
    """Mean torque, PTO stiffness k x {1, 4, 16} against the rigid rotor (0.15 s launch): the
    largest crank-speed difference falls with k (11 section 2.1)."""
    rigid = _launch(load(overrides={"engine.pulses": False, "engine.pto_model": "rigid"}), 0.15, 300)
    e = load()["engine"]
    L = e["pto_length"] + e["crank_web_length"]
    err = []
    for s in (1.0, 4.0, 16.0):  # k ~ 1/L_eff; the crank-web term absorbs the shortening
        r = _launch(load(overrides={"engine.pulses": False,
                                    "engine.crank_web_length": L / s - e["pto_length"]}), 0.15, 300)
        a, b = r[:, CH.index("w_e")], rigid[:, CH.index("w_e")]
        err.append(np.abs(a - b).max() / np.abs(b).max())
    assert err[0] > err[1] > err[2]
    assert err[2] < 1e-3
