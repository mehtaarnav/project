"""Spec 09 targets T9.1-T9.9 (docs/theory/09_tuning_vehicle_belt.md, section 7).

Numbers come from docs/theory/targets.py (closed forms, no cvtsim import), read here by
running the script. T9.9 compares against a stored pre-library series (tests/data)."""
import gzip
import math
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from cvtsim.io.params import component, load, tune_names
from cvtsim.cvt.profiles import Ramp
from cvtsim.cvt.primary import (engagement_omega, engagement_omega_profile, flyweight_force_profile,
                                negative_stiffness_omega_profile)
from cvtsim.cvt.secondary import secondary_static_force_helix
from cvtsim.cvt.system import (build, relax, simulate, CH, primary_accel, secondary_accel,
                               belt_temperature_step)
from cvtsim.drivetrain.steady import steady_speed_linear, simulate_fixed_ratio_tire
from cvtsim.drivetrain.vehicle import tire_force

ROOT = pathlib.Path(__file__).resolve().parents[1]
RPM = 2.0 * math.pi / 60.0
# "new vehicle terms off": 06 rigid driveline, lossless, no aero, no wheel inertia, 06 governor taper
OFF = {"vehicle.model": "rigid", "vehicle.eta": 1.0, "vehicle.J_wheels": 0.0, "vehicle.CdA": 0.0,
       "engine.governor_band": 200.0,
       # inputs frozen at commit 7268716, when the reference series was generated
       "geometry.belt_pitch_length": 0.88009, "geometry.rp_mech_min": 0.0315, "vehicle.Crr": 0.015,
       "sim.dt": 2.5e-7,
       # the reference referenced spring and ramp at full open with a 2 mm gap; restate them at contact
       "primary.engagement_gap": 0.002, "primary.F_preload": 683.04 + 13208.0 * 0.002,
       "primary.ramp": {"type": "straight", "r_w0": 0.04445 + 0.5 * 0.002, "kappa": 0.5},
       "engine.rpm": [2400, 2600, 2800, 3000, 3200, 3400, 3600],
       "engine.torque": [25.08, 24.54, 23.59, 22.51, 20.88, 19.66, 18.30]}


@pytest.fixture(scope="module")
def tg():
    out = subprocess.run([sys.executable, str(ROOT / "docs/theory/targets.py")],
                         capture_output=True, text=True, check=True).stdout
    return {l.split()[0]: float(l.split()[1]) for l in out.splitlines() if l.strip()}


def _rk4(rhs, y, h, n):
    for _ in range(n):
        k1 = rhs(y); k2 = rhs(y + 0.5 * h * k1); k3 = rhs(y + 0.5 * h * k2); k4 = rhs(y + h * k3)
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return y


def test_T9_1_straight_ramp_table_reproduces_05(tg):
    u = np.linspace(0.0, 0.0275, 12)
    p = load(overrides={"primary.ramp": {"type": "table", "table": np.column_stack([u, 0.04445 + 0.5 * u]).tolist()}})
    assert p["derived"]["ramp"].kind == "pchip"
    assert engagement_omega_profile(p) / RPM == pytest.approx(tg["engage_rpm_kappa0.5"], rel=1e-5)
    assert engagement_omega(p) / RPM == pytest.approx(tg["engage_rpm_kappa0.5"], rel=1e-5)
    for uu in (0.0, 0.007, 0.02):
        assert negative_stiffness_omega_profile(uu, p) / RPM == pytest.approx(tg["negstiff_rpm_kappa0.5"], rel=1e-5)


def test_T9_2_curved_ramp_flyweight_force(tg):
    p = load(overrides={"primary.ramp": {"type": "poly", "coeffs": [0.04445, 0.4, 4.0]}})
    assert flyweight_force_profile(0.01, 3000 * RPM, p) == pytest.approx(tg["t9_2_F_fw_u10mm_3000rpm"], rel=1e-12)


def test_T9_3_free_primary_curved_ramp_energy():
    """(P1'),(P2') with T_eng = 0, no belt, damping or stops, RK4 for 1 s. The curved ramp has
    negative net stiffness with the brief spring at any speed that balances its preload, so the
    test uses a stiffer test spring (k 40 kN/m, F_p0 568 N): a stable equilibrium near u = 10 mm
    keeps the motion inside the profile, where the ramp is smooth."""
    p = load(overrides={"primary.ramp": {"type": "poly", "coeffs": [0.04445, 0.4, 4.0]},
                        "primary.k_spring": 40000.0, "primary.F_preload": 568.0})
    pr, rk, rc = p["primary"], p["derived"]["ramp"].xk, p["derived"]["ramp"].c
    n, m_w, k, F0 = pr["n_weights"], p["derived"]["m_w"], pr["k_spring"], pr["F_preload"]
    J_e, m_ph = p["engine"]["J_e"], pr["m_axial"]

    def parts(y):
        w, u, ud = y
        return primary_accel(u, ud, w, 0.0, -(F0 + k * u), J_e, n, m_w, m_ph, rk, rc)

    def rhs(y):
        wd, udd = parts(y)[:2]
        return np.array([wd, y[2], udd])

    def energy(y):
        J, M = parts(y)[4:]
        return 0.5 * J * y[0] ** 2 + 0.5 * M * y[2] ** 2 + F0 * y[1] + 0.5 * k * y[1] ** 2

    y0 = np.array([2000.0 * RPM, 0.015, 0.0])
    us, y = [], y0
    for _ in range(100):
        y = _rk4(rhs, y, 1e-4, 100)
        us.append(y[1])
    assert 0.0 < min(us) and max(us) < p["primary"]["engagement_gap"] + p["primary"]["stroke"]
    assert max(us) - min(us) > 2e-3  # it moved
    assert abs(energy(y) - energy(y0)) / energy(y0) < 1e-6


def test_T9_4_two_angle_helix(tg):
    p = load(tune="progressive_helix")
    hel = p["derived"]["helix"]
    for xmm in (5, 15, 20):
        assert hel(xmm / 1e3)[0] == pytest.approx(tg[f"t9_4_h_{xmm}mm"], rel=1e-12)
    assert hel(0.015)[1] == pytest.approx(tg["t9_4_hprime_15mm"], rel=1e-12)
    assert secondary_static_force_helix(0.015, 20.0, p) == pytest.approx(tg["t9_4_Fs_15mm_20Nm"], rel=1e-12)


@pytest.mark.parametrize("helix", ["two_angle_45_35", "progressive_linear"])
def test_T9_5_free_secondary_progressive_helix_energy(helix):
    """Secondary 2x2 system with h'' terms, no damping or stops; a constant belt thrust F_a
    centres the motion inside the stroke, across the angle change. RK4 for 1 s."""
    spec = component("helixes", helix)
    p = load(overrides={"secondary.helix": spec})
    s, hel = p["secondary"], p["derived"]["helix"]
    J_s, J_m, m_m, k_l, k_t = s["J_shaft"], s["J_movable"], s["m_axial"], s["k_linear"], s["k_torsion"]
    F_l0, th0, F_a = p["derived"]["F_l0"], p["derived"]["theta_t0"], 528.0

    def parts(y):
        w, x, xd = y
        h, hp, hpp = hel(x)
        Q_x = F_a - F_l0 - k_l * x - k_t * (th0 + h) * hp
        return secondary_accel(xd, 0.0, Q_x, J_s, J_m, m_m, hp, hpp), h

    def rhs(y):
        (wd, xdd, *_), _ = parts(y)
        return np.array([wd, y[2], xdd])

    def energy(y):
        (_, _, a11, a12, a22), h = parts(y)
        w, x, xd = y
        U = (F_l0 - F_a) * x + 0.5 * k_l * x * x + 0.5 * k_t * (th0 + h) ** 2
        return 0.5 * (a11 * w * w + 2 * a12 * w * xd + a22 * xd * xd) + U

    y0 = np.array([300.0, 0.016, 0.0])
    xs, y = [], y0
    for _ in range(100):
        y = _rk4(rhs, y, 1e-4, 100)
        xs.append(y[1])
    assert max(xs) > 0.0105 and min(xs) < 0.0095 and min(xs) > 0.0  # crossed the angle change
    assert abs(energy(y) - energy(y0)) / energy(y0) < 1e-6


def test_T9_6_fixed_ratio_steady_speed(tg):
    v_q = steady_speed_linear(40.0, 0.05, 1.0, 8.1, 0.2921, 0.06, 250.0, 0.5 * 1.2 * 0.9, eta=0.9)
    assert v_q == pytest.approx(tg["t9_6_v_quadratic_noslip"], rel=1e-12)
    p = load(tune="default", overrides={"vehicle.Crr": 0.06, "vehicle.CdA": 0.9, "vehicle.eta": 0.9,
                                        "vehicle.mu_t": 0.7, "vehicle.s0": 0.08, "vehicle.b_over_L": 0.6})
    eng_w = np.array([0.0, 2000.0])
    r = simulate_fixed_ratio_tire(eng_w, 40.0 - 0.05 * eng_w, 0.05, 1.0, p, 60.0, 1e-3)
    assert r["v"][-1] == pytest.approx(tg["t9_6_v_with_slip"], rel=1e-4)
    assert r["s"][-1] == pytest.approx(tg["t9_6_tire_slip"], rel=1e-3)


def test_T9_7_tire_tracks_demand_then_spins(tg):
    p = load(tune="default")
    v = p["vehicle"]
    args = (v["mass"], v["tire_radius"], v["mu_t"], v["s0"], v["v_min"], v["b_over_L"], v["h_cg"] / v["wheelbase"])
    # below traction: wheel on a dyno at 10 m/s, axle torque 200 N m -> F_x = T/r_t, closed-form slip
    T, J, w, h = 200.0, v["J_wheels"], 10.0 / v["tire_radius"], 1e-3
    for _ in range(20000):
        F_x, N_r, s = tire_force(w, 10.0, 0.0, *args)
        w += h * (T - F_x * v["tire_radius"]) / J
    assert F_x == pytest.approx(T / v["tire_radius"], rel=1e-9)
    assert s == pytest.approx(tg["t9_7_slip_T200_v10"], rel=1e-9)
    assert N_r * v["mu_t"] == pytest.approx(tg["t9_7_traction_limit_N"], rel=1e-12)
    # above traction: the tire force a no-slip launch would need, m a with
    # a = eta T R G/r_t / (m + eta J (G/r_t)^2), exceeds mu_t N_r(a): the wheel spins
    T, R, J_e = 60.0, 3.0, p["engine"]["J_e"]
    G, r_t, m = v["final_drive"], v["tire_radius"], v["mass"]
    J = J_e * R * R + p["secondary"]["J_shaft"] + v["J_wheels"] / G ** 2
    a = v["eta"] * T * R * G / r_t / (m + v["eta"] * J * (G / r_t) ** 2)
    assert m * a > v["mu_t"] * m * (9.81 * v["b_over_L"] + a * v["h_cg"] / v["wheelbase"])
    r = simulate_fixed_ratio_tire([0.0, 1e4], [T, T], J_e, R, p, 2.0, 1e-4)
    assert r["s"].max() > 0.5
    assert np.all(np.abs(r["F_x"]) <= v["mu_t"] * r["N_r"] * (1 + 1e-12))


def test_T9_8_belt_heating_step_response(tg):
    th = load()["derived"]["belt_thermal"]
    assert th["C_b"] == pytest.approx(tg["t9_8_C_b"], rel=1e-12)
    assert th["hA"] == pytest.approx(tg["t9_8_hA"], rel=1e-12)
    assert th["tau"] == pytest.approx(tg["t9_8_tau"], rel=1e-12)
    P, n = 100.0, 1000
    T = th["T_amb"]
    for _ in range(n):
        T = belt_temperature_step(T, P, th["C_b"], th["hA"], th["T_amb"], th["tau"] / n)
    rise = (T - th["T_amb"]) / (P / th["hA"])
    assert rise == pytest.approx(tg["t9_8_rise_fraction_at_tau"], rel=1e-9)


def _launch(p, t_end=0.15):
    m = build(p)
    relax(m)
    return simulate(m, t_end, nout=200)


def test_T9_9_default_tune_reproduces_prelibrary_launch():
    """The default tune with the new vehicle terms off reproduces the pre-library (commit 7268716)
    launch series, every stored channel, to 1e-6 of its range."""
    with gzip.open(ROOT / "tests/data/prelibrary_launch_0p15s.csv.gz", "rt") as f:
        f.readline()
        names = f.readline().strip().split(",")
        ref = np.loadtxt(f, delimiter=",")
    res, _ = _launch(load(tune="default", overrides=OFF))
    assert res.shape[0] == ref.shape[0]
    for j, k in enumerate(names):
        a, b = res[:, CH.index(k)], ref[:, j]
        if k == "E_mech":
            # potential energy is defined up to a constant: the spring reference moved from full
            # open to belt contact (F_p0 g + k g^2/2 = 1.392496 J), so compare energy changes
            a, b = a - a[0], b - b[0]
        scale = max(np.abs(b).max(), 1e-12)
        assert np.abs(a - b).max() <= 1e-6 * scale, k


def test_default_tune_energy_channels_and_snapshots():
    """All 09 terms on: energy closes, new channels are consistent, snapshots follow the contract."""
    p = load(tune="default")
    res, snaps = _launch(p)
    c = lambda k: res[:, CH.index(k)]
    E, Win, Wd = c("E_mech"), c("W_in"), c("W_diss")
    assert np.max(np.abs((E - E[0]) - (Win - Wd))) < 1e-4 * Win[-1]
    assert np.all(np.isfinite(res))
    cb = np.cos(p["derived"]["beta"])
    assert np.allclose(c("cap_p"), 2 * p["belt"]["mu"] * c("r_eff_p") * c("Fa_p") / cb)
    assert np.allclose(c("w_wheel"), c("w_s") / p["vehicle"]["final_drive"])
    assert np.allclose(c("P_ground"), c("F_tire") * c("v"))
    assert np.all(c("T_belt") >= p["derived"]["belt_thermal"]["T_amb"] - 1e-9)
    assert c("v")[-1] > 0.0 and c("theta_e")[-1] > 0.0
    F, N = len(snaps["t"]), int(p["sim"]["belt_nodes"])
    assert F == 400
    assert snaps["belt"].shape == (F, N, 2) and snaps["tension"].shape == (F, N)
    for k in ("x_p", "x_s", "theta_e", "theta_s", "r_w"):
        assert snaps[k].shape == (F,)
    assert np.all(np.diff(snaps["t"]) > 0) and snaps["t"][0] == 0.0
    # frames sample the same trajectory as the series
    for k in ("x_p", "x_s", "theta_e", "theta_s", "r_w"):
        assert snaps[k][0] == c(k)[0]
    assert snaps["t"][-1] <= c("t")[-1]


def test_tune_library_sources_and_default():
    """Every component value has a source tag; the default tune equals the base clutch parameters."""
    from cvtsim.io.params import load_raw, COMPONENTS

    def leaves(node, path=""):
        if isinstance(node, dict) and "value" in node:
            yield path, node
        elif isinstance(node, dict):
            for k, v in node.items():
                if k not in ("description", "note"):
                    yield from leaves(v, f"{path}.{k}")

    for path, leaf in leaves(load_raw(COMPONENTS)):
        assert leaf.get("source"), path
    base, tuned = load(), load(tune="default")
    for sec in ("primary", "secondary"):
        for k, v in base[sec].items():
            assert tuned[sec][k] == v, (sec, k)
    for k in ("F_l0", "p_h", "theta_t0", "m_w"):
        assert tuned["derived"][k] == base["derived"][k]
    assert {"default", "heavy_flyweights", "stiff_primary", "progressive_helix"} <= set(tune_names())
    for name in tune_names():
        load(tune=name)
