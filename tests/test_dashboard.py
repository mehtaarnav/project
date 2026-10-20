"""Engineering dashboard pages build from synthetic channels. The arrays are shaped like a launch
but are NOT simulator output: `python tests/test_dashboard.py` renders them to docs/img/ as layout
examples only."""
import pathlib
import sys

import numpy as np
import pytest
import matplotlib.pyplot as plt

from cvtsim.drivetrain.engine import EngineCurve
from cvtsim.io.params import load
from cvtsim.viz import dashboard as D, style as S

NEW = ["theta_e", "theta_s", "r_w", "F_fw", "F_spr_p", "F_spr_s", "F_tors", "F_cam", "cap_p", "cap_s",
       "P_slip_p", "P_slip_s", "T_belt", "v", "w_wheel", "F_tire", "N_rear", "P_ground", "a_x"]
OLD = ["t", "w_e", "w_s", "x_p", "x_s", "rstar_p", "rstar_s", "Fa_p", "Fa_s", "tq_p", "tq_s_fix",
       "tq_s_mov", "T_eng"]


def synthetic(scale=1.0, n=800, t_end=6.0):
    """Illustrative launch-like channels; scale stretches the shift (a second synthetic 'tune')."""
    p = load()
    t = np.linspace(0, t_end, n)
    s = 1 / (1 + np.exp(-(t - 1.2 * scale) / 0.25))          # shift progress 0..1
    rpm = 1500 + 1000 * np.clip(t / 0.5, 0, 1) + 100 * s + 0 * t
    rpm = np.where(t > 0.5, 2500 + 200 * np.tanh((t - 0.5) / 0.6) + 400 * s, rpm)
    w_e = rpm * np.pi / 30
    ratio = 3.0 - 2.1 * s
    w_s = w_e / ratio * np.clip((t - 0.4) / 0.2, 0, 1)
    v = w_s * p["vehicle"]["tire_radius"] / p["vehicle"]["final_drive"]
    T_eng = EngineCurve.from_params(p).torque(w_e)
    xs = p["derived"]["geometry"]["x_p_stop"]
    x_p = xs * s * (t > 0.45)
    Fa_p = np.where(t > 0.45, 900 + 1500 * s, 0.0)
    F_fw = Fa_p + 400 + 900 * s
    F_spr_p = 400 + 900 * s
    Fa_s = 1100 + 600 * (1 - s)
    tq_p = T_eng * 0.95 * (t > 0.45)
    cap_p = 2 * 0.9 * 0.05 * Fa_p / 0.97 * (1 + 0.3 * s) + 1e-3
    tq_s = tq_p * ratio
    a_x = np.gradient(v, t)
    P_slip = 40 * np.exp(-((t - 0.8 * scale) / 0.5) ** 2)
    d = dict(t=t, w_e=w_e, w_s=w_s, x_p=x_p, x_s=0.02 * (1 - s), rstar_p=0.03 + 0.04 * s,
             rstar_s=0.09 - 0.04 * s, Fa_p=Fa_p, Fa_s=Fa_s, tq_p=tq_p, tq_s_fix=tq_s * 0.5, tq_s_mov=tq_s * 0.5,
             T_eng=T_eng, F_fw=F_fw, F_spr_p=F_spr_p, F_spr_s=500 + 6e3 * (1 - s) * 0.1,
             F_tors=300 + 200 * s, F_cam=400 + 300 * (1 - s), cap_p=cap_p, cap_s=cap_p * ratio * 1.25,
             P_slip_p=P_slip, P_slip_s=0.6 * P_slip, T_belt=40 + 35 * (1 - np.exp(-t / 3)), v=v,
             F_tire=250 * a_x, N_rear=1470 + 250 * a_x * 0.5, P_ground=0.88 * T_eng * w_e * np.clip((t - 0.4) / 0.6, 0, 1),
             a_x=a_x)
    CH = list(d)
    return np.column_stack([d[k] for k in CH]), CH, p


def _n_axes(fig):
    return len(fig.axes)


def drop(res, CH, keys):
    keep = [i for i, k in enumerate(CH) if k not in keys]
    return res[:, keep], [CH[i] for i in keep]


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


PAGES = {"engine": (D.engine, 4), "clutch_forces": (D.clutch_forces, 4),
         "traction": (D.traction, 6), "performance": (D.performance, 12)}


@pytest.mark.parametrize("name", PAGES)
def test_page_builds(name):
    fn, n = PAGES[name]
    res, CH, p = synthetic()
    assert _n_axes(fn(res, CH, p)) == n


@pytest.mark.parametrize("name", PAGES)
def test_page_without_new_channels(name):
    fn, n = PAGES[name]
    res, CH, p = synthetic()
    res, CH = drop(res, CH, set(NEW))
    assert _n_axes(fn(res, CH, p)) == n


@pytest.mark.parametrize("name", PAGES)
def test_page_with_only_time_channels(name):
    fn, _ = PAGES[name]
    res, CH, p = synthetic()
    res, CH = drop(res, CH, set(CH) - {"t", "w_e", "w_s"})
    fn(res, CH, p)


def test_compare_axes_and_missing():
    a, CH, p = synthetic()
    b, _, _ = synthetic(scale=1.5)
    fig = D.compare([("base", a), ("slow shift", b)], CH, p)
    assert _n_axes(fig) == 4
    a2, CH2 = drop(a, CH, set(NEW) | {"rstar_p", "rstar_s"})
    b2, _ = drop(b, CH, set(NEW) | {"rstar_p", "rstar_s"})
    assert _n_axes(D.compare([("a", a2), ("b", b2)], CH2, p)) == 4
    D.compare([("only", a)], CH, p)


def test_compare_many_runs():
    a, CH, p = synthetic()
    D.compare([(f"r{i}", synthetic(1 + 0.1 * i)[0]) for i in range(6)], CH, p)


def test_key_metrics_against_hand_values():
    t = np.linspace(0, 10, 2001)
    v = 3.0 * t                                    # constant 3 m/s^2
    CH = ["t", "v", "w_e"]
    res = np.column_stack([t, v, np.full_like(t, 300.0)])
    m = D.key_metrics(res, CH, load())
    assert m["t_20"] == pytest.approx(20 / 3.6 / 3.0, abs=1e-6)
    assert m["t_40"] == pytest.approx(40 / 3.6 / 3.0, abs=1e-6)
    assert m["t_100ft"] == pytest.approx(np.sqrt(2 * 30.48 / 3.0), abs=1e-4)   # x = a t^2 / 2
    assert m["t_150ft"] == pytest.approx(np.sqrt(2 * 45.72 / 3.0), abs=1e-4)
    assert m["a_peak"] == pytest.approx(3.0, abs=1e-9)


def test_light_theme_builds_and_restores():
    res, CH, p = synthetic()
    S.use("light")
    try:
        assert S.theme == "light" and S.WHITE == S.THEMES["light"]["WHITE"]
        assert _n_axes(D.engine(res, CH, p)) == 4
        fig = D.performance(res, CH, p)
        assert fig.get_facecolor()[:3] != (0.0, 0.0, 0.0)
    finally:
        S.use("dark")
    assert S.RED == S.THEMES["dark"]["RED"]


def test_cmap_monotone_lightness():
    for th in ("dark", "light"):
        S.use(th)
        rgb = S.cmap()(np.linspace(0, 1, 64))[:, :3]
        lum = rgb @ [0.2126, 0.7152, 0.0722]
        d = np.diff(lum)
        assert (d >= -1e-6).all() or (d <= 1e-6).all()
    S.use("dark")


if __name__ == "__main__":
    out = pathlib.Path(__file__).resolve().parents[1] / "docs" / "img"
    out.mkdir(parents=True, exist_ok=True)
    res, CH, p = synthetic()
    res2, _, _ = synthetic(scale=1.5)
    res3, _, _ = synthetic(scale=0.7)
    for theme in ("dark", "light"):
        S.use(theme)
        figs = {"engine": D.engine(res, CH, p), "clutch_forces": D.clutch_forces(res, CH, p),
                "traction": D.traction(res, CH, p), "performance": D.performance(res, CH, p),
                "compare": D.compare([("baseline", res), ("slow shift", res2), ("fast shift", res3)], CH, p)}
        for k, f in figs.items():
            f.savefig(out / f"dashboard_{k}_{theme}.png")
            plt.close(f)
    print("wrote", out, file=sys.stderr)
