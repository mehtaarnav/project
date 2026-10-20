"""Uncertainty framework (docs/theory/11, section 4): T11.7 and the machinery around it."""
import numpy as np
import pytest

from cvtsim import runs, uq
from cvtsim.io.params import load

INPUTS = uq.load_inputs()
NOM = {e["path"]: e["nominal"] for e in INPUTS}
COEF = {e["path"]: (i + 1) * (-1) ** i for i, e in enumerate(INPUTS)}  # synthetic sensitivities a_i
C0 = 1000.0


def linear_model(overrides):
    """y1 = C0 + sum a_i x_i ; y2 = 2 y1 (a second metric); module level so a process pool can pickle it."""
    y = C0 + sum(a * overrides.get(k, NOM[k]) for k, a in COEF.items())
    return dict(m1=y, m2=2.0 * y, m_nan=np.nan, **{k: y for k in uq.DECISION_METRICS})


def test_t11_7_tornado_exact():
    """Linear model: tornado entries equal a_i (x - x_nom) / y_nom exactly; sorted by span."""
    scr = uq.screen(linear_model, INPUTS, metrics=("m1", "m2", "m_nan"), workers=1)
    y0 = C0 + sum(a * NOM[k] for k, a in COEF.items())
    assert scr["base"]["m1"] == pytest.approx(y0, rel=1e-14)
    spans = [r["span"] for r in scr["rows"]]
    assert spans == sorted(spans, reverse=True)
    for r in scr["rows"]:
        a = COEF[r["path"]]
        assert r["rel_low"]["m1"] == pytest.approx(a * (r["low"] - NOM[r["path"]]) / y0, rel=1e-12, abs=1e-15)
        assert r["rel_high"]["m1"] == pytest.approx(a * (r["high"] - NOM[r["path"]]) / y0, rel=1e-12, abs=1e-15)
        assert r["rel_high"]["m2"] == pytest.approx(r["rel_high"]["m1"], rel=1e-12, abs=1e-15)
        assert np.isnan(r["rel_low"]["m_nan"])
        assert r["span"] == pytest.approx(abs(a) * (r["high"] - r["low"]) / y0, rel=1e-12)


def test_parallel_matches_serial():
    a = uq.screen(linear_model, INPUTS[:4], metrics=("m1",), workers=1)
    b = uq.screen(linear_model, INPUTS[:4], metrics=("m1",), workers=2)
    assert [r["name"] for r in a["rows"]] == [r["name"] for r in b["rows"]]
    assert [r["span"] for r in a["rows"]] == [r["span"] for r in b["rows"]]


def test_lhs_stratification():
    rng = np.random.default_rng(3)
    n, d = 30, 5
    u = uq.lhs(n, d, rng)
    assert u.shape == (n, d) and (u >= 0).all() and (u < 1).all()
    for j in range(d):
        assert sorted(np.floor(u[:, j] * n).astype(int)) == list(range(n))


def test_bands_use_top_k_and_ranges():
    b = uq.bands(linear_model, INPUTS, k=3, n=20, seed=1, metrics=("m1",), workers=1)
    spans = {r["name"]: r["span"] for r in b["screen"]["rows"]}
    top3 = sorted(spans, key=lambda n: -spans[n])[:3]
    assert b["top"] == top3
    by = {e["name"]: e for e in INPUTS}
    for j, name in enumerate(b["top"]):
        assert (b["samples"][:, j] >= by[name]["low"]).all() and (b["samples"][:, j] <= by[name]["high"]).all()
    p5, p50, p95 = b["stats"]["m1"]
    assert p5 <= p50 <= p95


def test_missing_paths_skipped_and_nominals_match_parameters():
    p = load()
    ok, skipped = uq.available(INPUTS, p)
    assert {"belt.mu", "build.kn_scale", "vehicle.CdA", "sheave.stop_k"} <= {e["path"] for e in ok}
    assert all(e["path"] not in {s["path"] for s in skipped} for e in ok)
    assert all("requires spec-11 physics" in e["source"] for e in skipped)
    assert "clutch.mu_helix" in {e["path"] for e in skipped} or uq.path_exists(p, "clutch.mu_helix")
    for e in ok:
        if e["path"].startswith("build."):
            continue
        d = p
        for k in e["path"].split("."):
            d = d[k]
        assert d == pytest.approx(e["nominal"], rel=1e-12), e["path"]
        assert e["low"] < e["nominal"] < e["high"]
        assert e["source"]


def test_cli_smoke(tmp_path, monkeypatch):
    monkeypatch.setattr(runs, "ROOT", tmp_path)
    uq.main(["bands", "--k", "3", "--n", "6", "--workers", "1"], run_fn=linear_model)
    out = list(tmp_path.glob("*-uq-bands"))
    assert len(out) == 1
    for f in ("tornado.png", "bands.png", "uq.md", "metrics.json", "config.yaml"):
        assert (out[0] / f).exists()
    assert "| belt_mu |" in (out[0] / "uq.md").read_text()
