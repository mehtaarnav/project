"""Uncertainty bands (docs/theory/11_engine_friction_uq.md, section 4).

    python -m cvtsim.uq screen [--t-end 2] [--workers N]
    python -m cvtsim.uq bands --k 5 --n 30 [--t-end 2] [--seed 0]

`screen` runs every input of io/uncertainty.yaml at its low and high value (others nominal) and
ranks the inputs by normalized span; `bands` runs a Latin-hypercube sample over the top k inputs
and reports the median and 5-95 % band of each metric. Runs execute in a process pool; each
default run is a tracked run (group "uq"), and the sweep is summarized in a uq run holding
tornado.png, bands.png and a markdown table.

A run_fn maps {parameter path: value} to a dict of scalar metrics. Inputs whose path is not
in the parameters yet (spec 11 physics not merged) are skipped; "build.kn_scale" is a build() kwarg.
"""
import argparse
import concurrent.futures as cf
import os
import pathlib
import uuid

import numpy as np
import yaml

from cvtsim import runs
from cvtsim.io.params import _num, load, load_raw

UNCERTAINTY = pathlib.Path(__file__).with_name("io") / "uncertainty.yaml"
BUILD_PREFIX = "build."
DECISION_METRICS = ("engage_rpm", "t_shift_s", "t_100ft_s", "v_end", "dist_m", "T_tight_max",
                    "T_belt_max_C", "bearing_peak_N", "vm_peak_MPa")
FT100 = 30.48  # m


def load_inputs(path=UNCERTAINTY):
    with open(path) as f:
        out = []
        for e in yaml.safe_load(f)["inputs"]:
            e = dict(e)
            for k in ("nominal", "low", "high"):
                e[k] = float(_num(e[k]))
            out.append(e)
    return out


def path_exists(p, path):
    if path.startswith(BUILD_PREFIX):
        return True
    d = p
    for k in path.split("."):
        if not isinstance(d, dict) or k not in d:
            return False
        d = d[k]
    return True


def available(inputs, p=None):
    """Split inputs into (usable, skipped) by whether their parameter path exists."""
    p = p if p is not None else load()
    ok = [e for e in inputs if path_exists(p, e["path"])]
    return ok, [e for e in inputs if e not in ok]


# ---------------------------------------------------------------- default run

def _t_100ft(t, dist):
    return float(np.interp(FT100, dist, t)) if dist[-1] >= FT100 else np.nan


def default_run(overrides, t_end=2.0, group="uq", record=True):
    """A short launch with the given parameter overrides; returns the decision metrics."""
    from cvtsim.cvt.system import CH, build, relax, simulate
    from cvtsim.fe.shaft_loads import shaft_histories
    from cvtsim.tune import metrics
    ov = dict(overrides)
    kn = ov.pop("build.kn_scale", 1.0)
    p = load(tune="default", overrides=ov)
    model = build(p, kn_scale=kn)
    relax(model)
    res, _ = simulate(model, t_end, nout=1500, n_frames=2)
    m = metrics(res, p)
    t, v = res[:, CH.index("t")], res[:, CH.index("v")]
    dist = np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(t))])
    m["t_100ft_s"] = _t_100ft(t, dist)
    sh = shaft_histories(res, p, CH)
    m.update(bearing_peak_N=max(sh[k].max() for k in sh["bearing_cols"]),
             vm_peak_MPa=max(sh["vm_prim"].max(), sh["vm_sec"].max()) / 1e6)
    if record:
        run = runs.Run(f"uq-{uuid.uuid4().hex[:8]}", dict(base=load_raw(), tune="default"),
                       overrides=dict(overrides, t_end=t_end), group=group)
        run.series(CH, res)
        run.log(**m)
        run.finish()
    return m


class DefaultRun:
    """Picklable default run_fn with a fixed t_end."""
    def __init__(self, t_end=2.0):
        self.t_end = t_end

    def __call__(self, overrides):
        return default_run(overrides, self.t_end)


def _map(run_fn, cases, workers):
    workers = workers if workers is not None else max(1, (os.cpu_count() or 2) - 1)
    if workers <= 1 or len(cases) <= 1:
        return [run_fn(c) for c in cases]
    with cf.ProcessPoolExecutor(max_workers=min(workers, len(cases))) as ex:
        return list(ex.map(run_fn, cases))


# ---------------------------------------------------------------- screen

def screen(run_fn, inputs, metrics=DECISION_METRICS, workers=None):
    """One-at-a-time low/high screen.

    Returns dict(base=metrics, rows=[dict(name, path, low, high, rel_low, rel_high, span)]) with
    rel_* = {metric: (value - nominal)/nominal}, span = max over metrics of |rel_high - rel_low|,
    rows sorted by span (largest first). Non-finite or zero-nominal metrics give NaN and do not
    count towards the span.
    """
    cases = [{}]
    for e in inputs:
        cases += [{e["path"]: e["low"]}, {e["path"]: e["high"]}]
    out = _map(run_fn, cases, workers)
    base = out[0]
    rows = []
    for i, e in enumerate(inputs):
        lo, hi = out[1 + 2 * i], out[2 + 2 * i]
        rel = lambda r: {k: ((r[k] - base[k]) / base[k]) if np.isfinite(r[k]) and np.isfinite(base[k]) and base[k] != 0
                         else np.nan for k in metrics}
        rl, rh = rel(lo), rel(hi)
        spans = [abs(rh[k] - rl[k]) for k in metrics if np.isfinite(rl[k]) and np.isfinite(rh[k])]
        rows.append(dict(name=e["name"], path=e["path"], low=e["low"], high=e["high"],
                         rel_low=rl, rel_high=rh, span=max(spans) if spans else 0.0))
    rows.sort(key=lambda r: -r["span"])
    return dict(base={k: base[k] for k in metrics}, rows=rows, metrics=list(metrics))


# ---------------------------------------------------------------- bands

def lhs(n, d, rng):
    """Latin hypercube on [0, 1]^d: each column has exactly one point per stratum [i/n, (i+1)/n)."""
    u = (np.arange(n)[:, None] + rng.random((n, d))) / n
    for j in range(d):
        u[:, j] = u[rng.permutation(n), j]
    return u


def bands(run_fn, inputs, k=5, n=30, seed=0, metrics=DECISION_METRICS, workers=None, scr=None, rng=None):
    """Latin-hypercube bands over the top-k inputs of a screen (computed here unless `scr` given).

    Returns dict(top=[input names], samples=(n, k) values, results=[metric dicts],
    stats={metric: (p5, median, p95)}, screen=scr).
    """
    scr = scr or screen(run_fn, inputs, metrics, workers)
    by_name = {e["name"]: e for e in inputs}
    top = [by_name[r["name"]] for r in scr["rows"][:k]]
    rng = rng or np.random.default_rng(seed)
    u = lhs(n, len(top), rng)
    lo, hi = np.array([e["low"] for e in top]), np.array([e["high"] for e in top])
    x = lo + u * (hi - lo)
    out = _map(run_fn, [{e["path"]: float(x[i, j]) for j, e in enumerate(top)} for i in range(n)], workers)
    stats = {}
    for m in metrics:
        v = np.array([r[m] for r in out], dtype=float)
        v = v[np.isfinite(v)]
        stats[m] = tuple(np.percentile(v, [5, 50, 95])) if v.size else (np.nan,) * 3
    return dict(top=[e["name"] for e in top], samples=x, unit=u, results=out, stats=stats, screen=scr)


# ---------------------------------------------------------------- reporting

def screen_table(scr):
    ms = scr["metrics"]
    lines = ["| input | low | high | span | " + " | ".join(ms) + " |", "|---|---|---|---|" + "---|" * len(ms)]
    lines.append("| nominal | | | | " + " | ".join(f"{scr['base'][m]:.4g}" for m in ms) + " |")
    for r in scr["rows"]:
        cells = [f"{r['rel_low'][m]:+.1%} / {r['rel_high'][m]:+.1%}" if np.isfinite(r["rel_low"][m]) and np.isfinite(r["rel_high"][m])
                 else "n/a" for m in ms]
        lines.append(f"| {r['name']} | {r['low']:.4g} | {r['high']:.4g} | {r['span']:.1%} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def bands_table(b):
    lines = ["| metric | p5 | median | p95 |", "|---|---|---|---|"]
    lines += [f"| {m} | {s[0]:.4g} | {s[1]:.4g} | {s[2]:.4g} |" for m, s in b["stats"].items()]
    return "\n".join(lines) + "\n"


def tornado_figure(scr, max_panels=6, top=12):
    import matplotlib.pyplot as plt
    from cvtsim.viz import style
    style.use("light")
    ms = [m for m in scr["metrics"] if np.isfinite(scr["base"][m])]
    ms = sorted(ms, key=lambda m: -max([abs(np.nan_to_num(r["rel_high"][m] - r["rel_low"][m])) for r in scr["rows"]] + [0]))[:max_panels]
    ncol = min(3, max(len(ms), 1))
    nrow = -(-max(len(ms), 1) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 0.34 * top * nrow + 1.6), squeeze=False)
    lo_c, hi_c = style.SERIES[0], style.SERIES[1]
    for ax, m in zip(axes.flat, ms):
        rows = sorted(scr["rows"], key=lambda r: -abs(np.nan_to_num(r["rel_high"][m] - r["rel_low"][m])))[:top][::-1]
        y = np.arange(len(rows))
        ax.barh(y, [np.nan_to_num(r["rel_low"][m]) * 100 for r in rows], color=lo_c, label="low")
        ax.barh(y, [np.nan_to_num(r["rel_high"][m]) * 100 for r in rows], color=hi_c, label="high")
        ax.set_yticks(y, [r["name"] for r in rows])
        ax.axvline(0, color=style.MUTED, lw=0.8)
        ax.set_xlabel("change vs nominal (%)")
        style.title(ax, f"{m} (nominal {scr['base'][m]:.4g})")
    axes.flat[0].legend()
    for ax in list(axes.flat)[len(ms):]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def bands_figure(b):
    import matplotlib.pyplot as plt
    from cvtsim.viz import style
    style.use("light")
    ms = [m for m in b["stats"] if np.isfinite(b["stats"][m][1])]
    ncol = min(3, max(len(ms), 1))
    nrow = -(-max(len(ms), 1) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 2.6 * nrow), squeeze=False)
    for ax, m in zip(axes.flat, ms):
        v = np.array([r[m] for r in b["results"]], dtype=float)
        p5, p50, p95 = b["stats"][m]
        ax.hist(v[np.isfinite(v)], bins=min(12, max(3, len(v) // 3)), color=style.SERIES[0])
        ax.axvspan(p5, p95, color=style.SERIES[1], alpha=0.2, label="5-95 %")
        ax.axvline(p50, color=style.SERIES[1], lw=1.5, label="median")
        style.title(ax, m)
    axes.flat[0].legend()
    for ax in list(axes.flat)[len(ms):]:
        ax.axis("off")
    fig.suptitle("Latin-hypercube bands over inputs: " + ", ".join(b["top"]), fontsize=9)
    fig.tight_layout()
    return fig


def summarize(name, scr, b=None, extra=None):
    """Write the uq summary run: markdown table plus tornado.png (and bands.png)."""
    import matplotlib.pyplot as plt
    run = runs.Run(name, dict(inputs=load_inputs(), **(extra or {})), group="uq")
    md = "# Uncertainty screen\n\nEntries are (low - nominal)/nominal / (high - nominal)/nominal.\n\n" + screen_table(scr)
    fig = tornado_figure(scr)
    run.figure(fig, "tornado")
    plt.close(fig)
    if b is not None:
        md += f"\n# Bands (top {len(b['top'])}: {', '.join(b['top'])}; n = {len(b['results'])})\n\n" + bands_table(b)
        fig = bands_figure(b)
        run.figure(fig, "bands")
        plt.close(fig)
        run.log(**{f"{m}_{q}": v for m, s in b["stats"].items() for q, v in zip(("p5", "p50", "p95"), s)})
    (run.dir / "uq.md").write_text(md)
    run.finish()
    runs.write_index()
    return run.dir, md


def main(argv=None, run_fn=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.uq")
    ap.add_argument("cmd", choices=("screen", "bands"))
    ap.add_argument("--t-end", type=float, default=2.0)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--only", nargs="+", default=None, help="screen only these input names")
    a = ap.parse_args(argv)
    inputs, skipped = available(load_inputs())
    if a.only:
        unknown = set(a.only) - {e["name"] for e in inputs}
        if unknown:
            ap.error(f"unknown input(s): {', '.join(sorted(unknown))}")
        inputs = [e for e in inputs if e["name"] in a.only]
    if skipped:
        print("skipped (path not in parameters yet):", ", ".join(e["path"] for e in skipped))
    run_fn = run_fn or DefaultRun(a.t_end)
    scr = screen(run_fn, inputs, workers=a.workers)
    b = bands(run_fn, inputs, a.k, a.n, a.seed, workers=a.workers, scr=scr) if a.cmd == "bands" else None
    d, md = summarize(f"uq-{a.cmd}", scr, b, dict(t_end=a.t_end, k=a.k, n=a.n, seed=a.seed, only=a.only))
    print(d)
    print(md)


if __name__ == "__main__":
    main()
