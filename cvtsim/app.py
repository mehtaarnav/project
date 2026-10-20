"""Local parts picker: choose every tunable part from dropdowns and launch a tracked run.

    python -m cvtsim.app            # then open http://localhost:8765

Each launch writes a tune file to local_tunes/ and runs `python -m cvtsim.tune run` in the
background. Finished runs are listed with their metrics, dashboard and mechanism viewer.
Standard library only; serves the repository root, so the viewer and run files load directly.
"""
from __future__ import annotations

import argparse
import html
import json
import subprocess
import sys
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs

import yaml

from cvtsim.io.params import COMPONENTS, SLOTS, load_raw

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_tunes"
RUNS = ROOT / "runs"
LABELS = {"primary_spring": "Primary spring", "flyweights": "Flyweights", "ramp": "Ramp",
          "secondary_spring": "Secondary spring", "torsion_spring": "Torsion spring",
          "helix": "Helix", "surface": "Surface"}
SHOWN = ("t_100ft_s", "engage_rpm", "t_shift_s", "v_end_kmh", "T_tight_max", "bearing_peak_N",
         "pto_side_peak_N", "vm_peak_MPa", "n_goodman_prim", "n_goodman_sec")
DEFAULTS = dict(primary_spring="brief", flyweights="brief", ramp="straight_45", secondary_spring="brief",
                torsion_spring="brief", helix="straight_45", surface="offroad")
JOBS: dict[str, subprocess.Popen] = {}


def parts():
    """{slot: [(name, description)]} from components.yaml."""
    comps = load_raw(COMPONENTS)
    return {slot: [(n, e.get("description", "")) for n, e in comps[sec].items()]
            for slot, (sec, _) in SLOTS.items()}


def torsion_holes(name):
    e = load_raw(COMPONENTS)["torsion_springs"][name]
    return e["torsion_holes_deg"]["value"]


def write_tune(form):
    """Tune YAML from the submitted form; returns its path."""
    comps = {slot: form[slot] for slot in SLOTS}
    label = form.get("label", "").strip() or "custom"
    label = "".join(ch if ch.isalnum() or ch in "-_" else "-" for ch in label)[:40]
    tune = dict(description=f"Built in cvtsim.app: {label}", components=comps,
                adjust=dict(primary_shim=dict(value=float(form["shim_mm"]) / 1e3, unit="m", source="app input"),
                            torsion_hole_deg=dict(value=float(form["hole_deg"]), source="app input")))
    LOCAL.mkdir(exist_ok=True)
    path = LOCAL / f"{label}.yaml"
    path.write_text(yaml.safe_dump(tune, sort_keys=False))
    return path


def launch(form):
    path = write_tune(form)
    cmd = [sys.executable, "-W", "ignore", "-m", "cvtsim.tune", "run", str(path), "--t-end", str(float(form["t_end"]))]
    log = open(LOCAL / f"{path.stem}.log", "w")
    JOBS[path.stem] = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    return path.stem


def recent_runs(n=15):
    out = []
    for d in sorted(RUNS.glob("*-tune-*"), reverse=True)[:n]:
        try:
            m = json.loads((d / "metrics.json").read_text())
        except (OSError, ValueError):
            continue
        out.append((d.name, m, (d / "viewer.json").exists(), (d / "dashboard.png").exists()))
    return out


def _select(name, options, chosen):
    opts = "".join(f'<option value="{html.escape(v)}" title="{html.escape(desc)}"'
                   f'{" selected" if v == chosen else ""}>{html.escape(v)}</option>' for v, desc in options)
    return f'<select name="{name}">{opts}</select>'


def page():
    ps = parts()
    rows = "".join(f"<label>{LABELS[s]}{_select(s, ps[s], DEFAULTS[s])}</label>" for s in SLOTS)
    holes = torsion_holes(DEFAULTS["torsion_spring"])
    hole = _select("hole_deg", [(f"{h:g}", "") for h in holes], f"{holes[0]:g}")
    jobs = "".join(f"<li>{html.escape(k)}: {'running' if p.poll() is None else 'exit ' + str(p.returncode)}</li>"
                   for k, p in JOBS.items())
    head = "".join(f"<th>{k}</th>" for k in SHOWN)
    body = ""
    for name, m, has_view, has_dash in recent_runs():
        cells = "".join(f"<td>{m[k]:.4g}</td>" if isinstance(m.get(k), (int, float)) else "<td></td>" for k in SHOWN)
        links = (f'<a href="/runs/{name}/dashboard.png">dashboard</a> ' if has_dash else "") + \
                (f'<a href="/viewer/index.html?run=../runs/{name}/viewer.json">viewer</a>' if has_view else "")
        body += f"<tr><td>{name}</td>{cells}<td>{links}</td></tr>"
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>cvtsim parts picker</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
body{{font:14px/1.4 system-ui,sans-serif;margin:24px;background:#f6f6f4;color:#1b1b1b}}
form{{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:12px;max-width:1000px}}
label{{display:flex;flex-direction:column;gap:4px;font-weight:600}}
select,input{{font:inherit;padding:6px;border:1px solid #bbb;border-radius:4px;background:#fff}}
button{{font:inherit;padding:8px 16px;border:0;border-radius:4px;background:#1f6feb;color:#fff;cursor:pointer}}
table{{border-collapse:collapse;margin-top:24px;font-size:13px}}td,th{{padding:4px 8px;border-bottom:1px solid #ddd;text-align:right}}
td:first-child,th:first-child{{text-align:left}}
</style></head><body>
<h1>cvtsim parts picker</h1>
<form method="post" action="/run">{rows}
<label>Torsion hole (deg){hole}</label>
<label>Primary shim (mm)<input name="shim_mm" type="number" step="0.1" value="0"></label>
<label>Launch length (s)<input name="t_end" type="number" step="0.5" min="0.5" value="8"></label>
<label>Run label<input name="label" value="custom"></label>
<div><button type="submit">Run launch</button></div></form>
<p>A launch of 8 s takes about 30 minutes. Hover an option for its description.</p>
<ul>{jobs}</ul>
<table><tr><th>run</th>{head}<th></th></tr>{body}</table>
</body></html>"""


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=str(ROOT), **k)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            return self._send(page())
        return super().do_GET()

    def do_POST(self):
        if self.path != "/run":
            return self.send_error(404)
        n = int(self.headers.get("Content-Length", 0))
        form = {k: v[0] for k, v in parse_qs(self.rfile.read(n).decode()).items()}
        try:
            launch(form)
        except (KeyError, ValueError) as e:
            return self._send(f"<p>Bad input: {html.escape(str(e))}</p><a href='/'>back</a>", 400)
        self.send_response(303)
        self.send_header("Location", "/")
        self.end_headers()

    def _send(self, text, code=200):
        data = text.encode()
        self.send_response(code)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, *a):
        pass


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.app")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    print(f"cvtsim parts picker on http://localhost:{a.port}")
    srv.serve_forever()


if __name__ == "__main__":
    main()
