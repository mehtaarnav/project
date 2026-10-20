"""Parts picker (cvtsim.app): form -> tune file -> same parameters as the library tune."""
import threading
import urllib.request
from http.server import ThreadingHTTPServer
from urllib.parse import urlencode

from cvtsim import app
from cvtsim.io.params import SLOTS, load

FORM = dict(app.DEFAULTS, shim_mm="0", hole_deg="30", t_end="8", label="pytest-default")


def _flat(d, pre=""):
    out = {}
    for k, v in d.items():
        if k == "tune":
            continue
        if isinstance(v, dict):
            out.update(_flat(v, pre + k + "."))
        else:
            out[pre + k] = v
    return out


def test_form_tune_reproduces_default(tmp_path, monkeypatch):
    monkeypatch.setattr(app, "LOCAL", tmp_path)
    path = app.write_tune(FORM)
    a, b = _flat(load(tune=str(path))), _flat(load(tune="default"))
    assert a.keys() == b.keys()
    obj = [k for k in a if hasattr(a[k], "__dict__")]  # profile objects: compare their samples
    diff = [k for k in a if k not in obj and str(a[k]) != str(b[k])]
    assert diff == []
    for k in obj:
        u = [0.0, 0.005, 0.01]
        assert [a[k](x) for x in u] == [b[k](x) for x in u] if callable(a[k]) else vars(a[k]).keys() == vars(b[k]).keys()


def test_page_lists_every_part():
    html = app.page()
    for slot, opts in app.parts().items():
        assert f'name="{slot}"' in html
        for name, _ in opts:
            assert f'value="{name}"' in html
    assert set(app.DEFAULTS) == set(SLOTS)


def test_post_launches_a_run(tmp_path, monkeypatch):
    monkeypatch.setattr(app, "LOCAL", tmp_path)
    seen = {}

    class FakePopen:
        def __init__(self, cmd, **kw):
            seen["cmd"] = cmd

        def poll(self):
            return None
    monkeypatch.setattr(app.subprocess, "Popen", FakePopen)
    srv = ThreadingHTTPServer(("127.0.0.1", 0), app.Handler)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    try:
        url = f"http://127.0.0.1:{srv.server_port}"
        r = urllib.request.urlopen(url + "/run", data=urlencode(FORM).encode())
        assert r.status == 200 and "pytest-default: running" in r.read().decode()
        assert seen["cmd"][4:6] == ["cvtsim.tune", "run"] and seen["cmd"][6].endswith("pytest-default.yaml")
        assert (tmp_path / "pytest-default.yaml").exists()
    finally:
        srv.shutdown()
