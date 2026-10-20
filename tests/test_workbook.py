"""Team workbook reader: a workbook in the CVT_Inputs template layout reproduces the brief."""
import pytest

openpyxl = pytest.importorskip("openpyxl")
from cvtsim.io.workbook import configuration  # noqa: E402


def _template(path):
    wb = openpyxl.Workbook()
    sheets = {
        "Configurations": [["Name", "Primary Spring", "Secondary Spring", "Secondary Pretension",
                            "Flyweights", "Secondary Cam", "Exclude", "Notes"],
                           ["Boombox", "D", "A", 30, 275, "A", None, None]],
        "Primary Spring": [["Spring", "SPRING RATE", "Uncompressed Length", "Minimum Preload "],
                           ["Units", "lbf/in", "in", "in"], ["D", 75.42, 4.254, 1.683]],
        "Secondary Spring": [["Spring", "LINEAR SPRING RATE", "TORSIONAL SPRING RATE",
                              "Uncompressed Length", "Minimum Linear Preload"],
                             ["Units", "lbf/in", "ft-lbf/deg", "in", "in"], ["A", 19.15, 0.05059, 6.98, 4.48]],
        "Flyweights": [["Label", "Weight"], ["Units", "g"], [275, 275]],
        "Secondary Cam": [["Cam", "Ramp Aangle", "Ramp Radius"], ["Units", "deg", "in"], ["A", 45, 1.5]],
        "Constants": [["Constants", "Variable", "Values", "Units"], ["", "mC", 38.7, "g"],
                      ["", "nFly", 3, "dim"], ["", "alpha", 45, "deg"], ["", "thetaP", 45, "deg"],
                      ["", "rFly", 1.75, "in"], ["", "GR", 8.1, "dim"], ["", "dTire", 23, "in"],
                      ["", "L", 0.895, "m"], ["", "muT", 0.3, "dim"], ["", "theta", 12, "deg"],
                      ["", "xc", 0.232, "m"], ["", "xPi", 0.353, "in"], ["", "rPmax", 3.6, "in"],
                      ["", "rPmin", 1.04, "in"], ["", "rSmax", 3.82, "in"], ["", "rSmin", 1.67, "in"]],
        "Dyno Data": [["Engine Speed", "Torque"], ["rpm", "ft-lbf"], [2400, 18.5], [3600, 13.5]],
    }
    wb.remove(wb.active)
    for name, rows in sheets.items():
        ws = wb.create_sheet(name)
        for r in rows:
            ws.append(r)
    wb.save(path)


def test_boombox_reproduces_brief(tmp_path):
    path = tmp_path / "inputs.xlsx"
    _template(path)
    o = configuration(path, "Boombox")
    assert o["primary.k_spring"] == pytest.approx(13208.0, rel=1e-4)
    assert o["primary.F_preload"] == pytest.approx(683.04, rel=1e-4)
    assert o["secondary.k_linear"] == pytest.approx(3353.66, rel=1e-4)
    assert o["secondary.k_torsion"] == pytest.approx(3.930, rel=1e-3)
    assert o["secondary.F_preload_equiv"] == pytest.approx(435.63, rel=1e-4)
    assert o["primary.kappa"] == pytest.approx(0.5, rel=1e-12)
    assert o["engine.torque"][0] == pytest.approx(25.08, rel=1e-3)
