"""Phase-0 targets script must run and reproduce the headline numbers quoted in the specs."""
import subprocess, sys, pathlib

def test_targets_script():
    root = pathlib.Path(__file__).resolve().parents[1]
    out = subprocess.run([sys.executable, str(root / "docs/theory/targets.py")],
                         capture_output=True, text=True, check=True).stdout
    vals = {l.split()[0]: float(l.split()[1]) for l in out.splitlines() if l.strip()}
    assert abs(vals["engage_rpm_kappa0.5"] - 1725.66) < 0.01
    assert abs(vals["negstiff_rpm_kappa0.5"] - 2262.57) < 0.01
    assert abs(vals["capstan_ratio_pi"] - 3.79189) < 1e-5
    assert abs(vals["geo_x_p_stop"] - 0.0178192) < 1e-6
