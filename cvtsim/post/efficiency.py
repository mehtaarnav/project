"""Power accounting and CVT efficiency map (docs/theory/13 sec 2). SI; magnitudes used so the
efficiency is defined for either speed sign. `res` is (n_out, len(CH)) from simulate()."""
import numpy as np


def _col(res, CH, name):
    return res[:, CH.index(name)]


def power_terms(res, CH):
    """Per-instant power terms [W]: P_eng, P_belt_in, P_out, P_slip, and the residual
    R_bal = P_belt_in - P_out - P_slip (inertia and belt hysteresis; diagnostic)."""
    w_e = _col(res, CH, "w_e")
    w_s = _col(res, CH, "w_s")
    P_eng = _col(res, CH, "T_eng") * w_e
    P_in = np.abs(_col(res, CH, "tq_p")) * np.abs(w_e)
    P_out = np.abs(_col(res, CH, "tq_s_fix") + _col(res, CH, "tq_s_mov")) * np.abs(w_s)
    P_slip = _col(res, CH, "P_slip_p") + _col(res, CH, "P_slip_s")
    return dict(P_eng=P_eng, P_belt_in=P_in, P_out=P_out, P_slip=P_slip,
                R_bal=P_in - P_out - P_slip)


def cvt_efficiency(res, CH, p):
    """eta_cvt = P_out / P_belt_in where P_belt_in > p['P_min'], NaN elsewhere."""
    t = power_terms(res, CH)
    P_in = t["P_belt_in"]
    eta = np.full(P_in.shape, np.nan)
    ok = P_in > p["P_min"]
    eta[ok] = t["P_out"][ok] / P_in[ok]
    return eta


def quasi_steady_mask(res, CH, p):
    """True where eta is defined and P_slip <= P_belt_in. Instants with slip power above the
    input (engagement transients: slip is paid from stored belt and inertia energy) are
    excluded from the map and from the [0, 1] range check (docs/theory/13 sec 2)."""
    t = power_terms(res, CH)
    return np.isfinite(cvt_efficiency(res, CH, p)) & (t["P_slip"] <= t["P_belt_in"])


def efficiency_map(res, CH, p):
    """Mean eta_cvt over quasi-steady instants in bins of g = rstar_s/rstar_p and primary rpm
    n_p = 30 w_e / pi.

    Returns (mean, counts, g_edges, rpm_edges); mean is NaN in empty bins, shape
    (len(g_edges)-1, len(rpm_edges)-1).
    """
    eta = cvt_efficiency(res, CH, p)
    ok = quasi_steady_mask(res, CH, p)
    g = _col(res, CH, "rstar_s") / _col(res, CH, "rstar_p")
    rpm = _col(res, CH, "w_e") * 30.0 / np.pi
    g_e = np.asarray(p["ratio_edges"], dtype=float)
    r_e = np.asarray(p["rpm_edges"], dtype=float)
    sums, _, _ = np.histogram2d(g[ok], rpm[ok], bins=[g_e, r_e], weights=eta[ok])
    counts, _, _ = np.histogram2d(g[ok], rpm[ok], bins=[g_e, r_e])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(counts > 0, sums / counts, np.nan)
    return mean, counts, g_e, r_e
