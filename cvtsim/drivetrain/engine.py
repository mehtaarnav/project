"""Engine full-throttle torque curve (docs/theory/06_drivetrain.md, section 2)."""
import numpy as np

RPM_TO_RAD = 2.0 * np.pi / 60.0


class EngineCurve:
    """Piecewise-linear full-throttle torque vs rpm.

    Below the first chart point the torque is held constant (unverified region).
    Above the last chart point the torque ramps linearly to zero at rpm_cut.
    """

    def __init__(self, rpm, torque, rpm_cut=3800.0):
        self.rpm = np.asarray(rpm, dtype=float)
        self.torque_tab = np.asarray(torque, dtype=float)
        self.rpm_cut = float(rpm_cut)
        self.rpm_min = float(self.rpm[0])
        self.rpm_max = float(self.rpm[-1])

    @classmethod
    def from_params(cls, p):
        e = p["engine"]
        return cls(e["rpm"], e["torque"], e["rpm_cut"])

    def torque(self, omega):
        """Engine torque [N m] at crank speed omega [rad/s]. Vectorized."""
        r = np.asarray(omega, dtype=float) / RPM_TO_RAD
        base = np.interp(r, self.rpm, self.torque_tab)
        taper = np.clip((self.rpm_cut - r) / (self.rpm_cut - self.rpm_max), 0.0, 1.0)
        out = base * taper
        return float(out) if out.ndim == 0 else out

    def flag_unverified(self, omega):
        """True where the engine speed is below the chart (torque not verified)."""
        r = np.asarray(omega, dtype=float) / RPM_TO_RAD
        res = r < self.rpm_min
        return bool(res) if res.ndim == 0 else res
