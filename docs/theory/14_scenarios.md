# 14 — Event scenarios (`cvtsim.cvt.system.schedule`, `cvtsim/io/scenarios/`)

A launch alone does not produce the largest loads. Hills, backshifts and landings do. A scenario is
a table of knots `[t, throttle, grade, traction]`, interpolated linearly in time. Values hold
constant before the first knot and after the last.

## 1. Model terms
* **Throttle** `θ ∈ [0, 1]` scales the governed map torque. Below full throttle the engine also
  absorbs a motoring torque (pumping and friction):
      T_e = θ T_map(ω) − (1 − θ) T_mot tanh(ω / 10 rad/s)
  `T_mot` is `engine.motoring_torque` (ASSUMED, range in the parameter file). The tanh term only
  removes the sign discontinuity at zero speed. Firing pulses keep the cycle mean `T_e`
  (11 section 1).
* **Grade** `φ` (uphill positive) adds `m g sin φ` to the road load. The tire normal load stays the
  level-ground value (assumption A14.1).
* **Traction factor** `λ ≥ 0` multiplies the rear-tire normal load and the tire force:
  `N_r → λ N_r`, `F_x → λ F_x`.
  * `λ = 0`: airborne. The wheels spin up freely under the driveline.
  * `λ > 1`: landing load. The wheels are forced back to ground speed, and the driveline carries
    the impact.

At `θ = 1, φ = 0, λ = 1` every term reduces exactly to the launch model of 09.

## 2. Supplied scenarios
| Name | Event | Knots |
|---|---|---|
| `launch` | full-throttle launch, level | constant |
| `hill_climb` | level launch, then a 20° climb from 2 s | grade 0 → 20° over 2–3 s |
| `backshift` | throttle lift at 4 s for 1.5 s, then full throttle | θ 1 → 0 → 1 |
| `landing` | airborne 3.0–3.6 s at full throttle, then a 3 g landing | λ 1 → 0 → 3 → 1 |

The grade, airtime and landing load are reporting choices. Replace them with course or logged data.

## 3. Targets
| ID | Case | Target |
|---|---|---|
| T14.1 | `launch` scenario vs no scenario | identical results (bit for bit) |
| T14.2 | airborne (`λ = 0`) on a 10° grade, from rest | `a_x = −(Crr m g tanh(v/V_s) + c_a v|v| + m g sin 10°)/m` at every output sample, rel 1e-9 |
| T14.3 | closed throttle | output `T_e = −T_mot tanh(ω_e/10)` at every output sample, rel 1e-9 |
| T14.4 | throttle step, grade and traction change in one run | energy closure < 1e-4 of input work |

## 4. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A14.1 | Tire normal load not reduced by `cos φ` on a grade | ↑ traction overestimated by up to 6 % at 20° |
| A14.2 | Motoring torque proportional to `1 − θ`, constant in speed | ? (no data) |
| A14.3 | A landing is a traction-factor pulse. No suspension dynamics; the vertical load history is prescribed | ? (replace with a logged load) |
