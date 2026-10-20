# 06 — Drivetrain: engine, vehicle, road load (`cvtsim.drivetrain`)

## 1. Sources
* Engine torque: brief, Kohler CH440 team-restricted chart (2400–3600 rpm).
  No public source checked, so **UNVERIFIED**.
* Road load (rolling resistance `C_rr m g`) and reflected inertia: T. D.
  Gillespie, *Fundamentals of Vehicle Dynamics*, SAE, 1992, Ch. 2 & 4
  (chapter numbers **UNVERIFIED**). The formulas are derived below.

## 2. Engine
Full-throttle torque `T_e(ω)`: piecewise-linear interpolation of the chart

| rpm | 2400 | 2600 | 2800 | 3000 | 3200 | 3400 | 3600 |
|---|---|---|---|---|---|---|---|
| N·m | 25.08 | 24.54 | 23.59 | 22.51 | 20.88 | 19.66 | 18.30 |

* Below 2400 rpm: **no data**. Held constant at 25.08 N·m (A6.1, flagged
  `UNVERIFIED_BELOW_2400`).
* Above 3600 rpm (governed engine): torque ramps linearly to 0 at 3800 rpm
  (A6.2, assumption).
* Engine friction is not separated from the chart. The chart is treated as net torque.
* Rotating inertia `J_e`: PLACEHOLDER 0.05 kg·m² (crank, flywheel, fixed primary
  half). Not supplied. Sensitivity is reported.

## 3. Vehicle reflected to the secondary shaft
Final drive ratio `G = 8.1` (secondary : wheel speed), tire radius `r_t`.
Rigid driveline: `v = ω_s r_t / G`. Kinetic energy `½ m v² = ½ (m r_t²/G²) ω_s²`, so

    J_veh = m r_t² / G² = 250·0.2921²/8.1² = 0.325112 kg·m²               (D1)

Wheel/axle rotational inertia is not supplied and is neglected (A6.3).
Road-load torque at the secondary:

    T_load = (C_rr m g · tanh(v/v0) + ½ ρ_air C_dA v²) r_t / G             (D2)

with `v0 = 0.05 m/s` (smooth sign) and `C_dA = 0` by default. Aerodynamic data
was not supplied (A6.4). Driveline efficiency is 1 (A6.5).

## 4. Steady state (fixed ratio test)
At fixed CVT ratio `R = ω_e/ω_s`, the speed balance `T_e(ω_e) R = T_load` with
`ω_e = R G v / r_t`. For the linear test engine `T(ω) = a − bω` and
`F = C_rr m g + c_a v²` (`c_a = ½ρC_dA`), multiply by `G/r_t`:

    c_a v² + b k² v + (C_rr m g − a k) = 0,   k = R G / r_t               (D3)

## 5. Symbols
| Symbol | Meaning | Value |
|---|---|---|
| `m` | vehicle + driver mass | 250 kg (brief) |
| `G` | final drive | 8.1 (brief) |
| `r_t` | tire radius | 0.2921 m (brief) |
| `C_rr` | rolling resistance | 0.015 (brief) |
| `J_e` | engine rotating inertia | PLACEHOLDER 0.05 kg·m² |
| `J_sec` | secondary shaft + fixed half | PLACEHOLDER 2.0e-3 kg·m² |

## 6. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A6.1 | T_e held at 25.08 N·m below 2400 rpm | ? — likely ↑ launch torque (CH440 torque typically falls at low rpm) |
| A6.2 | Governor: linear cut 3600→3800 rpm | ? |
| A6.3 | Wheel/axle/gear inertia neglected | ↑ acceleration, ↓ inertial shaft torques downstream |
| A6.4 | No aerodynamic drag | ↑ top speed; negligible during launch |
| A6.5 | Lossless gearbox/chain | ↑ wheel torque ~3–5 % |
| A6.6 | Rigid driveline (no torsional compliance, no tire slip) | ↑ high-frequency torque peaks transmitted; ↓ oscillatory amplification |

## 7. Implementation contract (`cvtsim/drivetrain/`)
* `engine.py`: `EngineCurve(rpm, torque, rpm_cut=3800)` with `.torque(omega)` (vectorized) and `.flag_unverified(omega) -> bool` (ω < 2400 rpm).
* `vehicle.py`: `reflected_inertia(m, r_t, G)`, `road_load_torque(omega_s, p)`, `vehicle_speed(omega_s, p)`.
* `steady.py`: `steady_speed_linear(a, b, R, G, r_t, Crr, m, c_a)` (D3), and `simulate_fixed_ratio(engine_torque_fn, J_e, R, p, t_end, h)`.
  The ODE is `(J_e R² + J_veh + J_sec) ω̇_s = T_e(R ω_s) R − T_load(ω_s)`.

## 8. NUMERIC TEST TARGETS (`tests/test_drivetrain.py`)
| Test | Target | Tol |
|---|---|---|
| T6.1 `EngineCurve.torque` at 2400, 3000, 3600 rpm | 25.08, 22.51, 18.30 | 1e-12 |
| T6.2 at 2900 rpm | (23.59+22.51)/2 = 23.05 | 1e-12 |
| T6.3 at 2000 rpm: 25.08 and `flag_unverified` True | – | – |
| T6.4 at 3700 rpm: 9.15; at 3900 rpm: 0 | – | 1e-12 |
| T6.5 `reflected_inertia(250, 0.2921, 8.1)` | 0.325112 kg·m² | rel 1e-6 |
| T6.6 (D3) with a=40, b=0.05, R=1, G=8.1, r_t=0.2921, C_rr=0.015, m=250, c_a=0.5 | v = 21.744 m/s (engine 5757.91 rpm) | rel 1e-4 |
| T6.7 `simulate_fixed_ratio` with the same linear engine, J_e=0.05, 60 s | final speed → 21.744 m/s | rel 1e-3 |
| T6.8 rolling only: road-load torque at v=10 m/s | 0.015·250·9.81·0.2921/8.1 = 1.32662 N·m (tanh≈1) | rel 1e-6 |
