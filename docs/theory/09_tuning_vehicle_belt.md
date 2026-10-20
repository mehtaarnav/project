# 09 — Tunable components, realistic vehicle, belt heating, outputs

Extends 05 (clutches), 06 (drivetrain) and 03 (belt). Notation and frames as in README.

## 1. Primary with a profiled ramp

Travel from belt contact `u = x_p ∈ [−g, x_p,stop]` (the team CVT_model_v4 convention: the spring preload and flyweight radius are given at belt contact). `g` is the travel from full open to contact; `u < 0` extrapolates the ramp linearly. Flyweight radius `r_w(u)` is a
tabulated, strictly increasing profile, interpolated with a monotone cubic (PCHIP). Then
`κ(u) = dr_w/du` and `κ'(u) = d²r_w/du²`. A straight ramp is `r_w = r_w0 + κ0 u`
(κ0 = 0.5 for the double 45° face, 05 §3).

    𝒯 = ½ (J_e + n m_w r_w²) ω² + ½ (m_ph + n m_w κ²) u̇²

Lagrange (derivation: `d/dt(M u̇) = M ü + 2 n m_w κ κ' u̇²`, and `∂𝒯/∂u = n m_w r_w κ ω² + n m_w κ κ' u̇²`):

    (J_e + n m_w r_w²) ω̇ + 2 n m_w r_w κ u̇ ω = T_eng + τ_belt − τ_gov                  (P1')
    (m_ph + n m_w κ²) ü + n m_w κ κ' u̇² − n m_w r_w κ ω² = −F_spr(u) − c_p u̇ − F_a,p + F_stop   (P2')

Spring `F_spr(u) = F_p0 + k_p u`. The preload `F_p0` is a library value plus an optional
shim offset `Δ`, which adds `k_p Δ`. Static flyweight force: `F_fw(u, ω) = n m_w r_w(u) κ(u) ω²`.
Engagement: `F_fw(0, ω) = F_p0`.

## 2. Secondary with a profiled helix

Opening travel `x ∈ [0, s_s]`. The helix is specified by its angle `α(x)` from the shaft
axis at cam radius `r_c`. A multi-angle helix is piecewise constant. Relative twist of the
movable half: `h(x) = ∫₀ˣ dx'/(r_c tan α(x'))`, so `h' = 1/(r_c tan α)`. A straight 45°
helix reproduces `p_h = r_c tan α` of 05 §4. With `Ω_m = ω − h' ẋ`:

    𝒯 = ½ J_s ω² + ½ J_m Ω_m² + ½ m_m ẋ²,   U = F_l0 x + ½ k_l x² + ½ k_t (θ_t0 + h(x))²

Lagrange gives (the `h''` terms come from `Ω̇_m = ω̇ − h'ẍ − h''ẋ²`; derivation:
`∂𝒯/∂ẋ = −J_m h' Ω_m + m_m ẋ`, whose time derivative carries `−J_m h'' ẋ Ω_m`, cancelled by
`∂𝒯/∂x = −J_m h'' ẋ Ω_m`, leaving `+J_m h' h'' ẋ²`; energy rate is `ω Q_θ + ẋ Q_x` exactly):

    [ J_s+J_m     −J_m h'       ] [ω̇]   [ Q_θ + J_m h'' ẋ²        ]
    [ −J_m h'   m_m + J_m h'²   ] [ẍ] = [ Q_x − J_m h' h'' ẋ²     ]

    Q_θ = τ_fix + τ_mov − T_load
    Q_x = −h' τ_mov + F_a,s − (F_l0 + k_l x) − k_t (θ_t0 + h) h' − c_s ẋ + F_stop

**Angle breaks.** A sharp break makes `h''` a Dirac (an impulse the explicit integrator
cannot carry). In code (`cvtsim/cvt/profiles.py`) `h'` changes over a length `w` centred on
the break with the smootherstep `10t³ − 15t⁴ + 6t⁵`, so `h''` and `h'''` stay continuous.
The blend is symmetric, so `h(x)` outside `[x_b − w/2, x_b + w/2]` equals the sharp-break
integral exactly. `w` = 2 mm is ASSUMED (a real helix has a finite transition too).
Ramp tables use PCHIP, which is C1: `κ'` jumps at the knots, and an integrator crossing a
knot makes an O(h) energy error per crossing (negligible at the launch step 2.5e-7 s).

Friction on the movable flank uses the speed `Ω_m`. Static closing force:
`F_s(x, T) = F_l0 + k_l x + [T/2 + k_t(θ_t0 + h(x))] h'(x)`. Torsion pretwist `θ_t0`
is selected from discrete hole positions. The 435.63 N "equivalent preload" convention
applies only to the brief's default tune.

## 3. Vehicle

The rear wheels are rigidly geared to the secondary, `ω_w = ω_s/G`. The vehicle speed `v` is a
separate state, and the tire couples the two by slip. The secondary inertia `J_s` therefore
contains `J_shaft + J_w/G²` and no longer the reflected vehicle mass. The load torque on the
secondary replaces (D2) of 06:

    T_ld = F_x r_t / (G η)   if F_x ω_s ≥ 0 (driving),   T_ld = η F_x r_t / G   otherwise
    m v̇  = F_x − F_roll − F_aero
    F_x  = μ_t N_r tanh(s/s0),   s = (ω_w r_t − v)/max(|v|, v_min)
    N_r  = m g (b/L_wb) + m v̇ h_cg / L_wb       (v̇ lagged one step)
    F_roll = C_rr m g tanh(v/0.05),   F_aero = ½ ρ C_dA v |v|

Dissipated power: driveline loss `T_ld ω_s − F_x ω_w r_t ≥ 0`, tire slip
`F_x (ω_w r_t − v) ≥ 0`, road load `(F_roll + F_aero) v`. Kinetic energy adds `½ m v²`.

Governor: `τ_gov = T_chart · clip((rpm − (rpm_cut − Δ))/Δ, 0, 1)` with `rpm_cut` = 3800 and
`Δ` = 100 rpm; the chart is held at 18.30 N·m from 3600 rpm up to the cut (UNVERIFIED).
`Δ` = 200 reproduces the 06 taper 3600→3800 exactly.

`vehicle.model: rigid` keeps the 06 driveline (`v = ω_s r_t/G`, vehicle mass reflected into
`J_s`, load torque from the road load through `η`). "New vehicle terms off" (T9.9) means
`model: rigid, η = 1, J_w = 0, C_dA = 0, Δ = 200 rpm`. The brief does not supply `C_dA`, so the
base parameter file keeps 0 and the tunes apply the ASSUMED 0.9. The base and `default` tune keep
the brief's `C_rr` = 0.015 (surface `brief`); the 09 presets are selectable surfaces.

| Input | Default | Status |
|---|---|---|
| `C_rr` | 0.06 (dirt) | ASSUMED; presets pavement 0.03, dirt 0.06, grass 0.10, sand 0.20 |
| `C_dA` | 0.9 m² | ASSUMED (range 0.8–1.2) |
| `η` | 0.90 | ASSUMED (gearbox + chain 0.88–0.95) |
| `J_w` (4 wheels) | 1.2 kg·m² | ASSUMED (0.25–0.45 each) |
| `μ_t` | 0.7 | ASSUMED (dirt 0.6–0.8) |
| `s0` | 0.08 | ASSUMED slip at about 76 % of peak force |
| `v_min` | 0.5 m/s | numerical |
| `b/L_wb`, `h_cg`, `L_wb` | 0.60, 0.50 m, 1.50 m | PLACEHOLDER |
| governor | 3800 rpm, 100 rpm linear cut | team rule / ASSUMED |

## 4. Belt heating

Lumped belt temperature: `C_b dT_b/dt = f_b (P_slip,p + P_slip,s) + P_hyst − hA (T_b − T_amb)`
(integrated exactly for a constant source over each step; `T_b` in K; heat leaves the mechanical
energy through `W_diss`, so `T_b` does not enter the energy balance),
where `C_b = m_belt c_p` and `m_belt = ρA L`. `P_hyst` is the belt Kelvin–Voigt dissipation,
and `f_b = 0.5` is the belt's share of flank slip heat. `c_p = 1700 J/(kg·K)`, `h = 40 W/(m²K)`,
`A = L (b + 2 h_c)`, and `T_amb = 40 °C` are all ASSUMED. Step response:
`T_b(t) = T_amb + (P/hA)(1 − e^{−t/τ})` with `τ = C_b/hA`.

## 5. Tune library

`cvtsim/io/tuning/components.yaml` holds named parts with sources: `primary_springs`,
`flyweights`, `ramps` (profile tables), `secondary_springs`, `torsion_springs`,
`helixes` (angle tables), `vehicle_surfaces`. A tune `cvtsim/io/tuning/tunes/<name>.yaml`
selects one of each plus scalar adjustments. `default.yaml` reproduces the brief's
parameters exactly. CLI:

    python -m cvtsim.tune run default [--t-end 8]
    python -m cvtsim.tune compare tuneA tuneB ...

Each tune run is a tracked run (group `tune`) with its dashboard.

## 6. Output contract (channels written by the simulator)

Existing `CH`, plus: `theta_e, theta_s` (shaft angles), `r_w, F_fw, F_spr_p, F_spr_s,
F_tors, F_cam, cap_p, cap_s` (traction capacity `2 μ r_eff F_a / cos β`), `P_slip_p,
P_slip_s, T_belt, v, w_wheel, F_tire, N_rear, P_ground, a_x`. Definitions: `F_spr_s = F_l0 + k_l x`,
`F_tors = k_t (θ_t0 + h) h'` (axial share of the torsion spring), `F_cam = −h' τ_mov` (same as the
existing `Fcam`, kept for the contract), `T_belt` in K, `P_ground = F_x v`. In the rigid model
`F_tire = (F_roll + F_aero) + m a_x`. Belt snapshots:
node positions and element tensions at `n_frames` (default 400) evenly spaced output times.

## 7. Test targets (computed in closed form by targets.py)

| ID | Check |
|---|---|
| T9.1 | Straight-ramp table with κ0 = 0.5 reproduces 05 targets (engagement 1725.66 rpm) |
| T9.2 | Curved ramp `r_w = r_w0 + 0.4u + 4u²`: `F_fw(u=0.01, 3000 rpm)` = closed form |
| T9.3 | Free primary with curved ramp (no damping or stops): energy conserved < 1e-6 over 1 s (RK4). With this ramp the brief spring has no stable equilibrium (net stiffness < 0 wherever the flyweights balance the preload), so the test uses a test spring (40 kN/m, 568 N) that keeps the motion inside the smooth profile |
| T9.4 | Two-angle helix (45° to 10 mm, then 35°): `h(x)` and `F_s(15 mm, 20 N·m)` = closed form |
| T9.5 | Free secondary with progressive helix: energy conserved < 1e-6 |
| T9.6 | Fixed-ratio steady speed with `C_rr`, `C_dA`, `η`: no-slip quadratic root `c_a v² + η b k² v + (C_rr m g − η a k) = 0`; with the tire the wheel turns at `v(1+s)`, `s = s0 artanh(F/(μ_t N_r))`, which lowers `v` by 1.7 % (18.659 vs 18.977 m/s for a=40, b=0.05, R=1, dirt, η 0.9, C_dA 0.9). The simulation is checked against the slip-corrected root (brentq in targets.py) |
| T9.7 | Tire: below traction `F_x` tracks demand; above `μ_t N_r` the wheel spins (slip > 0.5) |
| T9.8 | Belt heating step response: `T_b(τ)` = `T_amb + 0.632 P/hA` |
| T9.9 | `default` tune reproduces the pre-library launch metrics to 1e-6 when the new vehicle terms are off (every stored channel of the 0.15 s launch from commit 7268716, `tests/data/`, to 1e-6 of its range) |
