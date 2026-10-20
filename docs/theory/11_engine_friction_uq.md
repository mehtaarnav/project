# 11 — Engine firing pulses, PTO torsion, clutch friction, uncertainty bands

Extends 05, 06 and 09. Every new input carries a nominal value and a range `[low, high]`
with a source. Results are reported as bands (§4).

## 1. Single-cylinder torque pulses

Crank angle `θ` (0 = TDC firing, 720° cycle). Crank radius `r = S/2`, rod length `l`,
`λ = r/l`, piston area `A = π B²/4`. Piston displacement from TDC:

    s(θ) = r (1 − cos θ) + l (1 − sqrt(1 − λ² sin² θ))
    s'(θ) = r sin θ (1 + λ cos θ / sqrt(1 − λ² sin² θ))          (kinematic lever arm)

**Gas torque.** Cylinder gauge pressure `p(θ)` comes from a single-zone cycle: isentropic
compression and expansion (`n = 1.30`), Wiebe heat release
`x_b = 1 − exp(−a ((θ − θ_s)/Δθ)^(m+1))` with `a = 5`, `m = 2`, start `θ_s = −15°`,
duration `Δθ = 50°`, and intake and exhaust at ambient. Volume:
`V(θ) = V_c + A s(θ)` with `V_c = V_d/(r_c − 1)`. Then `T_gas(θ) = p(θ) A s'(θ)`.

**Mean-torque calibration.** At each speed the pulse is scaled so that its cycle mean equals
the dyno (brake) torque:

    T_eng(θ, ω) = T_map(ω) · T_gas(θ)/⟨T_gas⟩ + T_rec(θ, ω)

so the work per cycle is the map's. `T_rec` is the reciprocating-inertia torque (zero mean):

    T_rec = − m_rec s'(θ) [ s''(θ) ω² ]      (m_rec = piston + pin + rings + 1/3 rod)

The Wiebe parameters set the pulse shape. Their effect on the PTO alternating torque is part
of the uncertainty screen (§4).

Cycle phases (code, `cvtsim/drivetrain/pulses.py`): expansion `θ ∈ [0, π]`, exhaust `[π, 2π]`,
intake `[2π, 3π]`, compression `[3π, 4π]`. The valves are taken to open and close at the dead
centres; exhaust blowdown at BDC is instantaneous (`s'(π) = 0`, so the torque stays continuous).

### 1.1 Heat release closure (spec correction: `Q_in` was unspecified)

The text above fixes the shape parameters but not the heat released per cycle `Q_in`, and the
shape `T_gas/⟨T_gas⟩` depends on it (the motored compression torque does not scale with load).
Derivation: on the closed part, `ψ = θ` or `θ − 4π`, `ψ ∈ [−π, π]`, the first law with internal
energy `U = pV/(n − 1)` is `δQ = dU + p dV = (V dp + n p dV)/(n − 1)`, i.e.

    d(p V^n) = (n − 1) V^(n−1) dQ,      dQ = Q_in x_b'(ψ) dψ
    p(ψ) = [p_amb V1^n + (n − 1) Q_in I(ψ)] / V(ψ)^n,   I(ψ) = ∫_{−π}^{ψ} V^(n−1) x_b' dψ'

`p` is linear in `Q_in`: `p = p_mot + Q_in p_Q`. The motored part `p_amb (V1/V)^n` is even in `ψ`
and `s'` is odd, so its torque has zero cycle mean. The work per unit heat follows from the first
law over the closed part (`W = Q − ΔU`, `ΔU = p_Q(π) V1/(n − 1)` per joule):

    W_Q = ∫_{θ_s}^{π} x_b' [1 − (V/V1)^(n−1)] dψ = 0.442648     (indicated efficiency; Otto 1 − r_c^(1−n) = 0.464)

Closing `Q_in` so that the indicated work equals the map work, `Q_in W_Q = 4π T_map`, makes the
spec's scale factor exactly one and gives

    T_gas(θ, T_m) = g_mot(θ) + T_m g_Q(θ),   g_mot = A s' (p_mot − p_amb),   g_Q = 4π A s' p_Q / W_Q

with `⟨g_mot⟩ = 0` and `⟨g_Q⟩ = 1`, so `⟨T_gas⟩ = T_m` for any mean torque. `T_m` is the
governed map torque `T_e` of 09 §3: the governor cuts combustion, the compression pulse remains.
At 25.08 N·m: `Q_in` = 712.0 J/cycle, peak pressure 36.4 bar, `T_gas` from −100 to +359 N·m.
Friction is neglected in the closure (brake = indicated work), which ↓ `Q_in` and so ↓ the
alternating torque slightly (mechanical efficiency ≈ 0.85–0.9 would raise the combustion part by
10–15 %). The map and governor read the instantaneous crank speed `ω_e` (?: a real governor
responds to a filtered speed).

**Table.** `g_mot`, `g_Q` are tabulated on `n` uniform nodes over 4π (`engine.pulse_table_n` =
49152, 0.0146°) and interpolated linearly in the jitted loop. `n` is a multiple of 48 so that nodes
fall on BDC and `θ_s`; `I(ψ)` uses 8-point Gauss–Legendre per cell. The tables are normalized by
their own trapezoid means (`g_Q` to 1, `g_mot` to 0), which are the exact integrals of the
piecewise-linear interpolant, so the cycle mean in the loop is `T_m` to rounding (T11.2). The
unnormalized trapezoid mean converges to the continuum `W_Q/(4π)` = 0.0352248 N·m/J as `h²`:

| `n` | resolution | `|mean/(W_Q/4π) − 1|` |
|---|---|---|
| 3072 | 0.234° | 1.5e-7 |
| 12288 | 0.0586° | 9.7e-9 |
| 24576 | 0.0293° | 2.4e-9 |
| 49152 | 0.0146° | 6.0e-10 (used) |

### 1.2 Reciprocating inertia (spec correction: energy)

`T_rec = − m_rec s' s'' ω²` alone, on a constant `J_f`, is not conservative: the reciprocating
mass stores `½ m_rec (s' ω)²`. Lagrange with `𝒯 = ½ (J_f + m_rec s'²) ω²`
(`∂𝒯/∂θ = m_rec s' s'' ω²`, `d/dt(∂𝒯/∂ω) = J(θ) ω̇ + 2 m_rec s' s'' ω²`) gives

    (J_f + m_rec s'(θ_e)²) ω̇_e = T_gas − T_pto + T_rec,     T_rec = − m_rec s' s'' ω_e²

so `T_rec` is the spec's term and the inertia carries the `s'² ω̇` part. Energy rate is
`T_gas ω_e − T_pto ω_e` exactly. `⟨T_rec⟩ = 0` over a revolution because `s' s'' = d(s'²/2)/dθ`.
`s'' = r cos θ + r λ [cos 2θ / D + λ² sin²θ cos²θ / D³]`, `D = sqrt(1 − λ² sin² θ)`.

`engine.pulses: false` is the mean-torque engine of 09: `T_eng = T_map`, no `m_rec`. Its cycle
mean inertia `⟨m_rec s'²⟩` = 2.75e-4 kg·m² (0.55 % of `J_f + J_c`) is omitted there (↓ inertia).

| Input | Nominal [range] | Source |
|---|---|---|
| bore `B`, stroke `S`, compression `r_c` | 89.0 mm, 69.0 mm, 8.0 | Kohler CH440 published specs (engine-specs.net, listings). Baja SAE U.S./Canada accept only Rehlko (Kohler) PA-CH440-3302 (with BSAE kit) and PA-CH440-3337 |
| `n`, Wiebe `a`, `m`, `p_amb` | 1.30, 5, 2, 101325 Pa | ASSUMED (text above; ISA sea level) |
| rod length `l` | 105 mm [95, 115] | ASSUMED (λ ≈ 0.33; service manual not reachable) |
| `m_rec` | 0.45 kg [0.35, 0.60] | ASSUMED for an 89 mm piston + pin + rings + 1/3 rod |
| Wiebe `θ_s`, `Δθ` | −15° [−25, −5], 50° [40, 65] | ASSUMED typical SI values |

## 2. PTO torsional compliance and inertia split

The engine rotor splits into the crank + flywheel `J_f` (angle `θ_e`, receives `T_eng`) and the
primary clutch `J_c` (angle `θ_c`, receives the belt torque and the flyweight terms of P1′).
They are joined by the PTO torsional spring and damper:

    T_pto = k_t,pto (θ_e − θ_c) + c_t,pto (ω_e − ω_c),    k_t,pto = G π d⁴ / (32 L_eff)

with `d = 25.4 mm` (PTO, consistent across listings) and `L_eff` = PTO length 88.6 mm plus an
equivalent crank-web length, ASSUMED 30 mm [10, 60]. `c_t,pto` gives 2 % modal damping
(ASSUMED). The keyway reduces `J_p`: use `d_eff = d − w_key/2`, a standard approximation.
`T_pto(t)` is the dynamic PTO torque used for fatigue.

| Input | Nominal [range] | Source |
|---|---|---|
| `J_f` (crank + flywheel + magneto) | 0.035 kg·m² [0.025, 0.050] | ASSUMED; ask Rehlko |
| `J_c` (primary fixed half + spider + cover) | 0.015 kg·m² [0.010, 0.025] | PLACEHOLDER; weigh and CAD the clutch |
| PTO `d`, length | 25.4 mm, 88.6 mm | CH440 1 in × 3.49 in straight PTO (listings, e.g. PA-CH440-3275; no public drawing of -3302/-3337 found) |
| key width `w_key` | 6.35 mm | DERIVED: ASME B17.1 standard key for a 1 in shaft (1/4 in); UNVERIFIED on the engine |
| shear modulus `G` | 76.9 GPa | DERIVED `E/(2(1+ν))` from `shafts.material` (PLACEHOLDER generic steel) |

**Equations of motion.** The crank obeys §1.2. The clutch rotor carries (P1′) of 09 with `J_e`
replaced by `J_c` and the engine torque by `T_pto`:

    (J_c + n m_w r_w²) ω̇_c + 2 n m_w r_w κ u̇ ω_c = T_pto + τ_belt,p          (P1″)

(P2′) is unchanged with `ω → ω_c`; belt friction on the primary uses `ω_c`. Energy adds
`½ (J_f + m_rec s'²) ω_e² + ½ (J_c + n m_w r_w²) ω_c² + ½ k_t,pto (θ_e − θ_c)²`, and
`c_t,pto (ω_e − ω_c)²` goes into `W_diss`. `J_e = J_f + J_c` = 0.05 kg·m² (the 05/09 placeholder)
is kept as a derived value for the one-rotor model.

**Derived values.** `d_eff` = 22.225 mm, `L_eff` = 118.6 mm, `k_t,pto` = 15 536 N·m/rad.
`c_t,pto = 2 ζ sqrt(k J_red)`, `J_red = J_f J_c/(J_f + J_c)`, ζ = 0.02: 0.511 N·m·s/rad.
`f = sqrt(k (1/J_f + 1/J_c))/(2π)` = 193.6 Hz; with the flyweights on the clutch
(`n m_w r_w²`) 187.3 Hz at full open, 186.0 Hz at belt contact and 183.1 Hz at the 0.8 stop.

**Time step.** At `h` = 1.25e-7 s the PTO mode has 41 000 steps per period (≥ 20 required) and
`h ω_n` = 1.5e-4, far inside the symplectic-Euler limit 2. The belt segments still set the step;
`sim.dt` is unchanged.

### 2.1 Rigid limit

`engine.pto_model: rigid` is the `k_t,pto → ∞` limit: one rotor `J_f + m_rec s'² + J_c + n m_w r_w²`
with `θ_c ≡ θ_e`; `T_pto` is then reported as the transmitted torque `T_gas + T_rec − J_f(θ) ω̇_e`.
A finite but large `k` cannot replace it in an explicit scheme (`h < 2/ω_n`). With pulses off, the
0.15 s launch approaches the rigid one as `k` grows (largest crank-speed difference over the run,
relative to the peak): 8.6e-4 at `k`, 4.4e-4 at 4`k`, 2.1e-4 at 16`k` (about `1/sqrt(k)`, i.e.
`1/ω_n`; slow test `test_pto_rigid_limit`). With pulses off, the rigid PTO and zero clutch friction
the model reproduces the pre-spec-11 launch exactly (`test_mean_torque_model_reproduced`).

## 3. Clutch friction (Coulomb, regularized)

All friction forces use `tanh(v/v_x)` with `v_x = 1 mm/s` (numerical), and their power is
dissipated (`W_diss`).

* **Helix button.** The helix angle `α` is the one of 09 §2, `h' = 1/(r_c tan α)`. (Spec
  correction: that relation and the tangent below put `α` from the circumferential direction,
  not from the shaft axis; the circumferential slide per unit axial travel is `r_c h' = cot α`.
  The two readings agree at 45°.) The path unit tangent in the
  (axial, circumferential) plane is `(sin α, cos α)`. The contact normal is `N = |τ_c|/(r_c sin α)`,
  where `τ_c = τ_mov + k_t(θ_t0 + h)` is the torque carried through the cam. The friction virtual
  work over `δs = δx/sin α` gives

      Q_x,f = − μ_h |τ_c| / (r_c sin² α) · tanh(ẋ_s/v_x)

  At 45° this is ±2 μ_h of the frictionless cam force, i.e. ±20 % at μ_h = 0.1.
  For a profiled helix the code uses `1/sin² α = 1 + (r_c h'(x))²`. A rigid rotation of the
  secondary does not slide the button, so the friction has no `θ_s` component.
* **Primary ramp rollers.** `Q_u,f = − μ_r |F_fw| tanh(u̇/v_x)` (rolling plus pin friction lumped),
  `F_fw = n m_w r_w(u) κ(u) ω_c²`.
* **Sheave bushings, both clutches.** `Q_f = − μ_b |F_hub| tanh(ẋ/v_x)`, where `F_hub` is the
  in-plane belt resultant on that sheave (flank contact only; the primary hub-bearing load is
  carried by the shaft, not the bushing). Assumption (?): only the movable half's share (about ½
  for symmetric flanks) passes through its bushing, but its tipping moment `F_a r/L_b` adds edge
  loads not modelled; the full resultant is a lumped stand-in, covered by the `μ_b` range.

The friction enters (P2′) as `+Q_u,f` and the secondary `Q_x` as `+Q_x,f + Q_f`; the powers
`−Q ẋ ≥ 0` are added to `W_diss`. Channels: `F_fric_p = Q_u,f + Q_f,p`, `F_fric_s = Q_x,f + Q_f,s`.
Regularization stiffness `c = μ|F|/v_x` reaches about 5e5 N·s/m; with the axial masses (≈ 1–1.8 kg)
the explicit limit `h < 2M/c` ≈ 4e-6 s is 30× the launch step.

| Input | Nominal [range] | Source |
|---|---|---|
| `μ_h` helix button | 0.10 [0.05, 0.20] | ASSUMED (polymer/bronze button on aluminium helix) |
| `μ_r` roller + pin | 0.05 [0.02, 0.10] | ASSUMED |
| `μ_b` bushing | 0.10 [0.05, 0.15] | ASSUMED (PTFE/bronze bushing) |

## 4. Uncertainty bands

`cvtsim/io/uncertainty.yaml` lists every uncertain input as `{nominal, low, high, source}`
under its parameter path. The belt friction μ, `v_eps`, `k_n`, the damping values, the inertias,
drag, rolling resistance and the inputs above all belong in it.

1. **Screen (one at a time).** Run each input at its low and high values, with all others nominal.
   For every output metric, report `(low − nominal)/nominal` and `(high − nominal)/nominal`
   as a tornado chart sorted by span.
2. **Bands.** For the top `k` inputs (default 5, chosen by the largest normalized span across the
   decision metrics), run a Latin-hypercube sample of `n` runs (default 30), uniform on each range.
   Report the median and the 5–95 % band of each metric.
3. **Decision metrics.** 100 ft time; engagement rpm; shift-complete time; peak and alternating
   PTO torque; peak PTO-side radial (hub) load; peak jackshaft bearing loads; peak belt tension;
   belt temperature.

Runs execute in parallel processes, each a tracked run (group `uq`).

## 5. Test targets

| ID | Check |
|---|---|
| T11.1 | Slider-crank kinematics: `s(180°) = 2r`; `s'` matches the finite difference of `s` to 1e-9 (and `s''` that of `s'`), relative to their peaks |
| T11.2 | Pulse mean: `⟨T_eng⟩` over 720° = `T_map(ω)` to 1e-9 at 2400, 3000, 3600 rpm. Table: unnormalized mean = `W_Q/(4π)` (targets.py, adaptive quadrature) to 1e-9 at `pulse_table_n`, `h²` convergence |
| T11.3 | `⟨T_rec⟩ = 0` over a revolution to 1e-12 (of its peak) |
| T11.4 | Energy closure with pulses, PTO spring and friction: < 1e-4 of W_in over 0.5 s (measured 1.7e-6) |
| T11.5 | Two-inertia torsional frequency `sqrt(k (1/J_f + 1/J_c))/(2π)` from a free-vibration run, to 1e-3. The run has the primary open, so `J_c` includes `n m_w r_w(−g)²`: 187.32 Hz (measured error 4.4e-5); the spec value 193.60 Hz is checked on `derived.pto.f_n` |
| T11.6 | Helix friction, quasi-static at 45°: closing force = (τ_c/r_c)(1 + 2μ_h) while opening (upshift) and (τ_c/r_c)(1 − 2μ_h) while closing (backshift); ratio (1+2μ)/(1−2μ) = 1.5 at μ_h = 0.1 |
| T11.7 | UQ screen: with a synthetic linear model, the tornado reproduces the analytic sensitivities exactly |
| T11.8 | `pulses: false`, `pto_model: rigid`, all clutch μ = 0 reproduce the pre-spec-11 launch (commit 3eea3d6, `tests/data/prespec11_launch_0p15s.csv.gz`), every channel, to 1e-9 of its range |
| T11.9 | (slow) rigid limit: crank-speed difference to the rigid rotor falls monotonically for `k` × 1, 4, 16 and is < 1e-3 at 16 |

## 6. Implementation contract (sections 1–3)

* `cvtsim/drivetrain/pulses.py`: `crank_kinematics(θ, r, λ) → (s, s', s'')`,
  `rec_torque(θ, ω, m_rec, r, λ) → (T_rec, m_rec s'²)`, `pulse_table(engine, n) → dict(gm, gq, mean_raw, W_Q, …)`,
  `gas_torque(θ, T_m, gm, gq)`, `work_per_heat(engine)`.
* `cvtsim/cvt/system.py`: state adds `w_c, theta_c`; `w_e, theta_e` are the crank.
  `helix_friction`, `coulomb` give the §3 forces. Parameters: `engine.{pulses, pto_model, J_f, J_c,
  bore, stroke, compression_ratio, rod_length, m_rec, wiebe_*, n_poly, p_amb, pulse_table_n, pto_d,
  pto_length, crank_web_length, key_width, pto_damping_ratio}`, `primary.mu_roller`,
  `secondary.mu_helix`, `sheave.{mu_bushing, v_fric}`; `derived.pto = {G, d_eff, L_eff, k, c, f_n}`.
* Channels (appended to 09 §6): `T_pto`, `theta_c`, `w_c`, `F_fric_p`, `F_fric_s`, `T_pulse`
  (`T_gas + T_rec` on the crank). `T_eng` stays the governed map (cycle-mean) torque `T_e`, so the
  dashboards' torque-vs-speed plots keep their meaning; `F_fw`, `slip_p` use `ω_c`.
  `W_in = ∫ T_gas ω_e dt` (the motored part integrates to zero over a cycle at constant speed).
