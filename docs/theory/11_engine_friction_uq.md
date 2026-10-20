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

| Input | Nominal [range] | Source |
|---|---|---|
| bore `B`, stroke `S`, compression `r_c` | 89.0 mm, 69.0 mm, 8.0 | Kohler CH440 published specs (engine-specs.net, listings) |
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

## 3. Clutch friction (Coulomb, regularized)

All friction forces use `tanh(v/v_x)` with `v_x = 1 mm/s` (numerical), and their power is
dissipated (`W_diss`).

* **Helix button.** The helix angle `α` is measured from the shaft axis. The path unit tangent in the
  (axial, circumferential) plane is `(sin α, cos α)`. The contact normal is `N = |τ_c|/(r_c sin α)`,
  where `τ_c = τ_mov + k_t(θ_t0 + h)` is the torque carried through the cam. The friction virtual
  work over `δs = δx/sin α` gives

      Q_x,f = − μ_h |τ_c| / (r_c sin² α) · tanh(ẋ_s/v_x)

  At 45° this is ±2 μ_h of the frictionless cam force, i.e. ±20 % at μ_h = 0.1.
* **Primary ramp rollers.** `Q_u,f = − μ_r |F_fw| tanh(u̇/v_x)` (rolling plus pin friction lumped).
* **Sheave bushings, both clutches.** `Q_f = − μ_b |F_hub| tanh(ẋ/v_x)`, where `F_hub` is the
  in-plane belt resultant on that sheave.

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
| T11.1 | Slider-crank kinematics: `s(180°) = 2r`; `s'` matches the finite difference of `s` to 1e-9 |
| T11.2 | Pulse mean: `⟨T_eng⟩` over 720° = `T_map(ω)` to 1e-9 at 2400, 3000, 3600 rpm |
| T11.3 | `⟨T_rec⟩ = 0` over a revolution to 1e-12 |
| T11.4 | Energy closure with pulses, PTO spring and friction: < 1e-4 of W_in over 0.5 s |
| T11.5 | Two-inertia torsional frequency `sqrt(k (1/J_f + 1/J_c))/(2π)` from a free-vibration run, to 1e-3 |
| T11.6 | Helix friction, quasi-static at 45°: closing force = (τ_c/r_c)(1 + 2μ_h) while opening (upshift) and (τ_c/r_c)(1 − 2μ_h) while closing (backshift); ratio (1+2μ)/(1−2μ) = 1.5 at μ_h = 0.1 |
| T11.7 | UQ screen: with a synthetic linear model, the tornado reproduces the analytic sensitivities exactly |
