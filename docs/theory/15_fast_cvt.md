# 15 — Fast CVT model for full events (`cvtsim.fast`)

The segmented-belt model (03–11) takes about 4 min of compute per second of driving, which rules
out full events such as a 4-hour endurance. The fast model keeps the same clutch parts, engine,
tire and vehicle, but replaces the belt with:
* a ratio set by belt length closure;
* a clamp-force balance;
* a stick–slip belt clutch.

A few constants that the reduced belt cannot derive are calibrated against the full model (§5).

## 1. State
Engine speed `ω_e`, secondary speed `ω_s`, vehicle speed `v`, primary travel `x_p`, and belt
temperature `T_b`. The secondary travel `x_s` follows from `x_p` through belt closure. The belt
is inextensible and massless.

## 2. Geometry (05 §2)
For `x_p ≥ 0`, measured from belt contact:
* `r_p = r_p,low + x_p / (2 tan β)`;
* `r_s` closes the belt length `L` at centre distance `C`;
* `x_s = 2 tan β (r_s,max − r_s)`;
* wrap angles `φ_p = π − 2 asin((r_s − r_p)/C)` and `φ_s = 2π − φ_p`.

For `x_p < 0`, the sheave is in the engagement gap. The belt rides the hub bearing and transmits
no torque.

`r_s(x_p)`, `φ_p` and `φ_s` are tabulated on 4001 points and interpolated linearly.

## 3. Clamp forces
**Primary.** The flyweight force on the ramp minus the spring, as in 05/09:
    F_fw = n m_w r_w(x_p) κ(x_p) ω_e²,    F_spr = F_p0 + k_p x_p,    F_p = F_fw − F_spr

**Secondary.** Linear spring, torsion spring and helix feedback (05 §4, 09 §2), with `h(x_s)` the
helix twist:
    F_s = F_l0 + k_l x_s + k_t (θ_0 + h) h' + s_h h' T_s
* `T_s` is the belt torque on the secondary.
* `s_h` is the share carried by the movable half. The full model gives 0.50–0.56 over a launch;
  the default is 0.5.

**Belt tension.** A V-belt in a groove of half-angle `β`, over wrap `φ`, with mean tension `T̄`,
presses on a sheave with axial force `T̄' φ / (2 t_w)`:
* `T̄' = T̄ − T_c`, where `T_c = ρ A v_b²` is the centrifugal tension.
* Without friction in the radial balance, `t_w = tan β`.
* Friction as the belt wedges radially raises it. With the 12° flank this term dominates: the full
  model gives `T̄' φ / (2 F_s)` = 0.56 to 0.97 over a launch, against tan β = 0.21.
* The fast model uses `t_w = tan β + μ (c_s + c_s2 x_p / x_p,stop)`, with `c_s` and `c_s2`
  calibrated.

The secondary sets the tension unless the primary cannot supply it:
    T̄' = min(2 t_w F_s / φ_s,  2 t_w F_p⁺ / (g φ_p)),    F_p⁺ = max(F_p, 0)

**Thrust ratio.** The primary thrust that holds the ratio is
    F_p* = g F_s φ_p / φ_s
* `g = 1` at no load.
* Under torque, the driving sheave needs more thrust: the full model shows `F_p/F_s` of 1.5–2.0.
* Calibrated form: `g = 1 + g_1 τ + g_2 τ²`, with torque ratio `τ = |F_b| / F_cap`.

**Belt pull capacity.** At gross slip, the capstan relation in a groove gives
    F_cap = 2 T̄' tanh(μ' φ_min / 2),    μ' = c_μ μ / sin β,    φ_min = min(φ_p, φ_s)
`c_μ` is a calibration factor, 1 by default.

## 4. Dynamics
**Belt clutch.** The belt pull `F_b` couples the shafts:
    J_e ω̇_e = T_eng − F_b r_p,    J_s ω̇_s = F_b r_s − T_load
* `J_e = J_f + J_c + n m_w r_w0²`.
* `J_s` and `T_load` are as in the full model, with tire slip, driveline efficiency, grade and
  traction factor (09, 14).
* **Stick:** the constraint `r_p ω_e = r_s ω_s` gives
      F_b,stick (r_p²/J_e + r_s²/J_s) = r_p T_eng/J_e + r_s T_load/J_s + ṙ_p ω_e − ṙ_s ω_s
* The belt sticks while `|F_b,stick| ≤ F_cap`. Otherwise it slips, with
  `F_b = F_cap sign(r_p ω_e − r_s ω_s)`. It re-sticks when the slip speed changes sign and the
  stick force is within capacity.

**Shift law.** The ratio changes at a rate set by the thrust imbalance and the belt speed:
    ẋ_p = K ω_e D(F_p − F_p*, F_f) / F_s
* `D` removes the friction band `F_f = μ_r F_fw + μ_b 2T̄`.
* `K` is calibrated (§5).
* The travel is clipped to `[0, x_p,stop]`.
* In the engagement gap, the movable half moves with its own mass and the flyweights' mass,
  `M = m_ph + n m_w κ²`, under `F_p − c_p ẋ_p` less the roller friction band (05 P2'). It reaches
  the belt at `x_p = 0`.

**Engine** torque, throttle, motoring torque, schedule, tire, road load and belt heating are as in
06, 09 and 14. Belt heat is `|F_b (r_p ω_e − r_s ω_s)|`.

**Integration:** explicit, with step `h = 2e-4 s`.

## 5. Calibration
* **Fitted constants:** `K`, `g_1`, `g_2`, `c_s`, `c_s2`.
* **Data:** full-model runs listed in `cvtsim/fast/calibration.yaml`.
* **Objective:** squared errors, each normalised, in:
  * engine rpm (after engagement);
  * primary travel;
  * vehicle speed.
* **Validation:** full-model runs not used in the fit, including a different spring and a
  throttle lift. A setup comparison is only trusted if the fast model passes on parts it was not
  fitted to.

Acceptance criteria, set before fitting:

| Quantity | Limit |
|---|---|
| Engine rpm after engagement, RMS error | < 5 % of its mean |
| Primary travel, RMS error | < 1.5 mm |
| Vehicle speed at the end | within 3 % |
| Engagement rpm | within 3 % |

## 6. Targets (closed forms, independent of the fit)
| ID | Case | Target |
|---|---|---|
| T15.1 | Geometry at `x_p = 0` and at `x_p,stop` | `r_s = r_s,low`, `r_s/r_p` = low ratio, and 0.8 at the stop (05) |
| T15.2 | Wrap angles at `r_p = r_s` | `φ_p = φ_s = π` |
| T15.3 | Capacity, `T̄' = 500 N`, `φ = π`, `μ = 0.4`, `β = 12°` | `F_cap = 2·500·tanh(0.4π/(2 sin 12°))` = 995.268 N |
| T15.4 | Locked belt, constant `T_eng`, no load | `ω̇_e = T r_s² / (J_e r_s² + J_s r_p²)` |
| T15.5 | Engagement against the full model | within 3 % (calibration check, section 5) |
| T15.6 | Energy | `ΔE_kin = W_in − W_diss`, within 1e-3 of input work over a launch |
| T15.7 | Step | `h` vs `h/2`: rpm and speed at 8 s within 0.5 % |

## 7. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A15.1 | Massless, inextensible belt; no creep | ? (slip below gross sliding not modelled) |
| A15.2 | Thrust ratio `g(τ)` fitted to the full model, not derived | ? (valid within the calibrated range) |
| A15.3 | Shift rate proportional to belt speed and thrust imbalance (the form of Ide's ratio-change law for metal belts, assumed for rubber) | ? |
| A15.4 | Constant flyweight inertia | ↓ small, `n m_w r_w²` varies about 2 % of `J_e` |
| A15.5 | No firing pulses or PTO torsion (mean torque map) | ↓ peak shaft torque; use the full model for peaks |
