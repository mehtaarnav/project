# 13 — Belt life and CVT efficiency (`cvtsim.post`)

Post-processing of a launch result `res` (array `(n_out, len(CH))`, columns named by
`cvtsim.cvt.system.CH`). Nothing here feeds back into the simulation. Notation and frames as in
09. SI inside; rpm only in the efficiency map axis.

## 1. Sources
* No external source is used for the belt-life constants. Every constant is an input in
  `cvtsim/io/belt_life.yaml` with a tag (PLACEHOLDER, ASSUMED, numerical choice).
* Linear damage accumulation (Miner) is the rule already used in 12 §2 (`cvtsim.fe.fatigue.miner_damage`).
  Edition and section UNVERIFIED there.
* Temperature acceleration of rubber fatigue, "life halves for every 10 K above a reference": a
  rule of thumb used in rubber engineering practice. **ASSUMED**. No citation is given; the
  value is not calibrated for this belt.
* Power-law fatigue `N_f = N_ref (ε_ref/ε_b)^k`: a generic Basquin-type form. The exponent and
  reference are **PLACEHOLDER**, uncalibrated for rubber.

## 2. Power accounting (per instant)
Channel names from `CH`. `ω_e = w_e`, `ω_s = w_s`.

    P_eng     = T_eng · w_e                                   (net engine torque, after governor)
    P_belt_in = |tq_p| · |w_e|                                (belt-side power into the primary)
    P_out     = |tq_s_fix + tq_s_mov| · |w_s|                 (belt-side power out of the secondary)
    P_slip    = P_slip_p + P_slip_s                           (flank slip loss, ≥ 0)
    R_bal     = P_belt_in − P_out − P_slip                    (diagnostic: inertia, belt hysteresis)

Magnitudes are used so that the efficiency is defined for either sign of the speed; the
forward-drive sign is the normal case. Reverse power flow is out of scope.

    η_cvt = P_out / P_belt_in        where P_belt_in > P_min,   undefined (NaN) otherwise

`P_min` is a numerical threshold (1 W) so that the ratio is not formed from transient noise at
standstill. `η_cvt` is reported raw, never clipped.

**Quasi-steady mask.** `P_slip > P_belt_in` (slip power above the input power) occurs at belt
engagement: the flank slip is paid from stored belt and inertia energy, so the instantaneous ratio
has no meaning. Such instants are excluded from the map and from the `[0, 1]` range check
(`quasi_steady_mask`). The launch run has one such instant in 0.15 s (eta = 2.73, P_in = 10.6 W,
P_slip = 123 W at t = 0.1365 s); it is reported, not hidden.

**Efficiency map.** 2-D histogram of the mean `η_cvt` over the bins of
`g = r*_s / r*_p` (geometric ratio, channels `rstar_s/rstar_p`) and `n_p` (primary rpm, from `w_e`).
Only quasi-steady instants (defined `η_cvt` and `P_slip ≤ P_belt_in`) are counted. Bin mean `η̄ = Σ η / N` in each bin with `N ≥ 1`;
empty bins are NaN.

## 3. Belt bending fatigue

### 3.1 Bending cycles
Each passage of a belt element over a sheave is one bending cycle. For sheave `k ∈ {p, s}` with
effective radius `r_k = r_eff,k` and belt speed `v_k = r_k |ω_k|`:

    θ_k = π − 2 asin((r_other − r_k)/C)       wrap angle of the open belt (C = centre distance)
    L_k = r_k θ_k                             belt length in contact per sheave
    ṅ_k = v_k / L_k                           bending cycles per second

For equal sheaves `θ = π` exactly. Valid while `|r_other − r_k| < C`; outside that the
expression is clamped to `[0, 2π]` (not reached in the launch run; reported if it happens).
Where a sheave has no belt contact (`r_eff,k ≤ 0`, e.g. before the first engagement) its rate is 0.

### 3.2 Bending strain
Neutral axis to outer fibre `h_n = f_h · t_flank`, with `t_flank` = `belt.flank_height` (12 mm
PLACEHOLDER) and `f_h = 0.5` (PLACEHOLDER: neutral axis at mid-height). Over a sheave of radius `r`:

    ε_b = h_n / r

### 3.3 Damage rate
Life in cycles at strain `ε_b` (power law, PLACEHOLDER constants):

    N_f(ε_b) = N_ref (ε_ref / ε_b)^k,     k = 5 (PLACEHOLDER)

Temperature factor from the belt temperature `T_b` (channel `T_belt`, K), life halves for every
`ΔT_h = 10 K` above `T_ref` (ASSUMED, §1):

    f_T(T) = 2^{(T − T_ref)/ΔT_h}           damage multiplier (f_T = 1 at T_ref)

Damage per instant (Miner, linear), summed over sheaves:

    Ḋ = f_T(T_b) Σ_k ṅ_k / N_f(ε_k) = f_T(T_b) Σ_k ṅ_k (ε_k/ε_ref)^k / N_ref

Accumulated damage over the run by the trapezoid rule on `t`: `D = ∫ Ḋ dt`.

### 3.4 Projected life
Duty-cycle assumption (**PLACEHOLDER** `duty`): the simulated launch is representative of a
fraction `duty` of operating time, and the rest is not damaging. Failure is at `D = 1`:

    t_life = T_run / (D · duty)      [s],      life_h = t_life / 3600

with `T_run = t_N − t_0` of the result. With `duty = 1` every operating second is a launch like
this one. The value is only as good as the uncalibrated `N_ref`, `ε_ref`, `k`.

## 4. Heat fade (diagnostic only, not fed back)
A temperature-dependent friction law would give

    μ(T) = μ0 (1 − c (T − T_ref)),          c ∈ [c_lo, c_hi] (ASSUMED range, 1/K)

with `μ0 = belt.mu` (0.3, ASSUMED brief) and `T_ref` as §3.3. The output is `μ(T_b)` over the run
for `c_lo` and `c_hi`. The simulation keeps constant `μ`; these curves only show the size of the
effect that a fade law would have. Sign: `c > 0` means friction falls as the belt heats.

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A13.1 | Miner linear damage, no load-sequence effect | ? |
| A13.2 | Bending strain `h_n/r` with `r` the effective wrap radius, plane-section, no rubber cord mechanics | ↑ strain at small sheaves (r ≈ 14 mm primary min), so ↑ damage |
| A13.3 | Temperature factor `2^{ΔT/10}` applies to the whole belt temperature from `T_belt` (lumped) | ? |
| A13.4 | Belt speed per sheave from its own effective radius (slip not averaged) | ? |
| A13.5 | `N_ref, ε_ref, k` uncalibrated | ? — life numbers are relative until calibrated |
| A13.6 | `duty` constant over life | ? |

## 6. Closed-form test targets (`tests/test_post.py`, computed by hand, no `cvtsim` import)

| ID | Function | Closed form on synthetic input | Expected |
|---|---|---|---|
| T13.1 | `wrap_angle_and_length` equal sheaves | `r = R = 0.05`, `C = 0.2` | θ = π, L = 0.05π = 0.1570796 m |
| T13.2 | `wrap_angle_and_length` small sheave | `r = 0.02`, `R = 0.05`, `C = 0.2` | θ = π − 2 asin(0.15) = 2.8404561, L = 0.0568091 m |
| T13.3 | `bending_strain` | `h_n = 0.006`, `r = 0.02` | ε = 0.3 |
| T13.4 | `damage_rate`, one sheave, `ε = ε_ref`, `T = T_ref` | `v = 2 m/s`, `L = 0.1 m` | Ḋ = 20/N_ref |
| T13.5 | `temperature_factor` | `T − T_ref = +10 K`, `−10 K`, `0` | 2, 0.5, 1 |
| T13.6 | `power_terms` and `cvt_efficiency` | `w_e = 100`, `tq_p = −10`, `T_eng = 12`, `tq_s_fix = 300`, `tq_s_mov = −200`, `w_s = 5`, `P_slip = 100 + 50` | P_eng = 1200, P_belt_in = 1000, P_out = 500, η = 0.5, R_bal = 350 |
| T13.7 | `cvt_efficiency` threshold | `P_belt_in = 0.5 W` with `P_min = 1 W` | η = NaN |
| T13.8 | `efficiency_map` | edges `g: [0, 2, 4]`, `n: [0, 3000, 6000]` rpm; instants (g = 1.0, 2000 rpm, η = 0.8), (g = 1.5, 2500 rpm, η = 0.6), (g = 3.0, 4500 rpm, η = 0.9) | bin (0,0) = 0.7, bin (1,1) = 0.9, bins (0,1) and (1,0) NaN |
| T13.9 | `friction_fade` | `μ0 = 0.3`, `c = 0.002/K`, `T − T_ref = 50 K` | μ = 0.27 |
| T13.10 | `belt_damage` constant rate | `r_p = r_s = 0.05 m`, `ω_e = ω_s = 40 rad/s`, `h_n = 0.0025 m` (ε = 0.05 = ε_ref), `T = T_ref`, `t` in `[0, 2] s`, `N_ref = 10^6`, `duty = 1` | ṅ_k = 40/π s⁻¹, Ḋ = 80/(π·10^6) = 2.5464791e-5 s⁻¹, D = 5.0929582e-5, life_h = 10.908308 |
| T13.11 | `belt_damage` duty | as T13.10, `duty = 0.5` | life_h doubles of T13.10 |
| T13.12 | smoke, real launch 0.15 s | `build(load())`, `relax`, `simulate(t_end=0.15, nout=100)` | all outputs finite; at least one quasi-steady instant; `η ∈ [0,1]` on quasi-steady instants |
| T13.13 | `quasi_steady_mask` | `P_in = 1000 W`; `P_slip = 1500 W` (excluded), `P_slip = 500 W` (kept) | mask = [False, True] |

Targets are computed in the test file from the stated expressions, not from the code.

## 7. Implementation contract
* `cvtsim/post/efficiency.py`: `power_terms(res, CH)`, `cvt_efficiency(res, CH, p)`,
  `efficiency_map(res, CH, p)`.
* `cvtsim/post/belt_life.py`: `wrap_angle_and_length`, `bending_strain`, `temperature_factor`,
  `damage_rate(res, CH, p)`, `belt_damage(res, CH, p)`, `friction_fade(res, CH, p)`.
* `p` is the dict returned by `belt_life.inputs(params)`: YAML inputs plus `C`, `flank_height`,
  `mu0` taken from the simulation parameters. All SI.
