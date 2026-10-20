# 12 — Fatigue and bearing life (`cvtsim.fe.fatigue`, `cvtsim.fe.bearings`)

## 1. Scope and sources
* ASTM E1049-85, *Standard Practices for Cycle Counting in Fatigue Analysis*. Three-point rainflow
  method. Standard number and year verified from memory; edition **UNVERIFIED**.
* J. E. Shigley, C. R. Mischke, R. G. Budynas, *Mechanical Engineering Design*, 10th ed., McGraw-Hill.
  Sections on endurance limit, Marin factors, Goodman, Basquin and Miner. Equation and table numbers
  are **UNVERIFIED**; each use is tagged below.
* R. E. Peterson, *Stress Concentration Factors*, Wiley. Notch-sensitivity constant (§2.3) as stated
  in Shigley: **UNVERIFIED**.
* ISO 281:2007, *Rolling bearings — Dynamic load ratings and rating life*. L10 and equivalent-load
  averaging. Edition **UNVERIFIED** against the standard text.
* SKF General Catalogue, 6205 dynamic load rating. **UNVERIFIED** (see `materials.yaml`).

## 2. Fatigue

### 2.1 Rainflow (ASTM E1049-85, three-point)
Turning points only (endpoints and local extrema). Scan points with a stack. For the last three
points, `X = |p_n − p_{n−1}|` (latest range), `Y = |p_{n−1} − p_{n−2}|`.
* If `X < Y`: read the next point.
* If `Y` contains the starting point (the first point is still in the stack and `len = 3`):
  count `Y` as a half cycle, mean `(p_0 + p_1)/2`, discard `p_0` only.
* Otherwise: count `Y` as a full cycle, mean `(p_{n−1} + p_{n−2})/2`, delete `p_{n−1}, p_{n−2}`,
  keep `p_n`.
* Residue: each adjacent pair left in the stack is a half cycle, range `|a − b|`.

Output: `(range, mean, count)`, count ∈ {0.5, 1}.

### 2.2 Rotating shaft, combined stress
A fixed-direction bending load on a rotating shaft gives a fully reversed bending stress each
revolution. Let `M_a` be the bending moment amplitude and `d` the diameter.
    σ_a,b = Kf · 32 M_a / (π d³),   σ_m,b = 0
Torsion `T(t)` is rainflowed. For each cycle with torque range `ΔT` and mean `T_m`:
    τ_a = 16 (ΔT/2) / (π d³),   τ_m = 16 T_m / (π d³)     (then Kfs applied)
Shigley von Mises combination (**UNVERIFIED** equation number):
    σ_a' = sqrt((Kf σ_a)² + 3 (Kfs τ_a)²)
    σ_m' = sqrt((Kf σ_m)² + 3 (Kfs τ_m)²) = sqrt(3) Kfs |τ_m|

Assumption A12.1: bending and torsion are in phase and the bending amplitude is constant
(bias: conservative for the bending part, since the peak bending moment is used for all cycles).

### 2.3 Notch sensitivity (Peterson)
    q = 1 / (1 + a / r),    a = 0.0254 (2070 / Su)^1.8  mm,   Su in MPa, r in mm
    Kf = 1 + q (Kt − 1)
**UNVERIFIED:** the constant `a` and its form for steel in bending, as stated in Shigley's Peterson
treatment. Not cross-checked against the original paper.

### 2.4 Endurance limit (Marin)
    Se = ka kb kc kd ke kf · Se',   Se' = 0.5 Su  (Su ≤ 1400 MPa)
    ka = a Su^b  (Su in MPa). Coefficients, **UNVERIFIED** (Shigley Table 6-2, memory):

| Finish | a | b |
|---|---|---|
| ground | 1.58 | −0.085 |
| machined or cold-drawn | 4.51 | −0.265 |
| hot-rolled | 57.7 | −0.718 |
| as-forged | 272 | −0.995 |

    kb = 1.24 d^(−0.107)  (d in mm), valid 2.79 ≤ d ≤ 51 mm
kc (reliability), kd (temperature), ke (stress concentration in the endurance test), kf (miscellaneous)
are function arguments, default 1. Caller supplies them; none is invented here.

### 2.5 Safety factors and damage
* Modified Goodman (Shigley): `n = 1 / (σ_a'/Se + σ_m'/Su)`.
* Yield check: `n_y = Sy / (σ_a' + σ_m')`.
* Basquin S-N line through `(1e3 cycles, 0.9 Su)` and `(1e6 cycles, Se)`, **UNVERIFIED anchor**
  (Shigley's approach). `S = A N^b`, `b = log10(Se/(0.9 Su)) / 3`, `A = 0.9 Su · 1e3^(−b)`.
* Miner's rule over a rainflow histogram: `D = Σ n_i / N_i`. Each cycle is converted to a
  fully reversed Goodman-equivalent amplitude `σ_eq = σ_a / (1 − σ_m / Su)`. Cycles with
  `σ_eq < Se` contribute zero (endurance cutoff, assumption A12.2; bias: non-conservative for
  variable-amplitude loads with many sub-Se cycles). Failure predicted at `D = 1`.

## 3. Bearing life (ISO 281)
    L10 = (C / P)^p · 1e6  revolutions,   p = 3 (ball), 10/3 (roller)
Variable spectrum with revolutions `n_i` at load `P_i`:
    P_eq = (Σ P_i^p n_i / Σ n_i)^(1/p)
Life in hours from the rpm history (time-averaged mean speed `n̄`):
    L10h = L10 / (60 n̄)

**Side load (PTO):** `ratio = max |F_r(t)| / F_allow`. `F_allow` is a function argument; the
value in `materials.yaml` is PLACEHOLDER until Rehlko supplies the allowable.

## 4. Closed-form test targets (computed without importing `cvtsim`)
Constants: ball p = 3. Shaft d = 25 mm, M = 100 N·m, Kt = 1 unless stated.

| ID | Case | Target |
|---|---|---|
| T12.1 | Rainflow, `[0, 10, −10, 10, −10, 0]` | full cycle (range 20, mean 0, 1.0); half (10, +5, 0.5); half (20, 0, 0.5); half (10, −5, 0.5) |
| T12.2 | Rainflow, sampled sine, 3 cycles of amplitude A = 10 from a peak `[10, −10, 10, −10, 10, −10, 10]` | total count 3.0 (2 full + 2 half), all ranges 20, all means 0 |
| T12.3 | Bending stress, M = 100 N·m, d = 25 mm, Kf = 1 | σ_a = 65.1899 MPa (equals shaft target T7.7) |
| T12.4 | Peterson, r = 1 mm, Su = 1000 MPa | a = 0.0940980 mm, q = 0.913995, Kf (Kt = 2) = 1.913995 |
| T12.5 | Peterson, r = 2 mm, Su = 1000 MPa | q = 0.955065, Kf (Kt = 2) = 1.955065 |
| T12.6 | Goodman, σ_a = 100, σ_m = 50, Se = 300, Su = 600 MPa | n = 2.4 |
| T12.7 | Yield, Sy = 500, σ_a = 100, σ_m = 50 MPa | n_y = 3.333333 |
| T12.8 | L10, constant load, C/P = 10, p = 3 | 1e9 rev; at 1000 rpm, 16666.67 h |
| T12.9 | P_eq, 75 % of revs at 1000 N, 25 % at 2000 N, p = 3 | 1401.0197 N |
| T12.10 | Marin ka (machined, Su = 600 MPa), kb (d = 25 mm) | ka = 0.827878, kb = 0.878703, Se = 218.2377 MPa (ka·kb·300) |
| T12.11 | Basquin, Su = 600, Se = 300 MPa | A = 972.000 MPa, b = −0.0850908; N(Se) = 1e6 |
| T12.12 | Miner, 1000 cycles at σ_a = 400 MPa, σ_m = 0, Su = 600, Se = 300 | N = 34017.44, D = 0.0293967 |
| T12.13 | PTO side load, history `[100, −500, 300]` N, F_allow = 1000 N | ratio = 0.5 |
| T12.14 | Shaft combined, M_a = 100 N·m, torque history `[0, 50, −50, 50, −50, 0]`, d = 25 mm, Kf = Kfs = 1 | full cycle (range 100 N·m, mean 0) gives σ_a' = 71.0390 MPa, σ_m' = 0 (equals shaft target T7.7 σ_vm) |

Tolerance: rel 1e-5 for closed forms above unless stated. Rainflow and counts: exact.
Numbers are closed forms; no tolerance is loosened.

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A12.1 | Constant bending amplitude over the load history (peak moment) | ↑ conservative |
| A12.2 | Sub-Se Goodman-equivalent cycles ignored in Miner | ↓ non-conservative |
| A12.3 | Marin kc, kd, ke, kf supplied by caller (default 1) | ? |
| A12.4 | Basquin anchor 0.9 Su at 1e3 cycles, UNVERIFIED | ? |
| A12.5 | Bearing C values are PLACEHOLDER until part numbers are confirmed | ? results illustrative |
| A12.6 | Rpm and load histories are uniformly sampled; mean rpm is the sample mean | ? |

## 6. Implementation contract
* `cvtsim/fe/fatigue.py`: `turning_points`, `rainflow`, `peterson_a_mm`, `notch_sensitivity`,
  `fatigue_kf`, `marin_ka`, `marin_kb`, `endurance_limit`, `shaft_cycle_stresses`, `goodman_n`,
  `yield_n`, `basquin_constants`, `cycles_to_failure`, `miner_damage`.
* `cvtsim/fe/bearings.py`: `l10_revolutions`, `equivalent_load`, `mean_rpm`, `l10_hours`,
  `pto_side_load_ratio`.
* Constants: `cvtsim/io/materials.yaml` (material, bearing, PTO allowable), each with a `source:` tag.
  The Marin `ka` table is held in `fatigue.py` as UNVERIFIED (see §2.4) pending a YAML move.
