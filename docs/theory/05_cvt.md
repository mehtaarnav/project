# 05 — CVT primary and secondary (`cvtsim.cvt`)

## 1. Sources
* Lagrange's equations with variable inertia (flyweights) and holonomic
  helix constraint: Haug (1989), and any analytical dynamics text, e.g. H.
  Goldstein, *Classical Mechanics*, 3rd ed., Ch. 1–2 (**UNVERIFIED** chapter
  numbers). All formulas are derived below.
* N. Srivastava & I. Haque, "A review on belt and chain continuously variable
  transmissions (CVT): Dynamics and control", *Mech. Mach. Theory* 44(1)
  (2009) 19–41. Verified record. Context: rubber-belt CVTs use
  centrifugal primaries and torque-sensing secondaries. Section numbers are
  **UNVERIFIED**.
* Parameter values: the brief (team data). Status is per item in
  `cvtsim/io/dr26_default.yaml`.

## 2. Geometry (derived)
Belt length of an open drive with pitch radii `r1` (primary), `r2`
(secondary), center distance `C`, and `φ = asin((r2 − r1)/C)`:

    L(r1, r2, C) = 2 C cos φ + π (r1 + r2) + 2 φ (r2 − r1)                (G1)

Derivation: two tangent spans of length `C cos φ`. Wrap angles are
`π − 2φ` (primary) and `π + 2φ` (secondary).

Sheave axial travel → free radius (04 §2, C1):
`r*_p = r_p,min + x_p/(2 tan β)` and `r*_s = r_s,max − x_s/(2 tan β)`.

Given ratio limits `i_low = 3.0`, `i_high = 0.8` (`i = r_s/r_p`), strokes
`s_p`, `s_s`, and `Δ_p = s_p/(2tanβ)`, `Δ_s = s_s/(2tanβ)`:

    r_s,max = 3.0 r_p,min,   r_s,min = 0.8 (r_p,min + Δ_p),   r_s,max − r_s,min = Δ_s
    ⇒ r_p,min = (Δ_s + 0.8 Δ_p)/(3.0 − 0.8)

C then comes from `L(r_p,min, r_s,max, C) = L(r_p,max, r_s,min, C)`.

| Derived (β = 45°) | Value |
|---|---|
| `r_p,min`, `r_p,max` | 10.1964 mm, 22.9364 mm |
| `r_s,max`, `r_s,min` | 30.5891 mm, 18.3491 mm |
| `C` | 251.495 mm |
| `L` | 632.776 mm |

**DISCREPANCY (see README):** these radii are implausibly small. The cause is the
45° half-angle combined with ~25 mm strokes. Either `β` or the strokes need
confirming from CAD. The code derives geometry from YAML inputs, so correcting
`beta` or the strokes there propagates automatically. Optional YAML overrides
`r_p_min`, `C` exist and must be labeled with a source when used.

## 3. Primary (centrifugal, flyweights on straight ramps)

Generalized coordinates: shaft angle `θ_p` (engine crank, rigid) and
closing travel `x_p ∈ [0, s_p]`. `n` flyweights of mass `m_w`
(weight + cover) sit at radius `r_w(x_p) = r_w0 + κ x_p`. Straight ramp ⇒
constant kinematic ratio `κ = dr_w/dx_p`. A roller wedged between two
opposed faces, each at `α` from the radial plane, separates them by `2 tan α · dr_w`,
so `κ = 1/(2 tan α) = 0.5` at α = 45°. This value reproduces the brief's "≈2,260 rpm"
(UNVERIFIED vs CAD; κ = 1 for a single ramp).

Kinetic energy (flyweight point masses co-rotating with the shaft, moving radially):

    𝒯 = ½ (J_e + n m_w r_w²) ω_p² + ½ (m_ph + n m_w κ²) ẋ_p²

`J_e` is engine + primary rotating inertia excluding flyweights. `m_ph` is the
movable half + spider axial mass. Potential: spring
`U = F_p0 x_p + ½ k_p x_p²`. Lagrange:

    (J_e + n m_w r_w²) ω̇_p + 2 n m_w r_w κ ẋ_p ω_p = T_eng + τ_belt,p − c_θ ω_p        (P1)
    (m_ph + n m_w κ²) ẍ_p − n m_w κ r_w ω_p² = −(F_p0 + k_p x_p) − c_p ẋ_p − F_a,p + F_stop   (P2)

`τ_belt,p` is the belt torque on both primary halves (both turn with the shaft).
`F_a,p ≥ 0` is the belt thrust pushing the halves apart (04 §3).
The `2 n m_w r_w κ ẋ ω` term is the Coriolis torque: the engine loses speed as the
weights fly out. That is real physics, not an error.

**Flyweight axial force** `F_fw = n m_w κ r_w(x_p) ω²` (P2, static part).
**Engagement** (static, `x_p = 0`, no belt load): `F_fw = F_p0`, so

    ω_eng = sqrt(F_p0 / (n m_w κ r_w0))                                    (P3)

**Net axial stiffness** without belt: `k_net = k_p − n m_w κ² ω²`. It is negative for

    ω > ω_neg = sqrt(k_p/(n m_w κ²))                                       (P4)

(the "flyweight force grows with closure" effect; the brief requires it to
arise from the kinematics). (P4) follows from (P2) via `∂F_fw/∂x_p = n m_w κ² ω²`.

## 4. Secondary (helix cam + torsion + linear spring)

Generalized coordinates: shaft angle `θ_s` and opening travel `x_s ∈ [0, s_s]`.
The movable half is constrained by an ideal helix of lead
`p_h = r_c tan φ_c` [m/rad]. Its rotation relative to the shaft is
`θ_rel = σ x_s / p_h`, with `σ = −1`: opening rotates the movable half backward,
so the drive torque pushes it closed (torque feedback).

    𝒯 = ½ J_s ω_s² + ½ J_m (ω_s + σ ẋ_s/p_h)² + ½ m_m ẋ_s²
    U = F_l0 x_s + ½ k_l x_s² + ½ k_t (θ_t0 + x_s/p_h)²

`J_s` is the shaft + fixed half + downstream reflected inertia. `J_m` and `m_m` are the movable half.
`θ_t0 = 30°` is the torsion preload twist. The mass matrix in `(ω_s, ẋ_s)` is
`[[J_s+J_m, σJ_m/p_h],[σJ_m/p_h, m_m + J_m/p_h²]]`. Generalized forces:

    Q_θ = τ_fixed + τ_mov − T_load(ω_s)
    Q_x = σ τ_mov/p_h + F_a,s − (F_l0 + k_l x_s) − k_t(θ_t0 + x_s/p_h)/p_h − c_s ẋ_s + F_stop      (S1)

Friction on the movable flank uses sheave speed `ω_s + σẋ_s/p_h`. On the fixed flank it uses `ω_s`.

**Static secondary closing force** (no inertia, symmetric flanks so
`τ_mov = Tq/2`):

    F_s(x_s, Tq) = F_l0 + k_l x_s + [Tq/2 + k_t(θ_t0 + x_s/p_h)]/p_h     (S2)

**Preload split.** The brief gives a total equivalent preload of 435.63 N. The
torsion part at `x_s = 0` is `k_t θ_t0/p_h = 3.930·0.5236/0.0381 = 54.009 N`.
We assume 435.63 N *includes* it, so `F_l0 = 381.621 N`. This decomposition is an
assumption (A5.4).

## 5. Symbols and values
| Symbol | Meaning | Value | Source |
|---|---|---|---|
| `n`, `m_w` | flyweights, mass each | 3, 0.3137 kg (0.275 + 0.0387) | brief |
| `r_w0` | flyweight radius at x_p=0 | 0.04445 m | brief |
| `α`, `κ` | ramp angle, kinematic ratio | 45°, 0.5 | brief / derived, UNVERIFIED |
| `k_p`, `F_p0` | primary spring rate, preload | 13208.0 N/m, 683.04 N | brief |
| `s_p`, `s_s` | strokes | 25.48 mm, 24.48 mm | brief |
| `k_l`, `F_l0` | secondary linear spring, preload | 3353.66 N/m, 381.621 N | brief / derived |
| `k_t`, `θ_t0` | torsion spring | 3.930 N·m/rad, 30° | brief |
| `r_c`, `φ_c`, `p_h` | cam radius, angle, lead | 0.0381 m, 45°, 0.0381 m/rad | brief / derived |
| `c_p`, `c_s` | axial dampers | 10 N·s/m | brief (assumed), sensitivity run |
| `m_ph`, `m_m`, `J_m` | sheave masses/inertia | PLACEHOLDER 0.8 kg, 0.8 kg, 1.5e-3 kg·m² | — |

## 6. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A5.1 | Straight ramp ⇒ constant κ | ? (real ramp profiles tune shift curve) |
| A5.2 | Frictionless ramps, rollers, helix, splines | ↓ hysteresis; ↑ shift responsiveness |
| A5.3 | Flyweights as point masses (no own rotary inertia) | ↓ primary inertia slightly |
| A5.4 | 435.63 N includes the torsion-spring axial share | ? (if not, ↓ secondary clamp by 54 N) |
| A5.5 | Equal torque split between flanks in (S2) | exact only for symmetric slip; dynamics use actual τ_mov |
| A5.6 | Hard stops modeled as stiff penalty (1e7 N/m, 5 % damping) | ? (stop impact loads depend on it) |

## 7. Implementation contract (`cvtsim/cvt/`)
* `geometry.py`: `belt_length(r1, r2, C)`, `derive_geometry(beta, s_p, s_s, i_low, i_high) -> dict(rp_min, rp_max, rs_min, rs_max, C, L)`,
  `r_star_primary(x_p, g)`, `r_star_secondary(x_s, g)`, `ratio_from_xp(x_p, g)` (belt-length closure, solves `r_s` from G1).
* `primary.py`: `flyweight_axial_force(x_p, omega, p)`, `spring_force(x_p, p)`, `engagement_omega(p)`, `negative_stiffness_omega(p)`,
  `primary_mass(x_p, p) -> (J_rot, m_ax)`, `primary_coriolis_torque(x_p, xdot_p, omega, p)`.
* `secondary.py`: `secondary_static_force(x_s, Tq, p)` (S2), `secondary_mass_matrix(p)`, `secondary_spring_force(x_s, p)` (all non-belt, non-cam-torque terms of Q_x).
* Parameter object `p` is a plain dict loaded from YAML (`cvtsim.io.params.load()`).

## 8. NUMERIC TEST TARGETS (`tests/test_cvt.py`)
| Test | Target | Tol |
|---|---|---|
| T5.1 `derive_geometry` (β=45°) | r_p,min 0.0101964, r_p,max 0.0229364, r_s,max 0.0305891, r_s,min 0.0183491, C 0.251495, L 0.632776 m | rel 1e-5 |
| T5.2 engagement rpm, κ = 0.5 (P3) | 1725.66 rpm | rel 1e-5 |
| T5.3 engagement rpm, κ = 1 | 1220.23 rpm | rel 1e-5 |
| T5.4 negative-stiffness onset κ=0.5 (P4) | 2262.57 rpm (brief: "about 2,260") | rel 1e-5 |
| T5.5 `k_net = d(F_spring − F_fw)/dx_p` (finite difference) is < 0 at 2300 rpm and > 0 at 2200 rpm | – | – |
| T5.6 secondary torsion axial preload | 54.009 N | rel 1e-4 |
| T5.7 `secondary_static_force(0, 20)` | 435.63 + 10/0.0381 = 698.097 N | rel 1e-5 |
| T5.8 `secondary_static_force(0.01, 0)` | 381.621 + 33.5366 + 3.930·(0.5236+0.26247)/0.0381 = 496.235 N | rel 1e-4 |
| T5.9 energy: flyweight-primary free system (P1–P2 with T_eng = 0, no belt, no damping, no stops, start x_p inside the stroke) integrated with rk4 1 s: `𝒯 + U` conserved | rel 1e-6 |
| T5.10 `ratio_from_xp(0)` = 3.0, `ratio_from_xp(s_p)` = 0.8 | rel 1e-4 |
