# 04 — V-belt contact with analytic cone sheaves (`cvtsim.contact`)

## 1. Sources
* G. Gerbert, *Traction Belt Mechanics: Flat Belts, V-belts, V-ribbed Belts*,
  Chalmers University of Technology, Göteborg, 1999 (WorldCat OCLC 247549000;
  existence verified). We could not open the text, so **every chapter/section
  reference to Gerbert is UNVERIFIED**. Gerbert is used for (a) the
  wedge-action / effective friction `μ/sin β` in a V-groove and (b) the
  qualitative facts that, under torque, a V-belt shows *slip that grows with
  load* and *radial sliding* in the groove. Our quantitative tests use only
  closed forms derived here.
* Penalty regularization of unilateral contact and Coulomb friction:
  P. Wriggers, *Computational Contact Mechanics*, 2nd ed., Springer, 2006,
  Ch. 2 & 5 (chapter numbers **UNVERIFIED**).

## 2. Geometry of an analytic cone sheave (no mesh contact search)
A sheave pair (fixed half + movable half) has axis `c = (c_x, c_y)` in the
belt plane. Both flanks make the half-angle `β` with the radial plane, so the
groove width grows outward as `w(r) = b + 2 (r − r*) tan β`. Here `r*` is the
**free radius**, the radius at which a belt of top width `b` touches both flanks
with zero compression. Moving one half axially by `Δx` changes the width
by `Δx`, so

    r*(x) = r*_0 ± x / (2 tan β)                                           (C1)

(+ for primary closing travel `x_p`, − for secondary opening travel `x_s`).

For a node at in-plane distance `r = |x_k − c|`, with unit `ê_r = (x_k − c)/r`:

    penetration normal to each flank:  δ = (r* − r) sin β   if  r < r*  and  r < r_outer ,  else 0     (C2)

This is exact for a rigid wedge. Interference across the groove is `2(r*−r)tanβ`.
Each flank's share is `(r*−r) tan β` axially, and its component along the flank normal is
`(r*−r) tan β cos β = (r*−r) sin β`.

## 3. Forces
Per node and per flank (two flanks, `k ∈ {fixed, movable}`):

* **Normal (penalty + damping, compressive only):**
  `N = max(0, k_n δ + c_n δ̇)`, with `δ̇ = (ṙ* − ṙ) sin β` and `ṙ = ê_r·v_k`.
* **Resultant of the two normals on the node:** radial, outward,
  `F_r = 2 N sin β ê_r`. Axial components `±N cos β` cancel on the node and
  push the two halves apart. **Axial thrust on the movable half:** `F_a = Σ_nodes N cos β`.
* **Friction (regularized Coulomb, in-plane):** flank surface velocity
  `u_k = ω_k ẑ × (x − c)`. Slip `s_k = v − u_k` has both tangential and radial
  components (radial sliding, Gerbert):

      f_k = −μ N s_k / sqrt(|s_k|² + v_ε²)                                  (C3)

  The force on the sheave half is `−f_k`. Its torque about `+Z` is
  `τ_k = ((x−c) × (−f_k))_z`. Dissipated power `P_k = −f_k·s_k ≥ 0`.
  Friction acts only in the X–Y plane. Its axial component is neglected (A4.3).

**Energy consistency.** Penalty energy `U = Σ_nodes 2·½ k_n δ²`. Then
`−∂U/∂x_k = 2 k_n δ sin β ê_r`, which matches `F_r`. Also `∂U/∂r* = Σ 2 N sin β`
(with `c_n = 0`). With (C1), `|∂U/∂x| = Σ 2 N sin β/(2 tan β) = Σ N cos β = F_a`.
Sheave axial work and belt potential energy are therefore conjugate, and the
energy balance in Phase 2 closes.

### 3.1 V-groove capstan (derivation of test targets)
Belt element `ds = r dθ` on a groove at constant radius, sliding everywhere
(fully developed slip). Let `n` be the normal force per unit length per flank.
* Radial equilibrium (centrifugal neglected): `T dθ = 2 n sin β ds`.
* Tangential equilibrium: `dT = 2 μ n ds`.

Divide: `dT/T = (μ/sin β) dθ`, so

    T1/T2 = exp(μ θ / sin β)                                               (C4)

Axial thrust per flank: `F_a = ∫ n cos β ds = (cos β/(2 sin β)) ∫ T dθ`
`= (cos β/(2 sin β)) (T1 − T2) sin β/μ`. Hence

    F_a = cos β (T1 − T2) / (2 μ) = Tq cos β / (2 μ R)                     (C5)

with `Tq = (T1 − T2) R`. This is the brief's relation `F_a = T cos β/(2 μ R)`,
exact for full-wrap slip, independent of wrap angle.

With `μ = 0.3`, `β = 45°`, `θ = π`: `T1/T2 = 3.79189`.

## 4. Symbols
| Symbol | Meaning | Unit | DR26 |
|---|---|---|---|
| `β` | flank half-angle from radial plane | rad | 45° (brief; CAD shows ≈43.3° on some faces — DISCREPANCY flagged) |
| `μ` | belt–sheave friction | – | 0.3 (brief, assumed) |
| `k_n` | flank normal stiffness per node | N/m | from (B1) |
| `c_n` | flank normal damping per node | N·s/m | PLACEHOLDER: `2·0.05·sqrt(k_n m_node)` (5 % critical) |
| `v_ε` | friction regularization velocity | m/s | PLACEHOLDER 0.05 (convergence-checked) |
| `r*` | free radius | m | state (C1) |
| `r_outer` | sheave outer radius | m | PLACEHOLDER: r*_max + 5 mm |

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A4.1 | Rigid, perfectly conical flanks (no sheave deflection) | ↑ local contact force peaks; ? on ratio |
| A4.2 | Regularized Coulomb, constant μ (no Stribeck/rate effect) | ? — creep replaces stick (↑ slip at low load) |
| A4.3 | Axial friction component neglected | ↑ effective radial mobility of the belt |
| A4.4 | Penalty compliance acts as belt transverse compliance | ? — converges to rigid wedge as `k_n → ∞` |
| A4.5 | Contact evaluated at nodes only | ↓ resolution of contact-pressure peaks at entry/exit |

## 6. Implementation contract (`cvtsim/contact/sheave.py`)
`sheave_contact(x, v, cx, cy, r_star, r_star_dot, beta, k_n, c_n, mu, v_eps, r_outer, w_fixed, w_mov)`
→ `(f (N,2), Fa, tq_fixed, tq_mov, P_diss, U_pen, delta (N,), N_flank (N,))`.
This is a pure-NumPy kernel, wrapped with `maybe_njit`.
Also provide `capstan_ratio(mu, theta, beta)` and `axial_force_incipient(Tq, R, beta, mu)` (C4, C5).

## 7. NUMERIC TEST TARGETS (`tests/test_contact.py`)
**Capstan rig.** Pulley at origin, `r* = 0.05 m`, `β = 45°`, `μ = 0.3`, `v_ε = 1e-3`,
pulley spinning at `ω = +10 rad/s` (sliding everywhere, slip speed ≈ 0.5 m/s ≫ v_ε).
Open belt chain, `EI = 0`, `EA = 2.435e5`, wrapped over the upper half (θ = π).
The straight legs extend 0.1 m downward. End A is pinned. End B carries a constant
force `T2 = 100 N`, pulling downward. The pulley turns `+Z` and drags the belt toward A, so
the anchored end A is the tight side: `T1/T2 = e^{μπ/sinβ}`. Run with node
damping until static (≤ 1e-6 m/s nodal speed), or solve statics with Newton.

| Test | Target | Tol |
|---|---|---|
| T4.1 `capstan_ratio(0.3, π, π/4)` | 3.79189 | 1e-5 |
| T4.2 rig: anchor force / T2, `k_n = 1e7`, ≥ 60 nodes in wrap | 3.79189 | rel 2e-2 |
| T4.3 rig: axial thrust `F_a` vs (C5) with measured T1, T2 | `cos β (T1−T2)/(2μ)` | rel 2e-2 |
| T4.4 `axial_force_incipient(10, 0.05, π/4, 0.3)` | 235.702 N | 1e-3 |
| T4.5 penalty convergence: `k_n ∈ {1e5, 1e6, 1e7}`: max penetration scales ≈ 1/k_n (ratio of successive maxima in [5, 20]); error of T4.2 non-increasing | – | – |
| T4.6 energy consistency (μ=0, c_n=0): central finite difference of `U_pen` w.r.t. node positions equals `−f`, and `(∂U_pen/∂r*)/(2 tan β) = F_a` | rel 1e-6 |
| T4.7 friction dissipates: `P_diss ≥ 0` for random states | – | – |
| T4.8 no contact for `r > r*` or `r > r_outer`: forces exactly 0 | – | – |
| T4.9 slip trend (Gerbert, qualitative): two-pulley drive (equal r* = 0.05, C = 0.2, N=120 nodes, fixed `r*` both, belt rest length 0.2 % shorter than the path → pretension ≈ 490 N), driver at constant ω = 100 rad/s, driven pulley (inertia 1e-3 kg·m²) loaded with torque `Tq ∈ {2, 5, 10, 15}` N·m: relative speed loss `1 − ω_out/ω_in` at steady state is **strictly increasing** in `Tq` | marked `@pytest.mark.slow` | – |
