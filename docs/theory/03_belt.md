# 03 — Segmented belt (`cvtsim.belt`)

## 1. Sources
* Euler–Eytelwein (capstan) relation, extended to a V-groove: K. L. Johnson,
  *Contact Mechanics*, Cambridge Univ. Press, 1985. The belt-drive section's
  number is **UNVERIFIED**. The V-groove form is derived in 04_contact.md §3, so
  the test target does not rely on the citation.
* Beaded-string (discrete chain) frequencies: Lord Rayleigh, *The Theory of
  Sound*, Vol. 1, §120 (Lagrange's loaded string). Section number
  **UNVERIFIED**. The formula is re-derived in §4 below.
* Belt rubber shear modulus G = 30.15 MPa (hot) / 46.71 MPa (cold), small
  strain: Messick (2018), thesis, as given in the brief. The bibliographic
  record could **not** be verified online, so it is **UNVERIFIED**.

## 2. Model
Closed (or open, for tests) chain of `N` nodes `x_k ∈ R²` in the belt plane
(X–Y), lumped mass `m_k = ρA·L0`, with `L0 = L/N` (closed) or `L/(N−1)` (open).
Element `e = (k, k+1)` (closed: indices mod N):

* **Axial (Kelvin–Voigt):** `l = |x_{k+1}−x_k|`, `t̂ = (x_{k+1}−x_k)/l`,
  strain `ε = (l − L0)/L0`, strain rate `ε̇ = t̂·(v_{k+1}−v_k)/L0`,
  tension `T = EA (ε + η ε̇)`. Force on node k+1 is `−T t̂`, on node k it is `+T t̂`.
  The element carries compression too (A3.2).
* **Bending (discrete rotational springs at nodes):** turning angle at node k
  between `e = x_k − x_{k−1}` and `f = x_{k+1} − x_k`:
  `φ_k = atan2(e×f, e·f)`, where `e×f = e_x f_y − e_y f_x`. Energy
  `U_b = Σ ½ (EI/L0) φ_k²`. Gradient (with `ẑ×u = (−u_y, u_x)`):

      ∂φ_k/∂x_{k+1} = ẑ×f/|f|²,   ∂φ_k/∂x_{k−1} = ẑ×e/|e|²,
      ∂φ_k/∂x_k     = −(ẑ×f/|f|² + ẑ×e/|e|²)

  Force `F = −(EI/L0) φ_k ∂φ_k/∂x` on the three nodes. Derivation:
  `angle(u) = atan2(u_y,u_x)` has gradient `ẑ×u/|u|²`, and `φ = angle(f) − angle(e)`.
  An optional bending damper `M_d = (c_b/L0) φ̇_k` uses `φ̇ = ∇φ·v`.
  Open chains have no bending spring at end nodes.
* **Transverse compression stiffness** is used by contact (04). For a flank
  contact patch of height `h_c` (measured along the flank) per unit belt length, with lateral
  half-width `b/2` and rubber modulus `E_t = 3G` (incompressible, small
  strain), the normal stiffness per flank per node is

      k_n = E_t · h_c · L0 / (b/2)         [N/m]                         (B1)

  This is a column-in-compression estimate. Bulging and cords are ignored (A3.4).

Strain energy `U = Σ_e ½ EA L0 ε² + Σ_k ½ (EI/L0) φ_k²`. Dissipated power
`P_d = Σ_e EA η L0 ε̇² (+ bending)`. Both must be returned for energy checks.

## 3. Symbols
| Symbol | Meaning | Unit | DR26 value / status |
|---|---|---|---|
| `EA` | axial stiffness | N | 243.5e3 (brief, UNVERIFIED) |
| `ρ` | density | kg/m³ | 1150 (brief, assumed) |
| `A` | cross-section area | m² | PLACEHOLDER 3.0e-4 |
| `EI` | bending stiffness | N·m² | PLACEHOLDER 0.05 |
| `η` | axial Kelvin–Voigt time constant | s | PLACEHOLDER 2e-6 (≈0.4 critical on the highest segment mode, explicit-step stable; uncalibrated) |
| `G` | rubber shear modulus | Pa | 30.15e6 hot / 46.71e6 cold (Messick 2018, UNVERIFIED) |
| `b` | belt top width | m | PLACEHOLDER 0.030 |
| `h_c` | flank contact height | m | PLACEHOLDER 0.012 |
| `L` | belt pitch length | m | 0.632776 derived (05 §2) |

## 4. Derivations of test targets
* **Static stretch:** open chain, node 0 fixed, axial load `P` at the free end.
  Each element carries `T = P`, so the total elongation is `δ = P L/EA` (exact for
  the discrete model with linear `EA`).
* **Transverse vibration of a tensioned span** (EI = 0): N segments between
  pinned ends, pre-tension `T0` (from pre-stretch `ε0 = T0/EA`), node mass
  `m = ρA L0`. Linearize about the straight line. The transverse stiffness of
  each element is `T0/L0` (geometric stiffness). The mode-1 eigenvalue of the beaded string is

      ω1 = 2 sqrt(T0/(m l)) sin(π/(2N))                                     (B2)

  where `l = ℓ/N` is the stretched spacing and `m = ρA l/(1+ε0)`. The mass is fixed
  by the unstretched length, so (B2) equals the unstretched-spacing value times
  `sqrt(1+ε0)`.

  As N→∞ this tends to `(π/ℓ) sqrt(T0/ρA)` (continuum string, `f = (1/2ℓ)√(T0/ρA)`).
  With bending, the continuum gives `ω² = (T0/ρA)k² + (EI/ρA)k⁴` with `k = π/ℓ`.
* **Capstan at slip:** see 04_contact.md §3. Target `T1/T2 = exp(μθ/sinβ)`.

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A3.1 | Planar belt; no axial (Z) motion, twist or misalignment | ↓ belt side loads; ? on tension |
| A3.2 | Linear axial law in tension *and* compression (no slack-span buckling) | ↑ compressive forces in slack span (small; spans are in tension under drive) |
| A3.3 | Kelvin–Voigt damping, uncalibrated `η` | ? (sensitivity run in Phase 2) |
| A3.4 | `k_n` from a column estimate `E_t = 3G` (B1) | ? — probably ↑ stiffness (ignores bulging into the groove) |
| A3.5 | Lumped mass; rotary inertia of segments neglected | ↓ (negligible at L0 ≈ 4 mm) |
| A3.6 | Pitch line = node line; no neutral-axis shift in bending | ? |

## 6. Implementation contract (`cvtsim/belt/belt.py`)
All kernels are plain-NumPy, array-in/array-out functions with no Python
objects. Decorate them with `cvtsim.core.jit.maybe_njit`, which is `numba.njit(cache=True)` if
numba is importable and the identity otherwise. That lets Phase 2 call them from a jitted loop.
* `belt_forces(x, v, L0, EA, eta, EI, cb, closed) -> (f (N,2), U_strain, P_diss)`
* `belt_tensions(x, v, L0, EA, eta, closed) -> T (Ne,)` (element tensions, N)
* `init_belt_on_path(rp, rs, C, N) -> x (N,2)`: N nodes equally spaced by arc length on the
  open-belt path around circles `(0,0,rp)` and `(C,0,rs)`. Counter-clockwise order (the belt moves
  `+Z` sense). Returns positions and `L0`.
* `belt_length(r1, r2, C)` (exact open-belt formula, 05 §2).
* `node_mass(rhoA, L0)`; `flank_stiffness(G, h_c, L0, b)` per (B1).

## 7. NUMERIC TEST TARGETS (`tests/test_belt.py`)
| Test | Setup | Target | Tol |
|---|---|---|---|
| T3.1 static stretch | open chain, L = 1 m, N = 50, EA = 243.5e3, P = 500 N, solve statics | δ = P L/EA = 2.0533881e-3 m | rel 1e-9 |
| T3.2 string, discrete | span ℓ=0.2 m, N=20 segments, T0=100 N, ρA=0.345 kg/m, EI=0: lowest transverse eigenfrequency of the linearized system (numerical Jacobian of `belt_forces`, open chain with ends fixed) | 42.5278 Hz (B2, ε0 = T0/EA included) | rel 1e-4 |
| T3.3 string, continuum | same, N = 20 vs continuum 42.5628 Hz | rel 2e-3 | – |
| T3.4 string with bending | N = 100, EI = 0.05: compare to `sqrt((T0/ρA)k²+(EI/ρA)k⁴)/(2π)` | rel 2e-2 |
| T3.5 time-domain string | integrate T3.2 with velocity Verlet (h = 1e-6, small amplitude 1e-4 m); measured period from zero crossings | 1/42.5278 s | rel 1e-3 |
| T3.6 energy | closed belt, EI>0, η=0, random small initial velocities, Verlet 1e4 steps | `max|ΔE|/E0 < 1e-3` | – |
| T3.7 rigid rotation | closed belt on circle radius R rotating rigidly at ω, `η > 0`: damping force is zero | `|f_damp| < 1e-9` | – |
| T3.8 length | `init_belt_on_path(0.0101964, 0.0305891, 0.251495, 200)` total polyline length vs `belt_length` (chord error) | rel 1e-3 | – |
