# 02 — Time integration (`cvtsim.core.integrators`)

## 1. Sources
* J. Chung & G. M. Hulbert, "A time integration algorithm for structural
  dynamics with improved numerical dissipation: the generalized-α method",
  *J. Appl. Mech.* 60(2) (1993) 371–375, doi:10.1115/1.2900803. Verified.
  The parameter formulas in §3 are from that paper (its Eqs. for
  `α_m, α_f, γ, β` in terms of `ρ∞`). Equation numbers are **UNVERIFIED**.
* E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd
  ed., Springer, 2006, doi:10.1007/3-540-30666-8. Verified (record). Used for
  Störmer–Verlet symplecticity and bounded energy error (Ch. I and Ch. IX,
  chapter numbers **UNVERIFIED**). Also Hairer, Lubich & Wanner, *Acta
  Numerica* 12 (2003) 399–450 (verified record).

## 2. Problem class
Second-order system `M ẍ + f_int(x, ẋ, t) = f_ext(t)`, `x ∈ R^n`, constant
or state-dependent `M`. Tangents `K = ∂f_int/∂x`, `Cd = ∂f_int/∂ẋ`. If
these are not supplied, they are obtained by forward finite differences.

## 3. Generalized-α (Chung & Hulbert 1993)
Given the user parameter `ρ∞ ∈ [0, 1]` (high-frequency spectral radius):

    α_m = (2ρ∞ − 1)/(ρ∞ + 1),  α_f = ρ∞/(ρ∞ + 1),
    γ   = ½ − α_m + α_f,        β   = ¼ (1 − α_m + α_f)²

Newmark updates:

    x_{n+1} = x_n + h ẋ_n + h² [(½ − β) a_n + β a_{n+1}]
    ẋ_{n+1} = ẋ_n + h [(1 − γ) a_n + γ a_{n+1}]

Balance at the intermediate point:

    M a_{n+1−α_m} + f_int(x_{n+1−α_f}, ẋ_{n+1−α_f}) = f_ext(t_{n+1−α_f})
    with (·)_{n+1−α} = (1 − α)(·)_{n+1} + α(·)_n

Unknown `a_{n+1}`. Solve with Newton. The effective Jacobian is
`(1−α_m) M + (1−α_f)(γ h Cd + β h² K)`. Converge to
`|residual| ≤ 1e-10·max(1,|f_ext|)`, at most 20 iterations, and raise an error
if it does not converge. Initial `a_0` is solved from the equation of motion at `t_0`.

Properties (paper): second-order accurate, unconditionally stable for linear
problems, and the spectral radius → `ρ∞` as `ωh → ∞`. `ρ∞ = 1` gives
`α_m = α_f = ½`, `γ = ½`, `β = ¼`: no numerical dissipation.

## 4. Störmer–Verlet / symplectic Euler (explicit)
Velocity Verlet for `ẍ = a(x, ẋ, t)`:

    v_{½} = v_n + ½ h a(x_n, v_n, t_n)
    x_{n+1} = x_n + h v_{½}
    v_{n+1} = v_{½} + ½ h a(x_{n+1}, v_{½}, t_{n+1})     (velocity-dependent: lagged)

For conservative `a(x)` it is symplectic and second order. The energy error is
bounded with no secular drift (Hairer et al.). Linear stability on
`ẍ = −ω² x` requires `ω h < 2`.

Symplectic (semi-implicit) Euler: `v_{n+1} = v_n + h a(x_n, v_n, t_n)`,
`x_{n+1} = x_n + h v_{n+1}`. It is first-order, symplectic, and stable for `ωh < 2`.
The launch simulation uses this scheme (cheapest per step). Accuracy is verified there
by halving `h` (Phase 2).

## 5. Symbols
| Symbol | Meaning | Unit |
|---|---|---|
| `h` | time step | s |
| `ρ∞` | spectral radius at infinite frequency | – |
| `α_m, α_f, γ, β` | gen-α parameters | – |
| `ω` | natural circular frequency | rad/s |

## 6. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A2.1 | Lagged velocity in Verlet/symplectic Euler for damping/friction forces | ? (O(h); checked by `h`-halving) |
| A2.2 | `ρ∞ < 1` removes high-frequency energy numerically | ↓ high-frequency load content |

## 7. Implementation contract
`cvtsim/core/integrators.py`:
* `genalpha_params(rho_inf) -> (am, af, gamma, beta)`
* `GeneralizedAlpha(M, f_int, f_ext, rho_inf, K=None, Cd=None)` with `.run(x0, v0, t0, h, nsteps) -> (t, X, V, A)`
* `genalpha_amplification(omega_h, rho_inf, xi=0) -> 3×3 matrix` for the linear SDOF.
  State `(x, h·v, h²·a)`, derived from the update equations.
* `velocity_verlet(accel, x0, v0, t0, h, nsteps)`, `symplectic_euler(accel, x0, v0, t0, h, nsteps)` → `(t, X, V)`

## 8. NUMERIC TEST TARGETS (`tests/test_integrators.py`)
SDOF `ẍ + ω² x = 0`, `ω = 2π`, `x0 = 1`, `v0 = 0`. Exact: `x = cos ωt`. Errors are
`max_t |x(t) − cos ωt|` over `t ∈ [0, 1]`. Evaluating at `t = 1` alone (one period) cancels the
leading phase error and shows spurious superconvergence.

| Test | Target | Tolerance |
|---|---|---|
| T2.1 `genalpha_params(0.8)` | `(0.333333, 0.444444, 0.611111, 0.308642)` | abs 1e-6 |
| T2.2 gen-α (ρ∞=0.8) max error on [0,1] for h = 1/100, 1/200, 1/400: observed order `log2(e_h/e_{h/2})` | 2.0 | ±0.15 |
| T2.3 spectral radius of amplification matrix at `ωh = 1e5` (deviation from ρ∞ decays ≈ 22/(ωh)) for ρ∞ ∈ {0, 0.5, 0.8, 1} | = ρ∞ | abs 2e-3 |
| T2.4 spectral radius ≤ 1 + 1e-12 for `ωh ∈ logspace(-3, 4, 200)`, ρ∞ ∈ {0, 0.5, 1} | unconditional stability | – |
| T2.5 ρ∞=1: spectral radius = 1 for all `ωh` above | no dissipation | abs 1e-10 |
| T2.6 Verlet order (same SDOF, h = 1/100…1/400) | 2.0 | ±0.15 |
| T2.7 Verlet stability: 10⁴ steps at `ωh = 1.99` → `max|x| ≤ 20`; at `ωh = 2.01` → `max|x| > 1e3` | threshold `ωh = 2` | – |
| T2.8 Verlet energy, pendulum `θ̈ = −(g/L) sin θ`, θ0 = 1 rad, h = 0.01, 10⁵ steps | `max|E−E0|/E0 < 1e-3` **and** max over last 10 % ≤ 1.5 × max over first 10 % (no secular drift) | – |
| T2.9 gen-α (ρ∞=0.9) nonlinear: Duffing `ẍ + x + x³ = 0`, x0=1, h=1e-3, compare to `solve_ivp` (RK45, rtol 1e-10) at t=5 | abs err < 1e-4 | – |
| T2.10 symplectic Euler order on SDOF | 1.0 | ±0.15 |
