# 01 — Rigid-body dynamics and joints (`cvtsim.core`)

## 1. Sources
* J. Baumgarte, "Stabilization of constraints and integrals of motion in
  dynamical systems", *Comput. Methods Appl. Mech. Eng.* 1 (1972) 1–16,
  doi:10.1016/0045-7825(72)90018-7. Verified (bibliographic record). Used for
  the stabilized constraint equation (§4).
* E. J. Haug, *Computer-Aided Kinematics and Dynamics of Mechanical Systems,
  Vol. I*, Allyn & Bacon, 1989. Used for the constrained Newton–Euler form
  and the joint primitives (dot-1, dot-2, spherical). The bibliographic record
  was not re-checked online, and the chapter numbers are **UNVERIFIED**.
* Pendulum period with complete elliptic integral `K`: standard. E.g. L. D.
  Landau & E. M. Lifshitz, *Mechanics*, 3rd ed., §11, problem 1. Section
  number **UNVERIFIED**. The target is computed with `scipy.special.ellipk`.
* Torque-free symmetric top: Landau & Lifshitz, *Mechanics*, §33 (Euler's
  equations), §35 (free symmetric top). Section numbers **UNVERIFIED**. The
  closed form is derived below, so the tests do not depend on the citation.

## 2. Symbols

| Symbol | Meaning | Unit |
|---|---|---|
| `r_i` | body *i* center-of-mass position, world frame | m |
| `q_i` | unit quaternion (w,x,y,z), body→world rotation `R_i = R(q_i)` | – |
| `v_i` | CoM velocity, world | m/s |
| `w_i` | angular velocity, **world** frame | rad/s |
| `m_i`, `Jb_i` | mass; inertia tensor about CoM in body frame | kg, kg·m² |
| `Jw_i` | `R_i Jb_i R_iᵀ`, world-frame inertia | kg·m² |
| `s_i` | joint attachment point, body frame (from CoM) | m |
| `rho_i` | `R_i s_i`, world-frame lever arm | m |
| `Phi` | constraint vector (positions) | m or – |
| `Jc` | constraint Jacobian, `dPhi/dt = Jc u` | – |
| `lambda` | Lagrange multipliers (constraint forces) | N, N·m |
| `u` | stacked `[v_1, w_1, v_2, w_2, …]` | – |

## 3. Equations of motion (constrained Newton–Euler)

For each body, with `[a]×` the skew matrix such that `[a]× b = a × b`:

    m_i dv_i/dt            = f_i                       + (Jc_iᵀ λ)_trans
    Jw_i dw_i/dt           = tau_i − w_i × (Jw_i w_i)  + (Jc_iᵀ λ)_rot          (1)
    dr_i/dt = v_i ,   dq_i/dt = ½ (0, w_i) ⊗ q_i                                (2)

Here `f_i` and `tau_i` are the applied force and the torque about the CoM, both in the world
frame. Eq. (1) is the Newton–Euler equation with the gyroscopic term
`−w × Jw w`. In the world frame `d/dt(Jw w) = Jw dw/dt + w × (Jw w)`.
Stacked: `M u̇ = F(q,u) + Jcᵀ λ`, with `M = blockdiag(m_i I3, Jw_i)` and `F = [f_i, tau_i − w_i×Jw_i w_i]`.

**Index reduction with Baumgarte stabilization (Baumgarte 1972).** Replace
`Phi = 0` by `Phi̇̇ + 2 a Phi̇ + b² Phi = 0`, i.e. `Jc u̇ = gamma − 2 a Jc u − b² Phi`
with `gamma ≡ −(J̇c u)` (quadratic velocity terms, §5). Solve the KKT system

    [ M   −Jcᵀ ] [ u̇ ]   [ F                          ]
    [ Jc    0  ] [ λ ] = [ gamma − 2a Jc u − b² Phi    ]                       (3)

with a dense solve (systems here have < 100 equations). Default `a = b = 0`
(pure index-1). Stabilization is optional.

**Position projection (optional, recommended).** After each step, project
positions back onto the manifold with Gauss–Newton steps
`dq = −Jcᵀ (Jc Jcᵀ)⁻¹ Phi`. Translations are applied directly. Rotations
are applied as a small-rotation update `q ← exp(½ dθ) ⊗ q`. Iterate until `|Phi|∞ < 1e-12` (max 5
iterations). Then project velocities: `u ← u − M⁻¹Jcᵀ(Jc M⁻¹ Jcᵀ)⁻¹ Jc u`.

## 4. Joint primitives

Notation: body *i* (may be `ground`), body *j*. `d = r_j + rho_j − r_i − rho_i`.
`a = R_i a_i`, with `a_i` a body-*i* unit vector. `b = R_j b_j`.

| Primitive | `Phi` | Jacobian row blocks `[v_i, w_i, v_j, w_j]` | `gamma` (= −J̇c u) |
|---|---|---|---|
| spherical (3 rows) | `d` (sign: `r_i+rho_i−r_j−rho_j`) | `[I, −[rho_i]×, −I, [rho_j]×]` | `−(w_i×(w_i×rho_i) − w_j×(w_j×rho_j))` |
| dot-1 (1 row) | `a·b` | `[0, (a×b)ᵀ, 0, −(a×b)ᵀ]` | `−[(w_i×(w_i×a))·b + 2(w_i×a)·(w_j×b) + a·(w_j×(w_j×b))]` |
| dot-2 (1 row) | `a·d` | `[−aᵀ, (a×d − rho_i×a)ᵀ, aᵀ, (rho_j×a)ᵀ]` | `−[a·(w_j×(w_j×rho_j) − w_i×(w_i×rho_i)) + 2(w_i×a)·ḋ + d·(w_i×(w_i×a))]` with `ḋ = v_j+w_j×rho_j−v_i−w_i×rho_i` |

Derivations. Spherical: `d/dt(R s) = w × (R s)`. Differentiating again gives
`ẇ×rho + w×(w×rho)`. The second term is the quadratic part. Dot-1:
`d/dt(a·b) = (w_i×a)·b + a·(w_j×b) = w_i·(a×b) − w_j·(a×b)`. Dot-2:
`d/dt(a·d) = (w_i×a)·d + a·ḋ`. Use `(w×a)·d = w·(a×d)` and
`a·(w×rho) = w·(rho×a)`. The second derivatives keep the quadratic terms
`2 ȧ·ḋ`, `a·d̈_quad` and `d·ä_quad`.

**Joints (rows):**
* **Revolute** about body-*i* axis `a_i` (with body-*j* axis `b_j` aligned):
  spherical (3) plus dot-1(`u_i`, `b_j`) and dot-1(`w_i`, `b_j`), where
  `u_i`, `w_i` ⟂ `a_i`. 5 rows.
* **Prismatic** along `a_i`: dot-1 on 3 axis pairs to lock relative
  rotation, and dot-2 with two vectors `u_i`, `w_i` ⟂ `a_i` to lock transverse
  translation. 5 rows. The three rotation rows use pairs (`u_i`,`w_j`),
  (`u_i`,`a_j`), (`w_i`,`a_j`) with matched initial orthogonal body triads.
* **Screw / helix** with lead `p` [m/rad] along `a_i`: cylindrical (dot-1 ×2
  for axis alignment, dot-2 ×2 for transverse offset) plus the screw row

      Phi_s = a·d − p·θ_rel − h0 ,   θ_rel = unwrapped relative rotation about a
      Jc_s  = Jc_dot2(a, d) + p·[0, aᵀ, 0, −aᵀ]
      gamma_s = gamma_dot2(a, d) + p·(w_i×a)·(w_j − w_i)

  `θ_rel` comes from `atan2` of reference vectors ⟂ `a`, unwrapped
  continuously by stored state. The term `p(w_i×a)·(w_j−w_i)` follows from
  `θ̈_rel = a·(ẇ_j−ẇ_i) + (w_i×a)·(w_j−w_i)`. Sign: `Phi_s` increases with
  `d` along `a`. Positive `θ_rel` (j turns +about `a` relative to i)
  advances the body by `+p θ`.
* **Ground** is a body with infinite mass. Its blocks are omitted.

## 5. Assumptions

| # | Assumption | Bias |
|---|---|---|
| A1.1 | Rigid bodies (no deformation in core) | ↑ peak loads (no compliance to absorb impacts) |
| A1.2 | Ideal frictionless joints | ↓ dissipation, ↑ efficiency |
| A1.3 | Baumgarte/projection modifies dynamics by O(h) off-manifold correction | ? (small; checked by drift test) |

## 6. Implementation contract
`cvtsim/core/rigid.py`: `quat_mul`, `quat_to_R`, `quat_exp(rotvec)`, `RigidBody(m, Jb, r, q, v, w)`.
`cvtsim/core/joints.py`: classes `Spherical`, `Revolute`, `Prismatic`, `Screw`. Each has
`phi(sys)`, `jac(sys)`, `gamma(sys)`, `nrows`.
`cvtsim/core/mbs.py`: `System(bodies, joints, gravity=(0,0,-9.81), baumgarte=(a,b), project=bool)`, with
`accelerations()`, `step(h, method='symplectic_euler'|'rk4')`, `energy()` (kinetic + gravity potential),
and `constraint_violation()` (max |Phi|).
`symplectic_euler`: `u_{n+1} = u_n + h u̇(q_n,u_n)`, `r_{n+1} = r_n + h v_{n+1}`,
`q_{n+1} = quat_exp(h w_{n+1}) ⊗ q_n` (normalized).

## 7. NUMERIC TEST TARGETS (`tests/test_core.py`)

Pendulum: point-like body `m = 1 kg`, `Jb = 1e-6·I` (≈ point mass), on a revolute
joint to ground at distance `L = 1 m`, `g = 9.81`. Release from rest at θ0.

| Test | Target | Tolerance |
|---|---|---|
| T1.1 period θ0=10° | 2.00989 s | rel 1e-3 (h = 1e-4, rk4 or s-Euler with projection) |
| T1.2 period θ0=60° | 2.15287 s | rel 1e-3 |
| T1.3 period θ0=120° | 2.75409 s | rel 1e-3 |
| T1.4 energy, θ0=60°, 10 s, h=1e-3, rk4 + projection | `max|E−E0|/|E0|` | < 1e-6 |
| T1.4b projected symplectic Euler: energy error at 2 s, h = 1e-3 vs 5e-4 | ratio 2 (first-order drift; projection breaks symplecticity, so the error is not bounded) | ratio in [1.7, 2.3] |
| T1.5 constraint drift, θ0=60°, 10 s, h=1e-3, Baumgarte a=b=1/(10h)... see note | `max|Phi|` | ≤ 10·g·h² = 9.81e-5 m |
| T1.6 drift with projection | `max|Phi|` | ≤ 1e-10 m |
| T1.7 drift **without** stabilization or projection is larger than with projection | qualitative ordering | strict |
| T1.8 free symmetric top `I1=I2=1, I3=2`, body `w=(1,0,3)` (rk4, h=1e-3, 10 s): body-frame precession rate of (w1,w2) | `(I3−I1)/I1·w3 = 3.000 rad/s` | rel 1e-3 |
| T1.9 same: rate at which symmetry axis precesses about `L` in space | `|L|/I1 = 6.08276 rad/s` | rel 1e-3 |
| T1.10 same: `|L_world|` conserved | rel drift | < 1e-6 (rk4) |
| T1.11 screw joint, lead p=0.01 m/rad, body m=2 kg, axial inertia I=0.01 kg·m², axial force F=10 N, no gravity: axial accel | `F/(m + I/p²) = 10/(2+100) = 0.0980392 m/s²` | rel 1e-6 |
| T1.12 screw: `|z − p θ|` over 1 s | ≤ 1e-9 m (projection on) | – |
| T1.13 prismatic block on 30° incline, frictionless: accel along incline | `g sin30° = 4.905 m/s²` | rel 1e-9 |

Note for T1.5: use Baumgarte `a = b = 1/(10h)`. The bound comes from a local truncation error
of `½ g h²` per step, which the stabilizer damps within ~10 steps.
Pendulum period measurement: time between successive same-direction zero
crossings of θ, linearly interpolated, averaged over ≥ 3 periods.
