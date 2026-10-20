# 08 — Linear static FE, 10-node tetrahedra (`cvtsim.fe.tet10`)

## 1. Sources
* O. C. Zienkiewicz, R. L. Taylor & J. Z. Zhu, *The Finite Element Method:
  Its Basis and Fundamentals*, 7th ed., Butterworth-Heinemann, 2013. Used for
  quadratic tetrahedron shape functions in volume coordinates and Gauss rules.
  Chapter numbers are **UNVERIFIED**.
* T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic
  Finite Element Analysis*, Prentice-Hall 1987 / Dover 2000. Used for the patch test
  and isoparametric formulation. Section numbers are **UNVERIFIED**.
* B. M. Irons & A. Razzaque, "Experience with the patch test for convergence
  of finite elements", in *The Mathematical Foundations of the Finite Element
  Method* (ed. A. K. Aziz), Academic Press, 1972, 557–587. Pages are **UNVERIFIED**.
* G. Kirsch, "Die Theorie der Elastizität und die Bedürfnisse der
  Festigkeitslehre", *Zeitschrift des Vereines deutscher Ingenieure* 42 (1898)
  797–807. Verified (secondary record).

## 2. Element
Volume coordinates `L1..L4` (`L1 = 1 − ξ − η − ζ`, `L2 = ξ`, `L3 = η`, `L4 = ζ`).
Corner nodes 1–4, mid-edge nodes 5:(1,2) 6:(2,3) 7:(1,3) 8:(1,4) 9:(2,4) 10:(3,4).
Shape functions: corners `N_i = L_i(2L_i − 1)`, mid-edge `N_ij = 4 L_i L_j`.
Isoparametric `x = Σ N_a x_a`, `J = ∂x/∂ξ`, `B` is the 6×30 strain–displacement matrix
in Voigt order `(εxx, εyy, εzz, γyz, γxz, γxy)`.
Isotropic `D` from `E, ν`. Stiffness `K_e = ∫ Bᵀ D B dV`, computed with the 4-point rule
(weights 1/24, points `(a,b,b,b)` permutations, `a = 0.5854101966249685`,
`b = 0.1381966011250105`). This rule is exact for degree-2 polynomials, and therefore exact
for straight-edged tet10, where `B` is linear. Reject elements with `det J ≤ 0`.

Loads: consistent nodal forces of surface tractions on 6-node triangular
faces. Use a 6-point (degree-4) triangle rule, or 3-point for constant traction.
Dirichlet BCs are applied by row/column elimination.
Sparse assembly (`scipy.sparse.coo → csr`), solved with `spsolve`.

Stress recovery: evaluate `σ = D B u` of each element at its own nodes
(natural coordinates of the node), then average over elements sharing the node
(unweighted). Also provide stresses at Gauss points.

## 3. Meshers
* `box_mesh(Lx, Ly, Lz, nx, ny, nz)`: hexes split into 6 tets (Kuhn
  split, conforming), then quadratic mid-edge nodes deduplicated by sorted
  edge key.
* `plate_hole_quarter_mesh(a, W, t, n_r, n_θ, n_z, grading)`: quarter plate
  `[0,W]²` minus a hole of radius `a`. Structured O-grid built from two blocks
  (θ ∈ [0,45°] mapping to edge x=W and θ ∈ [45°,90°] mapping to edge y=W). Radial
  geometric grading `grading` concentrates elements at the hole. Hexes are split as above.
* `perturb_interior(mesh, amp, seed)`: random displacement of interior corner
  nodes. Mid-edge nodes are re-placed at edge midpoints (straight edges kept).

## 4. Kirsch solution (exact, used as boundary traction)
Uniaxial `σ∞` along `x`, hole radius `a`, polar `(r, θ)` from the hole center:

    σ_rr = σ∞/2 (1 − a²/r²) + σ∞/2 (1 − 4a²/r² + 3a⁴/r⁴) cos 2θ
    σ_θθ = σ∞/2 (1 + a²/r²) − σ∞/2 (1 + 3a⁴/r⁴) cos 2θ
    σ_rθ = −σ∞/2 (1 + 2a²/r² − 3a⁴/r⁴) sin 2θ

At `r = a`, `θ = 90°`: `σ_θθ = σ_xx = 3σ∞`. This is the stress concentration factor 3.
The solution is plane stress / plane strain with traction BCs and is independent of
`ν`. With `ν = 0` a thin 3D plate with `u_z = 0` on `z = 0` and free faces
has `σ_zz ≡ 0` and the same in-plane field **exactly**. So applying the exact Kirsch
traction on the outer edges `x = W` and `y = W` of a finite quarter plate makes
the exact 3D solution equal to Kirsch. The test therefore checks pure discretization error.
Symmetry BCs: `u_x = 0` on `x = 0`, `u_y = 0` on `y = 0`.

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A8.1 | Linear elasticity, small strain | ↓ near yield |
| A8.2 | Straight-edged elements (curved hole approximated by chords of quadratic edges placed on the circle at corners and midpoints) | ? — converges with refinement |
| A8.3 | Nodal averaging of stress | ↓ peak stress at coarse meshes |

## 6. Implementation contract (`cvtsim/fe/tet10.py`, `cvtsim/fe/mesh.py`)
* `Mesh(nodes (n,3), elems (m,10))`; `box_mesh`, `plate_hole_quarter_mesh`, `perturb_interior`.
  Mid-edge nodes on the hole are projected onto the circle. A8.2 then applies only to curved-edge
  Jacobians, so the test meshes are kept nearly straight.
* `stiffness(mesh, E, nu) -> csr`, `traction_load(mesh, faces, traction_fn) -> F`,
  `solve(K, F, fixed_dofs, fixed_vals) -> u`, `nodal_stress(mesh, u, E, nu) -> (n,6)`,
  `boundary_faces(mesh, predicate)`.

## 7. NUMERIC TEST TARGETS (`tests/test_fe.py`)
| Test | Setup | Target | Tol |
|---|---|---|---|
| T8.1 patch test | unit cube, `box_mesh` 3×3×3 perturbed (amp 0.15 of spacing), E=1, ν=0.3. Prescribe on all boundary nodes `u = A x + c` with a random A | interior `u` exact; all element Gauss stresses equal `D·sym(A)` | abs 1e-10 (u), rel 1e-8 (σ) |
| T8.2 patch, quadratic field | same mesh, `u` quadratic in x with matching body force: skipped. Instead: rigid-body modes give zero strain energy (6 zero eigenvalues of free `K` on a 1-element mesh) | 6 eigenvalues < 1e-10·max | – |
| T8.3 cantilever | L=10, h=w=1, E=1000, ν=0, root clamped, tip shear P=1 as uniform traction on end face. Meshes nx = 10, 20, 40 (ny=nz=2,2,4) | tip deflection → `PL³/3EI + PL/(κGA) = 4.024` (κ=5/6) | finest rel 2e-2; errors decrease monotonically |
| T8.4 Kirsch SCF | a=1, W=5, t=0.2, ν=0, σ∞=1, exact Kirsch traction on outer edges. 3 meshes (n_r,n_θ) = (8,8), (16,16), (32,32) with n_z=1 | `σ_xx(0, a) = 3.0` | finest within 3 %; error decreases monotonically with refinement |
| T8.5 Kirsch field | finest mesh, σ_θθ along the hole vs exact at 5 angles | rel 5 % | – |
| T8.6 element volume | single tet10 with straight edges | `Σ w det J = V` exact | 1e-14 |
