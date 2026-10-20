# 07 — Shafts: Timoshenko beams on bearings, stress (`cvtsim.fe.shaft`)

## 1. Sources
* J. S. Przemieniecki, *Theory of Matrix Structural Analysis*, McGraw-Hill,
  1968 (Dover reprint 1985). Source of the shear-deformable beam stiffness
  matrix with `Φ = 12EI/(κGAL²)`. Section number **UNVERIFIED**.
  The element is exact for nodal loads on uniform segments. Tests check this
  against the closed forms in §4.
* G. R. Cowper, "The shear coefficient in Timoshenko's beam theory",
  *J. Appl. Mech.* 33(2) (1966) 335–340, doi:10.1115/1.3625046. Verified
  record. Solid circle: `κ = 6(1+ν)/(7+6ν)`. Table number **UNVERIFIED**.
* W. D. Pilkey & D. F. Pilkey, *Peterson's Stress Concentration Factors*,
  3rd ed., Wiley, 2008, ISBN 978-0-470-04824-5. Verified record. The
  curve-fit coefficients for the shoulder fillet (§3) are reproduced
  **from memory of the chart fits, UNVERIFIED**. They must be checked
  against the book (bending: Chart 3.10/3.11 area; torsion: Chart 3.12 area;
  chart numbers also UNVERIFIED) before design use. The finite-width hole fit
  is cross-checked against Kirsch's exact value 3 (08_fe_tet10.md).

## 2. Beam model
Each shaft lies along its own local `z`. Nodes sit at stations `z_k`. Each node
has 6 DOFs: `u_x, u_y` (lateral), `ψ_x, ψ_y` (cross-section rotations,
defined so that `ψ_x ≈ du_x/dz` and `ψ_y ≈ du_y/dz`), `u_z` (axial), and `φ_z` (twist).
Each lateral plane is an independent planar Timoshenko problem
(`EI`, `κGA` per element). Element stiffness for DOFs `(u1, ψ1, u2, ψ2)`:

    K = EI/((1+Φ)L³) · [[ 12,    6L,       −12,   6L      ],
                        [ 6L, (4+Φ)L²,     −6L, (2−Φ)L²   ],
                        [−12,   −6L,        12,  −6L      ],
                        [ 6L, (2−Φ)L²,     −6L, (4+Φ)L²   ]],   Φ = 12EI/(κGAL²)

Axial bar `EA/L [[1,−1],[−1,1]]`. Torsion `GJ/L [[1,−1],[−1,1]]`, where `J = πd⁴/32`
(solid) or `π(d_o⁴ − d_i⁴)/32`.

**Bearings:** linear springs at nodes. Radial `k_r` acts in both planes, axial `k_a`
on the thrust bearing only, and torsion is grounded where specified (e.g.
output sprocket held by the load in a static solve). Moment stiffness is zero
(simple supports, A7.2). `k_r → ∞` (default 1e12 N/m) reproduces
rigid supports. **Bearing reaction** = force the bearing applies on the shaft,
`R = −k u`. The **bearing load** reported is the force the shaft applies on the bearing, `−R`.

**Internal forces** at a station come from element end forces `K_e u_e`: shear `V`,
bending moments `M_x, M_y`, torque `T`, axial `N`. Resultant moment
`M = sqrt(M_x² + M_y²)`.

## 3. Nominal and peak stress
Solid section of diameter `d`:

    σ_b = 32 M/(π d³),   τ_t = 16 T/(π d³),   σ_a = 4 N/(π d²)
    σ_vm,peak = sqrt((K_tb σ_b + K_ta σ_a)² + 3 (K_ts τ_t)²)

This conservatively combines peak factors at a common point (A7.3).

**Shoulder fillet** (`D/d`, fillet radius `r`, `t = (D−d)/2`), valid for
`0.1 ≤ t/r ≤ 2` (bending) and `0.25 ≤ t/r ≤ 4` (torsion). The fit is
`K_t = C1 + C2 (2t/D) + C3 (2t/D)² + C4 (2t/D)³`. **UNVERIFIED coefficients:**

| Load | C1 | C2 | C3 | C4 |
|---|---|---|---|---|
| bending | 0.947 + 1.206√(t/r) − 0.131 t/r | 0.022 − 3.405√(t/r) + 0.915 t/r | 0.869 + 1.777√(t/r) − 0.555 t/r | −0.810 + 0.422√(t/r) − 0.260 t/r |
| torsion | 0.905 + 0.783√(t/r) − 0.075 t/r | −0.437 − 1.969√(t/r) + 0.553 t/r | 1.557 + 1.073√(t/r) − 0.578 t/r | −1.061 + 0.171√(t/r) + 0.086 t/r |

Outside the validity range: raise `ValueError`.
**Finite-width plate with central hole, tension (net-section Kt), UNVERIFIED fit**:
`K_tn = 2 + 0.284(1 − d/H) − 0.600(1 − d/H)² + 1.32(1 − d/H)³`. As `d/H → 0` this
gives 3.004, close to Kirsch's 3 (cross-check).

## 4. Closed-form test targets (steel E = 200 GPa, ν = 0.3, d = 25 mm)
* `κ = 6(1+ν)/(7+6ν) = 0.886364`.
* Simply supported span L = 0.3 m, P = 1000 N at midspan:
  reactions 500 N each. `δ_EB = PL³/(48EI) = 1.46677e-4 m`.
  `δ_T = δ_EB + PL/(4κGA) = 1.48918e-4 m`.
* Overhung: bearings at z = 0 and z = a = 0.2 m, load P = +1000 N (+y) at
  z = a + c = 0.25 m. Bearing reactions on the shaft: `R1 = +P c/a = +250 N`,
  `R2 = −P (a+c)/a = −1250 N`. Tip deflection: `δ_EB = P c²(a+c)/(3EI) = 5.43249e-5 m`.
  Timoshenko: `δ = δ_EB + P c (a+c)/(a κGA) = 5.61923e-5 m`.
  (Shear adds `P c/(κGA)` in the overhang and `(Pc/a)·c/(κGA)` through the
  span-section rotation.)
* Torsion: T = 50 N·m over 0.3 m: twist `TL/(GJ) = 5.08481e-3 rad`.
* Stress: M = 100 N·m, T = 50 N·m, Kt = 1: `σ_b = 65.1899 MPa`,
  `τ = 16.2975 MPa`, `σ_vm = 71.0390 MPa`.

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A7.1 | Quasi-static shaft response to dynamic loads (shaft modes ≫ load frequencies) | ↓ if any load content approaches the first bending mode (checked: modes reported) |
| A7.2 | Bearings as simple supports (no moment stiffness) | ↑ span moments for stiff bearing pairs, ↓ bearing-edge moments |
| A7.3 | Peak factors combined at one point | ↑ (conservative) |
| A7.4 | Shaft geometry, bearing positions, Kt geometry: PLACEHOLDERS (not supplied) | ? — results illustrative until CAD values are entered |
| A7.5 | Sheave hub load applied at a single node (no hub stiffening) | ↑ local bending at hub |

## 6. Implementation contract (`cvtsim/fe/shaft.py`, `cvtsim/fe/peterson.py`)
* `Shaft(stations_z, diameters, E, nu, kappa=None)`: diameters per element (stepped). Methods:
  `add_bearing(z, k_r=1e12, k_a=0, k_t=0)`, `solve(loads) -> ShaftResult`.
  `loads` is a list of `(z, Fx, Fy, Fz, Mx, My, Tz)` at nodes.
  `ShaftResult` has `.u`, `.reactions` (per bearing: `Fx, Fy, Fz, Tz` on the shaft),
  `.internal(z) -> (Vx, Vy, Mx, My, T, N)`, and `.stress(z, Ktb=1, Kts=1, Kta=1) -> (sig_b, tau, sig_a, sig_vm)`.
  For load cases that cannot be equilibrated (e.g. torque with no torsional ground), the solver adds
  `k_t` at the specified torque-reaction station (the output sprocket).
  Stations must include every bearing and load `z` (the constructor merges them).
* `influence_matrix(...)`: linear response (reactions, internal forces at chosen stations) to unit
  loads. Used for fast superposition in Phase 2.
* `peterson.py`: `kt_shoulder_bending(D, d, r)`, `kt_shoulder_torsion(D, d, r)`, `kt_hole_plate_tension(d, H)`.
* `first_bending_frequency(shaft, mass_per_length, lumped=[(z,m),…])`: lumped-mass eigenvalue (consistent with A7.1).

## 7. NUMERIC TEST TARGETS (`tests/test_shaft.py`)
| Test | Target | Tol |
|---|---|---|
| T7.1 SS beam reactions | 500, 500 N | rel 1e-9 |
| T7.2 SS midspan deflection (Timoshenko) | 1.48918e-4 m | rel 1e-5 |
| T7.3 same with κGA → ∞ (κ=1e9) | 1.46677e-4 m | rel 1e-5 |
| T7.4 overhung reactions | R1 = +250, R2 = −1250 N | rel 1e-9 |
| T7.5 overhung tip deflection | 5.61923e-5 m (EB limit 5.43249e-5) | rel 1e-5 |
| T7.6 twist | 5.08481e-3 rad | rel 1e-6 |
| T7.7 stresses (M=100, T=50) | σ_b 65.1899 MPa, τ 16.2975 MPa, σ_vm 71.0390 MPa | rel 1e-5 |
| T7.8 `kt_shoulder_bending(1.2, 1.0, 0.1)` | 1.66575 (fit evaluation; UNVERIFIED vs chart) | 1e-4 |
| T7.9 `kt_hole_plate_tension(1e-6, 1)` | 3.004 | 1e-3 |
| T7.10 `kt_shoulder_bending(1.2,1.0,r)` and `kt_shoulder_torsion` strictly decreasing in r over r ∈ [0.05, 0.5]·d where valid; all values > 1 | qualitative | – |
| T7.11 equilibrium: Σ reactions + Σ loads = 0 (forces and moments) for a random load set | 1e-9 relative | – |
