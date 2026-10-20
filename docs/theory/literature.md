# Literature review and citation status

Verification means the bibliographic record (authors, title, venue, year,
volume/pages or publisher) was confirmed in an online catalog or publisher
record on 2026-10-20. We could not open full texts in this environment, so
**any section, equation, chart, or page claim not listed as verified is
UNVERIFIED**. Our tests never depend on a citation. Every target is derived
in the module spec and recomputed by `targets.py`.

| Reference | Status | Role in cvtsim |
|---|---|---|
| N. Srivastava, I. Haque, "A review on belt and chain continuously variable transmissions (CVT): Dynamics and control", *Mechanism and Machine Theory* 44(1) (2009) 19–41 | Verified (record) | Survey framing. Rubber-belt CVT models range from steady (Gerbert-type) belt mechanics to transient multibody/FE belt models, plus shift-dynamics models for control. cvtsim sits in the "transient multibody, discretized belt" class. |
| G. Gerbert, *Traction Belt Mechanics: Flat Belts, V-belts, V-ribbed Belts*, Chalmers Univ. of Technology, Göteborg, 1999 (WorldCat OCLC 247549000) | Verified (record). Content **UNVERIFIED** (text not accessed) | Steady V-belt theory: wedge friction `μ/sin β`, radial belt sliding in the groove, slip growing with load, thrust ratio. Used qualitatively only (04 §1). Its earlier "Force and slip behaviour in V-belt drives", *Acta Polytechnica Scandinavica* (1972), is also verified (record). |
| J. Chung, G. M. Hulbert, "A time integration algorithm for structural dynamics with improved numerical dissipation: the generalized-α method", *J. Appl. Mech.* 60(2) (1993) 371–375, doi:10.1115/1.2900803 | Verified | Implicit integrator (02). |
| A. A. Shabana, "An absolute nodal coordinate formulation for the large rotation and deformation analysis of flexible bodies", Tech. Rep. MBS96-1-UIC, Univ. Illinois at Chicago, 1996 | Verified (cited record) | **Future** belt option. ANCF gives a constant mass matrix and exact rigid-body motion, which suits a continuous belt with large rotation. Not implemented. |
| A. Tasora, M. Anitescu, "A matrix-free cone complementarity approach for solving large-scale, nonsmooth, rigid body dynamics", *Comput. Methods Appl. Mech. Eng.* 200 (2011) 439–453, doi:10.1016/j.cma.2010.06.030 | Verified | Alternative to our penalty contact: nonsmooth (complementarity) frictional contact. Allows large time steps with no penalty stiffness to tune. A future path if penalty stiffness limits the step size. |
| Project Chrono (projectchrono.org; Chrono::Engine, Tasora, Negrut, Serban et al.) | Verified (project site) | Open-source multibody + FEA (incl. ANCF beams, NSC/SMC contact). Candidate cross-validation code for the belt. Not used here. |
| MBDyn (Politecnico di Milano, Dept. of Aerospace Science and Technology; main developer P. Masarati; GPL since 2001) | Verified (project pages) | Open-source multibody with nonlinear beams. Alternative cross-check. Not used. |
| CVTech-IBC archived Baja belt selection table (via Baja SAE forum, search-engine extract; pages not reachable from this environment) | **UNVERIFIED (secondary)** | BD52-2167-S 835 mm @ 200 mm C/C, BD52-2172-S 860 mm @ 214 mm, BD52-2179-S 895 mm, BD52-2183-S 915 mm @ 242 mm. The 895 mm pitch matches the team workbook. |
| J. Martinez, D. Corona, S. Diaz, R. Boccardo, "Characterization of an Adjustable Centrifugal Continuously Variable Transmission for a Baja SAE Prototype", SAE Technical Paper 2010-01-0861 (2010), doi:10.4271/2010-01-0861 | **Verified (full text read)** | CVTech-IBC Baja CVT dynamometer characterization. Kit part numbers (flyweights 135/160/190 g; drive springs #1-#3; driven springs A/B), ratio range 0.241-1.357 (driven/drive speed), 235 mm centres, 992 mm belt, 0.97 s CVT stabilization time, measured 100 ft acceleration 4.10-4.54 s (Table 4), qualitative tuning trends (Table 3). Used for tune-library components and as the vehicle validation band (10_validation.md). |
| R. Spilde, Y. Liu, "System Design of Drivetrain for a Baja SAE Four-Wheel Drive Vehicle", SAE Technical Paper 2025-01-5069 (2025), doi:10.4271/2025-01-5069 | **Verified (full text read)** | Kohler Command Pro 14 Baja (CH440) 7.17 kW, 25.08 N·m at 2400 rpm (Fig. 2); measured top speed 36.7 mph = 59.06 km/h; other teams 33-39 mph; measured 100 ft acceleration 3.92 s (2024) and 5.28 s (2023); C_d ≈ 0.44 (their ref. 21), rotating tires up to 25 % of aero drag, aero 72.9 % of drive force at top speed (CFD); rolling resistance 0.01-0.02 dynamic, ≥ 0.02 static; CVT efficiency falls with speed. Used in 10_validation.md and vehicle defaults. |
| J. Baumgarte, "Stabilization of constraints and integrals of motion in dynamical systems", *CMAME* 1 (1972) 1–16, doi:10.1016/0045-7825(72)90018-7 | Verified | Constraint stabilization (01). |
| E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer, 2006 | Verified | Störmer–Verlet properties (02). |
| G. R. Cowper, "The shear coefficient in Timoshenko's beam theory", *J. Appl. Mech.* 33(2) (1966) 335–340 | Verified | Shear coefficient (07). |
| W. D. Pilkey, D. F. Pilkey, *Peterson's Stress Concentration Factors*, 3rd ed., Wiley, 2008 | Verified (record). Chart fits **UNVERIFIED** | Kt (07). |
| G. Kirsch (1898), *Z. VDI* 42, 797–807 | Verified (secondary record) | Plate-with-hole exact solution (08). |
| M. J. Messick, *An Experimentally-Validated V-Belt Model for Axial Force and Efficiency in a Continuously Variable Transmission*, M.S. thesis, Virginia Tech, 2018, hdl:10919/85055 | **Verified (full text read)** | CVTech BD52-2172-S, aramid-reinforced single-cog belt: G 46.71 MPa (25 °C) / 30.15 MPa (60 °C) by DMA, E'/E'' = 8.2 at 60 °C (§4.1); 27 Kevlar 29 yarns EA 220.5 kN + rubber 22.95 kN = 243.5 kN (§4.2-4.3); equivalent section 2.772e-4 m², EI 0.2602 N·m² with the neutral axis at the yarns (§4.4); measured wedge (radial) stiffness 28.98 MPa at 25 °C, 18.71 MPa at 60 °C on 12° blocks (§4.5); belt mass 375 g, 0.437 kg/m, cord line 2.85 mm below the outer surface, pitch length 857.2 mm (§4.6); CVTech sheave half-angle 12°; kinetic μ = 0.4 (Gerbert's 0.3-0.5, not measured; < 1 % efficiency effect over that range, §7.4); dyno-measured axial forces and efficiency (ch. 9). |
| Haug (1989); Przemieniecki (1968); Zienkiewicz–Taylor–Zhu (2013); Hughes (1987); Irons & Razzaque (1972); Johnson (1985); Wriggers (2006); Gillespie (1992); Landau & Lifshitz | Standard texts, records **not re-checked** → treat as UNVERIFIED | Background only. All formulas used are re-derived in the specs. |

## Where cvtsim sits
Srivastava & Haque's review splits CVT belt modeling into (i) steady
continuum belt mechanics (Gerbert-type), (ii) transient discretized-belt
multibody models, and (iii) reduced shift-dynamics models for control.
cvtsim is class (ii) with analytic cone sheaves and penalty/regularized
friction. That is the cheapest model that resolves *dynamic* sheave thrust, belt
span tensions and slip during a shift, which are the inputs to shaft and bearing
loads. Its main known weaknesses relative to the literature:
* No belt axial/twist DOFs. Shabana's ANCF or a Chrono-style 3D belt would add them.
* Penalty rather than complementarity contact (Tasora–Anitescu). That costs
  time-step and tuning, and causes some regularization creep.
* Uncalibrated rubber properties.
