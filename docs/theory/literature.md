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
| J. Baumgarte, "Stabilization of constraints and integrals of motion in dynamical systems", *CMAME* 1 (1972) 1–16, doi:10.1016/0045-7825(72)90018-7 | Verified | Constraint stabilization (01). |
| E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer, 2006 | Verified | Störmer–Verlet properties (02). |
| G. R. Cowper, "The shear coefficient in Timoshenko's beam theory", *J. Appl. Mech.* 33(2) (1966) 335–340 | Verified | Shear coefficient (07). |
| W. D. Pilkey, D. F. Pilkey, *Peterson's Stress Concentration Factors*, 3rd ed., Wiley, 2008 | Verified (record). Chart fits **UNVERIFIED** | Kt (07). |
| G. Kirsch (1898), *Z. VDI* 42, 797–807 | Verified (secondary record) | Plate-with-hole exact solution (08). |
| Messick (2018) thesis, rubber G hot/cold | **UNVERIFIED** (not found) | Belt modulus (03). |
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
