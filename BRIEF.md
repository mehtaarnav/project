Build `cvtsim`, a Python CVT dynamics + FE simulator, in this repo (mehtaarnav/project).

CONTEXT: I'm a UW-Madison Baja SAE student designing a rubber-belt mechanical CVT: a CVTech-style primary with flyweights on ramps, and a secondary with a helix cam plus torsion and linear springs, on a Kohler CH440 engine. A full rubber-belt transient FE model in Ansys was too slow, struggled to converge at the shift, and its belt is uncalibrated. I want our own purpose-built simulator that does NOT use Ansys. The goal is dynamic loads on shafts and bearings during shifting, then stresses. Budget: $100 of credits for the whole build. It runs overnight, unattended.

YOU ARE THE ORCHESTRATOR AND MATH AUTHORITY. Do the math yourself. Delegate implementation to subagents with model "haiku" to save budget. Track usage and keep about 15% in reserve; if you're running low, stop cleanly and write the report.

STACK: Python 3.11+, NumPy, SciPy, Numba (optional), pytest, matplotlib. GitHub Actions CI runs pytest on every push. Layout: cvtsim/{core,belt,contact,cvt,drivetrain,fe,io}, tests/, docs/theory/, examples/, REPORT.md. Make small commits, only with tests passing.

PHASE 0 (you write it, before any code): docs/theory/, one file per module. Each file has the equations with derivations, sources cited by paper or book and section, a symbol and unit table, frame conventions, assumptions with the direction each one biases results, and NUMERIC TEST TARGETS computed independently. Verify every citation; if you can't verify one, label it UNVERIFIED.
1. Rigid-body dynamics + joints (revolute, prismatic, screw/helix): constrained Newton–Euler with constraint stabilization. Tests: pendulum period (elliptic-integral reference) and energy conservation; torque-free top precession; constraint drift bound.
2. Time integration: generalized-alpha (Chung & Hulbert 1993) and/or symplectic explicit. Tests: convergence order, energy drift, stability limit.
3. Segmented belt: axial, bending and transverse-compression stiffness plus damping. Tests: capstan / Euler–Eytelwein tension ratio at slip; string natural frequency of a tensioned span; static stretch against EA.
4. V-belt contact with analytic cone sheaves (no mesh contact search): penalty normal force plus regularized Coulomb friction. Tests: F_a = T cos(beta)/(2 mu R) at incipient slip; convergence as penalty stiffness rises; thrust ratio and slip trends against Gerbert, Traction Belt Mechanics (cite precisely).
5. CVT primary (flyweights on ramps vs spring) and secondary (helix cam torque feedback plus springs). Tests: engagement rpm from the static force balance; secondary thrust from the cam relation.
6. Drivetrain: engine torque curve plus inertia, vehicle inertia through the final drive, road load. Test: steady-state speed balance.
7. Shafts: Timoshenko beams on bearings; stress with Peterson stress-concentration factors. Tests: hand statics reactions, beam-formula deflection.
8. Linear static FE with 10-node tetrahedral elements. Tests: cantilever tip deflection, Kirsch plate-with-hole stress concentration about 3 with mesh convergence, patch test.

Also docs/theory/literature.md, covering: Srivastava & Haque's review of belt/chain CVT dynamics (Mechanism and Machine Theory 2009); Gerbert; Chung & Hulbert; Shabana (ANCF, as a future belt option); Tasora/Anitescu and Project Chrono, and MBDyn. Verify each one.

PHASE 1: Haiku subagents in parallel, one per module group: core+integrator / belt+contact / cvt+drivetrain / shafts+fe. Each gets ONLY its spec files and must make every spec test pass. You review each result against the spec — signs, units, frames, energy — before merging. Reject anything that weakens a test.

PHASE 2: examples/launch.py runs a full-throttle launch from idle through upshift. It outputs the shift curve (engine rpm vs vehicle speed), ratio, sheave positions, belt tension, slip, bearing reactions, shaft torque and bending, and shaft stress over time, saved as CSV and PNG plots. Include time-step and penalty-stiffness convergence checks and an energy-balance check.

PHASE 3: REPORT.md with:
- Pass/fail for every test, with numbers.
- The launch results, with plots.
- An assumptions table with the direction of bias for each.
- What is uncalibrated: belt stiffness, friction, damping.
- What is not built yet: flexible bodies via Craig–Bampton, FE on real CAD parts.
- Budget used.

PARAMETERS (cvtsim/io/dr26_default.yaml). Keep each value's source note:
- Engine: Kohler CH440, team-restricted chart. Torque in N·m at rpm: 2400: 25.08, 2600: 24.54, 2800: 23.59, 3000: 22.51, 3200: 20.88, 3400: 19.66, 3600: 18.30. Below 2400 rpm the data is unverified, so flag it.
- Ratio range 3.0 (low) to 0.8 (overdrive). Do NOT use CVTech's 0.43.
- Primary: spring 13208.0 N/m, preload 683.04 N, 3 flyweights at 0.275 kg each plus cover 0.0387 kg each, flyweight radius 44.45 mm, ramp 45 deg (straight-ramp assumption), stroke 25.48 mm.
- Secondary: linear spring 3353.66 N/m, torsion spring 3.930 N·m/rad with 30 deg preload, equivalent preload 435.63 N, cam radius 38.1 mm, cam angle 45 deg, stroke 24.48 mm.
- Belt: CVTech BD52-2172-S proxy. Rubber shear modulus 30.15 MPa hot / 46.71 MPa cold (Messick 2018 thesis, small strain). Density 1150 kg/m³ (assumed). EA about 243.5 kN (unverified). Friction mu 0.3 (assumed). Sheave contact angle 45 deg as modeled; the CAD shows about 43.3 deg on some faces, so flag the discrepancy.
- Vehicle: mass 250 kg, final drive 8.1, tire radius 0.2921 m, rolling resistance coefficient 0.015.
- Damping 10 N·s/m on equivalent springs (assumed; run a sensitivity check).
- Physics flag: the flyweight force grows with sheave closure, so net primary axial stiffness turns negative above about 2,260 rpm. Model it through the real ramp kinematics, not as an equivalent spring.

RULES:
- Never fake or loosen a test.
- Never invent a number. Derive it, cite it, or flag it UNVERIFIED.
- No Ansys anywhere.

Done means: CI is green, launch.py runs end to end, and REPORT.md is complete and honest.
