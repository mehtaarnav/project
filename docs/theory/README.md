# cvtsim theory specification (Phase 0)

This directory is the math authority for `cvtsim`. Code must implement these
equations; tests must check the NUMERIC TEST TARGETS given here. Every target
was computed independently of `cvtsim` by `docs/theory/targets.py` (closed
forms only — run it to reproduce).

| File | Module | Package |
|---|---|---|
| [01_rigid_body.md](01_rigid_body.md) | rigid bodies + joints | `cvtsim.core` |
| [02_integrators.md](02_integrators.md) | time integration | `cvtsim.core.integrators` |
| [03_belt.md](03_belt.md) | segmented belt | `cvtsim.belt` |
| [04_contact.md](04_contact.md) | V-belt / cone-sheave contact | `cvtsim.contact` |
| [05_cvt.md](05_cvt.md) | primary + secondary clutches, geometry | `cvtsim.cvt` |
| [06_drivetrain.md](06_drivetrain.md) | engine, vehicle, road load | `cvtsim.drivetrain` |
| [07_shafts.md](07_shafts.md) | Timoshenko shafts, bearings, stress | `cvtsim.fe.shaft` |
| [08_fe_tet10.md](08_fe_tet10.md) | linear static 10-node tet FE | `cvtsim.fe.tet10` |
| [literature.md](literature.md) | literature review + verification status | — |

## Global conventions (all modules)

* **Units:** SI throughout (m, kg, s, N, N·m, rad, Pa). rpm only at I/O
  boundaries, converted with `omega = rpm*2*pi/60`.
* **World frame:** right-handed. Origin at the primary shaft axis in the belt
  plane. `+X` points from the primary axis to the secondary axis (the center
  distance is `C`). `+Y` completes the belt plane. `+Z` is parallel to both shaft axes.
  Positive rotation is about `+Z` (counter-clockwise seen from `+Z`). The
  engine turns the primary in the `+Z` sense.
* **Belt plane model:** belt nodes move in the X–Y plane only. Axial (Z)
  belt motion, from belt misalignment as one sheave half moves, is neglected.
* **Sheave axial coordinates:** `x_p` is primary *closing* travel of the
  movable half (0 = fully open, belt at minimum radius). `x_s` is secondary
  *opening* travel (0 = fully closed, belt at maximum radius). Both are bounded
  by stops at `[0, stroke]`.
* **Forces "on" a body** are the forces that body receives. Reactions follow
  Newton's third law explicitly in code; no sign is implied by context.
* **Bias notation:** each assumption lists the direction it biases results:
  ↑ = over-predicts, ↓ = under-predicts, ? = direction unknown.
* **Citation status:** a reference is *verified* when its bibliographic record
  was confirmed online during Phase 0 (2026-10-20). Section, equation, or chart
  numbers we could not open are marked **UNVERIFIED**.
* **Parameter status:** every numeric parameter in
  `cvtsim/io/dr26_default.yaml` carries a `source:` note. Values not supplied
  by the brief and not derivable are marked `PLACEHOLDER` and must be
  replaced from CAD/measurement. Results that depend on them are not
  quantitative predictions.

## Inconsistencies in the input data

1. **Sheave angle vs. stroke vs. ratio range.** With the half-angle
   `beta = 45°` (flank measured from the radial plane), axial travel `s` of one
   sheave half changes the belt radius by `s/(2 tan beta)`. That is 12.74 mm
   for the 25.48 mm primary stroke and 12.24 mm for the 24.48 mm secondary
   stroke. Requiring ratio 3.0 at one end and 0.8 at the other then
   *forces* `r_p` = 10.20 → 22.94 mm and `r_s` = 30.59 → 18.35 mm (05_cvt.md
   §2). A 10 mm primary pitch radius is not physically realizable with a
   CVTech-size belt. Either `beta` is not a 45° half-angle (a ~13° rubber-CVT
   half-angle would give ~55 mm of radial travel), or the strokes are not
   sheave strokes. The simulator follows the brief literally and flags this.
   Absolute belt tensions and shaft loads scale roughly with `1/r`, so treat
   them as uncertain by a factor of ~2–3 until this is resolved.
2. **Flyweight kinematic ratio.** A weight guided radially against one 45°
   ramp has `kappa = dr_w/dx_p = 1/tan(alpha) = 1`. That predicts engagement
   at 1220 rpm and negative net primary stiffness above 1131 rpm. The brief
   states "about 2,260 rpm". That is reproduced to 0.1 % (2262.6 rpm) by
   `kappa = 1/(2 tan alpha) = 0.5`, the kinematics of a roller wedged between
   two opposed 45° faces. We adopt `kappa = 0.5` (engagement 1725.7 rpm) and
   mark it UNVERIFIED against CAD (05_cvt.md §3).
