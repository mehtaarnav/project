# cvtsim theory specification

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
| [09_tuning_vehicle_belt.md](09_tuning_vehicle_belt.md) | profiled ramp/helix, tire vehicle, belt heating, tunes | `cvtsim.cvt.profiles`, `cvtsim.cvt.system`, `cvtsim.tune` |
| [11_engine_friction_uq.md](11_engine_friction_uq.md) | firing pulses, PTO torsion, clutch friction, uncertainty bands | `cvtsim.drivetrain.pulses`, `cvtsim.cvt.system` |
| [10_validation.md](10_validation.md) | measured launch and top-speed bands | — |
| [12_fatigue_bearings.md](12_fatigue_bearings.md) | rainflow, Goodman, Miner, ISO 281, launch durability | `cvtsim.fe.fatigue`, `cvtsim.fe.bearings`, `cvtsim.post.durability` |
| [13_belt_life_efficiency.md](13_belt_life_efficiency.md) | belt life and CVT efficiency | `cvtsim.post` |
| [14_scenarios.md](14_scenarios.md) | hill climb, backshift, landing | `cvtsim.cvt.system.schedule` |
| [15_fast_cvt.md](15_fast_cvt.md) | fast CVT model for whole events, calibration | `cvtsim.fast` |
| [16_events.md](16_events.md) | courses, driver, event results | `cvtsim.events` |
| [17_event_loads.md](17_event_loads.md) | event durability, full-model replays, sweeps | `cvtsim.event_loads` |
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
  movable half, measured from belt contact (negative over the engagement gap, down to
  `−engagement_gap` when fully open). `x_s` is secondary
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

1. **Sheave angle (resolved).** The brief's 45° half-angle forced pitch radii of 10–31 mm.
   The measured 12° (Messick 2018 Fig 4.8; stock primary 11.64°) reproduces the measured
   59.94 mm primary radius range from the 25.48 mm stroke exactly, and is now used. One residual
   mismatch remains: the secondary stroke implied by the measured radii is 23.22 mm vs 24.48 mm
   in the parameters (05_cvt.md §2).
2. **Flyweight kinematic ratio.** A weight guided radially against one 45°
   ramp has `kappa = dr_w/dx_p = 1/tan(alpha) = 1`. That predicts engagement
   at 1220 rpm and negative net primary stiffness above 1131 rpm. The brief
   states "about 2,260 rpm". That is reproduced to 0.1 % (2262.6 rpm) by
   `kappa = 1/(2 tan alpha) = 0.5`, the kinematics of a roller wedged between
   two opposed 45° faces. We adopt `kappa = 0.5` (engagement 1725.7 rpm) and
   mark it UNVERIFIED against CAD (05_cvt.md §3).
