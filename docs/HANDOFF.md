# Handoff — state of cvtsim on 2026-10-20

Read this first in a new session. Branch `claude/peaceful-gates-aeobqu` holds all work; `main` is
older and has not been updated since the first merge.

## Rules the owner set
* Commits authored `mehtaarnav <as4622913@gmail.com>`. No co-author lines or AI attribution
  anywhere: commits, code, docs. History was rewritten on 2026-10-20 to remove two such lines.
* Never invent a number. Every parameter has a `source:` note: a citation, `measured`, `DERIVED`,
  `ASSUMED`, `UNVERIFIED` or `PLACEHOLDER`.
* Never loosen a test to pass. Correct a target only with a written derivation.
* Small commits, only with passing tests.
* Plain, professional writing. No jargon tags in prose meant for people.
* Do not commit the SAE papers or the team workbook and MATLAB model (`reference/` and
  `cvtsim/io/CVT_Inputs.xlsx` are git-ignored).
* Licence: MIT.

## What exists
* **Model:**
  * segmented belt with cone-sheave contact;
  * flyweight primary and helix secondary, with friction;
  * engine firing pulses and PTO torsion;
  * tire and vehicle.
* **Scenarios:** launch, hill climb, backshift, landing (`docs/theory/14_scenarios.md`).
* **Durability:** fatigue by revolutions, bearing equivalent load and required rating, PTO side
  load (`docs/theory/12_fatigue_bearings.md` §7). Audited by hand with
  `examples/audit_durability.py`.
* **Uncertainty:**
  * `python -m cvtsim.uq screen` ranks the inputs;
  * `python -m cvtsim.uq bands --top ...` gives ranges over the assumed input ranges.
* **Tools:**
  * `python -m cvtsim.app`: parts picker in the browser;
  * `python -m cvtsim.tune`: tracked runs;
  * `python -m cvtsim.inputs --open`: inputs that still need data;
  * `docs/YOUR_CAR.md`: how another team models its car.
* **Viewer:** `viewer/index.html`. Its data format, including per-frame mechanism state and event
  times for animation, is in `docs/viewer_contract.md`.

## Latest results (8 s launch, run `runs/20261020-211533-tune-default`)
| Result | Value |
|---|---|
| Engagement | 1914 rpm |
| Shift complete | 5.38 s |
| Speed at 8 s | 56.9 km/h (measured band 53–63) |
| Peak belt tension | 1623 N |
| PTO side load | 1408 N (85 N·m at the mounting face) |
| Goodman factor, crank PTO / jackshaft | 1.44 / 1.07 (simultaneous stresses; recomputed after the review fixes, the run files still show the old 1.86 / 0.81) |
| Yield factor, crank PTO / jackshaft | 3.8 / 3.2 |
| PTO main bearing | equivalent load 1089 N from the belt only; engine gas and inertia loads are not included |

The jackshaft and crank numbers rest on placeholder shaft geometry and material.

Sensitivity ranking for the PTO loads, largest first:
1. flywheel inertia;
2. clutch inertia;
3. belt bending stiffness;
4. engagement gap;
5. belt friction.

## Fast model (Phase A of the lap simulator)
`cvtsim.fast` (spec 15) runs an 8 s launch in about 10 ms.
* Calibrated on the full-model launch and hill climb.
* Validated on full-model runs it was not fitted to (backshift scenario; CVTech 0651-1100 spring),
  passing all acceptance checks: rpm within 3.6 %, primary travel within 1.1 mm, final speed
  within 0.7 %, engagement within 0.6 %.
* One fit-run miss: hill-climb final speed is 3.3 % off against the 3 % limit.
* The calibration runs used two different friction smoothings: 0.05 m/s for the launch and hill
  climb, 0.0125 m/s for the backshift and spring runs.

Next phases:
* **B:** tracks and driver for acceleration (150 ft), hill climb and endurance.
* **C:** per-lap results and durability, and part sweeps.
* No brakes yet, by the owner's choice.

## Open work, in order
1. **Done:** the sensitivity bands (`runs/20261020-215336-uq-bands`, 30 samples, first 2 s) and a
   hill climb (`runs/20261020-213705-tune-default-hill_climb`) are committed.
   * The bands predate the review fixes and used a stale screen; rerun both before quoting them.
   * 5–95 % ranges: peak bearing load 2.8–5.6 kN, peak PTO side load 1.0–1.5 kN, jackshaft Goodman
     factor 0.50–0.99.
   * The upper tail comes from the engagement and lock-up transient.
   * Long runs only survive a container restart when started as session background tasks.
2. **Rewrite REPORT.md from the current baseline.** It still describes the model before firing
   pulses, friction and durability.
3. **Run each scenario once.** Add their peak loads to the report. The design load is probably the
   landing or the hill climb, not the launch.
4. **Merge to `main`.** Ask the owner first.
5. **Data still needed from the team or Rehlko:**
   * flywheel and clutch inertia;
   * crankshaft and jackshaft material;
   * jackshaft layout: steps, bearing seats, sprocket position;
   * bearing part numbers;
   * allowable PTO side load;
   * gearbox tooth counts (8.1 vs 8.0);
   * a logged launch for calibration.
6. **Handoff package corrections already applied:**
   * the 0.8 overdrive stop is labelled a brief requirement;
   * the hardware limit 0.4642 is a tune, `stock_unrestricted`;
   * the measured car (275 kg, 0.550 rear, 1.4224 m wheelbase) is a tune, `measured_car`;
   * the loader rejects unreachable ratios.

## Notes
* A remote branch `app` exists. Its work is merged into the working branch; it can be deleted.
* After the history rewrite, GitHub may keep listing Claude as a contributor for a while.
  Renaming `main` and renaming it back forces a refresh.
