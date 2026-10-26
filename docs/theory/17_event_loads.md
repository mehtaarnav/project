# 17 — Event loads, durability and replays (`cvtsim.event_loads`)

Turns an event on the fast model (16) into loads and durability, then replays the event's worst
moments in the full model. It also compares set-ups and final drives on an event.

## 1. Loads over the whole event (fast model)
The fast model gives:
* the belt pull `F_b`;
* the mean span tension `T̄'` (without the centrifugal part, which does not load the hub);
* the pitch radii.

The force on the primary hub, with the spans at angle `γ = asin((r_s − r_p)/C)` to the centre
line, is
    F_x = 2 T̄' cos γ,    F_y = F_b sin γ
The secondary carries the equal and opposite force.

The section torques are:
* primary PTO: `F_b r_p`;
* secondary: `F_b r_s`, split half to each sheave half (15 §3).

These feed the shaft model (07) and launch durability (12 §7) unchanged. The results are fatigue
damage and Goodman factors, bearing equivalent loads, L10 life in hours at the event's mean
speed, and the PTO side load.

**Limit (A17.1):** these are mean-torque loads. Firing pulses and PTO twist are absent, so the
torque range, and therefore the crank PTO fatigue, is underestimated. Bending per revolution and
bearing life are dominated by the belt pull, which the fast model carries.

`--hours` scales per-event damage to a duty, for example 4 h of endurance.

## 2. Replays of the worst moments (full model)
1. **Pick windows.** Take the `k` times where the fast model's PTO side load is highest, at least
   3 s apart.
2. **Replay each window in the full model** (03–11), driven by the fast model's throttle and
   grade history over the window:
   * **Peak after 3 s:** the full model is warm-started from the fast state 0.6 s before the peak
     (§2.1). It runs 1 s, and the first 0.2 s is discarded as start-up transient.
   * **Peak in the first 3 s:** the window is the launch, replayed from rest.
3. **Report for each window:** peak PTO torque (with pulses), PTO side load, belt tension, bearing
   load, and the Goodman and yield factors of both shafts.

### 2.1 Warm start
The belt is laid on the open path at the free pitch radii of `x_p` and `x_s`, moving at
`r_p ω_e`. Shaft speeds, crank angle, sheave positions, vehicle speed and belt temperature are
set from the state.

Check (full launch, CVTech 0651-1100, restarted at 4.0 s from its own state; 0.2–0.8 s window
against the uninterrupted run):

| Quantity | Warm start | Uninterrupted |
|---|---|---|
| Peak PTO torque | 216 N·m | 212 N·m |
| Peak hub force | 1407 N | 1372 N |
| Peak belt tension | 1418 N | 1407 N |
| Mean engine speed | 333.0 rad/s | 334.0 rad/s |

Instantaneous values differ because the belt-vibration phase is not restored. Peaks and
statistics agree within 3 %.

## 3. Sweep
`sweep` runs every tune × final drive on one event. It reports:
* time, mean rpm, belt slip and corner overspeed;
* PTO side load, PTO-side bearing equivalent load and jackshaft Goodman factor (mean-torque, §1).

## 4. Targets
| ID | Case | Target |
|---|---|---|
| T17.1 | Hub force at `r_p = r_s` | `F_x = 2 T̄'`, `F_y = 0` |
| T17.2 | Hub force, `r_s − r_p = C/2` | `γ = 30°`: `F_x = √3 T̄'`, `F_y = F_b/2` |
| T17.3 | Window picking | the `k` largest values, at least `sep` apart |
| T17.4 | Warm start against the uninterrupted run | window peaks within 5 % (§2.1; slow test) |
| T17.5 | Sweep with two final drives | each row uses its own final drive; mean rpm rises with the reduction |

## 5. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A17.1 | Event fatigue and bearing life from mean torque | ↓ crank PTO torsion range (replays give the peaks) |
| A17.2 | Replay windows chosen by the fast model's PTO side load | ? a peak that only the full model shows (belt vibration) can be missed |
| A17.3 | Replays use traction 1 | ↑ on rough sections the full model sees more grip than the course |
