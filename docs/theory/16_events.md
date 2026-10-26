# 16 — Events: courses and driver (`cvtsim.events`)

Runs a car around a whole event on the fast CVT model (15). The three events are acceleration,
hill climb and endurance. Courses live in `cvtsim/io/tracks/`.

## 1. Course
A course is a list of segments `[length, grade at start, grade at end, corner radius, traction]`.
* Grade is in degrees and varies linearly along a segment.
* A radius of 0 means a straight.
* Traction scales tire grip and load, as in 14.

The course is sampled every 1 m and repeated for multi-lap runs. Grade and traction enter the
model as in 14, but by distance travelled, not by time.

## 2. Driver
* **Corner speed:** `v_c = sqrt(μ_t λ g R)`.
* **Target speed:** a backward pass, so the driver lifts off early enough to reach each corner at
  its limit:
      v_t(s_i) = min(v_c(s_i), sqrt(v_t(s_i+1)² + 2 a_lift (s_i+1 − s_i))),    a_lift = max(a_coast + g sin φ, 0)
  There are no brakes. The only deceleration is lift-off (`driver.a_coast`, ASSUMED 0.5 m/s²).
  The fast model gives 0.6–0.7 m/s² on level ground. On a steep descent `a_lift` is 0, and the
  driver cannot slow down.
* **Throttle:** `θ_cmd = clip(K_d (v_t − v), 0, 1)`, followed by a first-order lag
  `driver.throttle_lag`.
* **End of run:** the car reaches the end of the course, or the time limit (`t_limit` per lap)
  runs out.

## 3. Results
For each event:
* time, or lap time;
* average and peak speed;
* engine rpm (mean, 10th and 90th percentile while moving);
* fraction of time the belt slips;
* belt temperature;
* engine work;
* ratio range used;
* wheel-spin fraction;
* `corner_overspeed_kmh`, how much faster than the corner limit the car arrived. With no brakes,
  any overspeed means the lap is faster than a real car could drive it.

## 4. Supplied courses
| Course | Content | Status |
|---|---|---|
| `acceleration` | 150 ft (45.72 m), flat | event definition; check the current rules (some events use 100 ft) |
| `hill_climb` | 15 m run-up, then 20° for 50 m | invented; with the default grip the car is traction-limited above about 22° |
| `endurance` | 1.49 km lap: corners of 8–30 m radius, climbs to 15°, descents to 12°, rough sections | invented, representative only |

## 5. Targets
| ID | Case | Target |
|---|---|---|
| T16.1 | Course of two segments, 10 m flat then 10 m from 0 to 10° | 21 samples; grade at s = 15 m is 5° |
| T16.2 | Target speed on the flat before a corner of radius R | `v_c = sqrt(μ g R)`; at distance d before it, `sqrt(v_c² + 2 a d)` |
| T16.3 | Acceleration event against the time scenario | time to 45.72 m equals the time the full-throttle launch scenario reaches 45.72 m, within one step |
| T16.4 | Energy over an endurance lap | closure within 1e-3 of engine work |

## 6. Assumptions
| # | Assumption | Bias |
|---|---|---|
| A16.1 | No brakes; deceleration only by lift-off | ↓ lap time too short where corner overspeed > 0 |
| A16.2 | 1-D car: corners only set the target speed; no lateral load transfer or cornering drag | ↓ lap time |
| A16.3 | Driver tracks the target with a proportional throttle | ? |
| A16.4 | Invented hill and endurance courses | results compare set-ups; they do not predict event times |
