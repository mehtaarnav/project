# Modelling your own car

cvtsim ships with the UW–Madison DR26 car. Any team with a rubber-belt CVT can model its own car
by copying one file and replacing the numbers.

## 1. Start from the DR26 file

```
cp cvtsim/io/dr26_default.yaml my_car.yaml
python -m cvtsim.inputs my_car.yaml --open     # every input that still needs data
```

Every entry has the form `{value: ..., unit: ..., source: "..."}`. Write where each number came
from in `source`. Start a note with `measured`, `DERIVED`, `ASSUMED` or `PLACEHOLDER`, so that
`cvtsim.inputs` can report what is still unknown.

Run your car:

```
python -m cvtsim.tune run default --car my_car.yaml                       # launch
python -m cvtsim.tune run default --car my_car.yaml --scenario hill_climb
```

## 2. What to measure, in order of effect

The uncertainty screen (`python -m cvtsim.uq screen`) ranks inputs by their effect on the results.
For the DR26 car, the order for shaft and bearing loads was the following:

| Input | Parameter | How to get it |
|---|---|---|
| Flywheel + crank inertia | `engine.J_f` | Bifilar or trifilar pendulum on the flywheel, or the engine maker |
| Clutch inertia (primary, without flyweights) | `engine.J_c` | CAD mass properties of the fixed half, spider and cover, or a pendulum test |
| Belt bending stiffness | `belt.EI` | Three-point bend test on a belt section |
| Engagement travel | `primary.engagement_gap` | Measure sheave travel from fully open to belt contact |
| Belt–sheave friction | `belt.mu` | Literature range 0.3–0.5 unless tested |

For launch time and top speed, also measure:
* vehicle mass with the driver: `vehicle.mass`;
* rear weight fraction: `vehicle.b_over_L`;
* CG height: `vehicle.h_cg`;
* tire radius: `vehicle.tire_radius`;
* gear reduction: `vehicle.final_drive` (count the teeth).

## 3. CVT geometry

| Parameter | Meaning |
|---|---|
| `geometry.belt_pitch_length`, `geometry.center_distance` | from the belt and clutch maker |
| `geometry.rp_mech_min`, `rp_mech_max`, `rs_mech_min`, `rs_mech_max` | smallest and largest belt pitch radius on each clutch, measured |
| `sheave.beta_deg` | sheave flank half-angle, measured |
| `ratio.high` | overdrive stop. The loader rejects a ratio the radii cannot reach and names the violated limit |

## 4. Parts you can swap

Springs, flyweights, ramps, helices and torsion springs live in
`cvtsim/io/tuning/components.yaml`. Add your parts there and they appear in the dropdowns of
`python -m cvtsim.app`. A tune file in `cvtsim/io/tuning/tunes/` selects one of each.

## 5. Shafts and bearings

Stress and bearing life need the real layout. In `shafts.*`, give:
* section stations and diameters;
* bearing positions;
* sheave and output positions;
* the critical shoulder (diameters and fillet radius).

In `durability.*`, give the shaft materials and bearing part numbers. Add material strengths and
bearing ratings to `cvtsim/io/materials.yaml` from the supplier's data sheet.

## 6. Engine

`engine.rpm` and `engine.torque` hold the dyno curve. Firing pulses need:
* bore and stroke;
* rod length;
* reciprocating mass.

Each is listed under `engine.*` with the value used for the Kohler CH440.
