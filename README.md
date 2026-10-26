# cvtsim

A dynamic simulator for rubber-belt CVTs in Baja SAE cars. It shows what the clutches, belt,
engine and car do during a launch, a hill climb, a backshift or a landing. It also gives the
loads on the engine PTO, the jackshaft and their bearings. No commercial solver is needed.

What it models:
* **Belt:** a flexible segmented belt with stretch, bending and damping, clamped by cone sheaves
  with friction. Slip, lock-up and belt heating come out of the model; they are not imposed.
* **Primary clutch:** flyweights on a profiled ramp against the spring, with roller and bushing
  friction.
* **Secondary clutch:** helix and torsion spring, with button friction.
* **Engine:** single-cylinder firing pulses on the dyno torque curve, and a torsionally flexible
  PTO between flywheel and clutch.
* **Car:** tire slip, weight transfer, rolling resistance, drag and grade.
* **Results:**
  * shaft bending and torque;
  * bearing loads;
  * PTO side load;
  * fatigue (rainflow, Goodman, Miner) and bearing life (ISO 281).

Every equation is written out with its sources in [docs/theory/](docs/theory/README.md), and is
tested against an independent hand calculation.

## Install

Python 3.11 or newer.

```
git clone <this repository> && cd <repository>
pip install -e .[fast,test]      # "fast" installs numba, which compiles the solver; install it
```

## Run

```
python -m cvtsim.app                                         # dropdowns for every part, in the browser
python -m cvtsim.tune run default                            # 8 s launch, about 30 min
python -m cvtsim.tune run default --scenario hill_climb      # or backshift, landing
python -m cvtsim.tune compare default cvtech_stock           # two set-ups, same event
python -m cvtsim.tune list                                   # tunes and scenarios available
```

Each run is saved under `runs/<date>-<name>/` and contains:
* the full parameter set and the git commit;
* metrics;
* the time series;
* a dashboard image;
* `viewer.json`.

Open the mechanism viewer to play a run back: start `python -m http.server`, then go to
`viewer/index.html?run=../runs/<run>/viewer.json`.

## Fast model

`cvtsim.fast` replaces the segmented belt with belt-length closure, a clamp-force balance and a
stick–slip belt clutch. It uses the same parts, engine, tire and scenarios as the full model. An
8 s launch takes about 10 ms instead of about 30 min, which makes whole events and part sweeps
practical.

Its few belt constants are fitted to full-model runs (`python -m cvtsim.fast.calibrate`). On
full-model runs it was not fitted to (a throttle lift, and a stiffer spring), it matches:
* engine speed within 3.6 %;
* primary travel within 1.1 mm;
* final speed within 0.7 %.

Use the full model for peak shaft and bearing loads; the fast model has no firing pulses and no
PTO twist. See [docs/theory/15_fast_cvt.md](docs/theory/15_fast_cvt.md).

```python
from cvtsim.io.params import load
from cvtsim.fast.model import build_fast, simulate_fast, FCH
res = simulate_fast(build_fast(load(tune="cvtech_stock"), scenario="hill_climb"), 8.0)
```

## Events

Run a set-up through a whole event on the fast model:

```
python -m cvtsim.events run cvtech_stock --event acceleration      # 150 ft, standing start
python -m cvtsim.events run cvtech_stock --event hill_climb
python -m cvtsim.events run cvtech_stock --event endurance --laps 106   # about 4 h, runs in about 30 s
python -m cvtsim.events compare default cvtech_stock cvtech_0651_1100 --event endurance
```

A course is a list of segments (length, grade, corner radius, traction) in
`cvtsim/io/tracks/`. The driver follows a target speed set by corner grip. There are no brakes
yet, so the driver slows only by lifting off. The hill-climb and endurance courses are invented
and representative; replace them with your own course. Each run reports:
* time;
* speed;
* engine rpm band;
* belt slip;
* ratio range;
* wheel spin;
* how much faster than the corner limit the car arrived.

See [docs/theory/16_events.md](docs/theory/16_events.md).

Loads and durability over a whole event, with the worst moments replayed in the full model, and
sweeps over set-ups and final drives:

```
python -m cvtsim.event_loads run cvtech_stock --event endurance --laps 106 --hours 4 --replays 3
python -m cvtsim.event_loads sweep default cvtech_stock cvtech_0651_1100 --event endurance --final-drive 8.1 8.96
```

The fast model's event loads are mean-torque loads. The replays add firing pulses and belt
dynamics at the moments where the loads peak. See
[docs/theory/17_event_loads.md](docs/theory/17_event_loads.md).

## Your car

The DR26 car is the example. Copy its parameter file, replace the numbers, and run with
`--car my_car.yaml`. [docs/YOUR_CAR.md](docs/YOUR_CAR.md) lists what to measure, in order of how
much it changes the results. To list every input that is still assumed:

```
python -m cvtsim.inputs my_car.yaml --open
```

## How far to trust the results

* **Launch time and top speed** match published Baja data within the measured spread
  ([REPORT.md](REPORT.md)).
* **Trends** are reliable: which part or tune raises or lowers a load, and when the peak occurs.
* **Absolute stress and life** are only as good as the shaft geometry, materials, inertias and
  bearing ratings you enter. `python -m cvtsim.inputs --open` shows which are still assumed.
* **Uncertainty:** `python -m cvtsim.uq screen` ranks the inputs by effect, and
  `python -m cvtsim.uq bands` gives result ranges over the assumed input ranges.

## Tests

```
pytest -q -m "not slow"     # about 6 min
pytest -q                   # everything, including the coupled-system checks, about 15 min
```

## License

MIT. See [LICENSE](LICENSE).
