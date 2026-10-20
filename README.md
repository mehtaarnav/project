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
