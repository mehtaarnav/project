# cvtsim

Rubber-belt CVT dynamics and FE simulator for the DR26 Baja car. It couples a segmented belt
with analytic cone sheaves, centrifugal and torque-sensing clutches, the engine and the vehicle,
then computes Timoshenko shaft loads and stresses and runs a 10-node tetrahedral FE solver.
No commercial solver is used.

* Build specification: [BRIEF.md](BRIEF.md)
* Theory, sources and test targets: [docs/theory/](docs/theory/README.md)
* Results and validation: [REPORT.md](REPORT.md)

```
pip install -e .[fast,test]
pytest -q
python examples/launch.py
```

## Choosing parts

```
python -m cvtsim.app        # open http://localhost:8765
```

Pick each spring, flyweight set, ramp, helix and surface from a dropdown, then start a launch.
Results appear in the run table with links to the dashboard and the mechanism viewer. To add a
part, add an entry to `cvtsim/io/tuning/components.yaml` and it appears in the dropdown.
