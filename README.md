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
