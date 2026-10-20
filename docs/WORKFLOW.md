# How work is done in cvtsim

This is the working standard for the simulator. It treats engineering simulation the way a
good ML team treats training: every result is a tracked run with a recorded config, and no
result counts until it has passed validation. Nothing changes the physics without a spec,
and nothing is trusted without a check that could have failed.

## 1. The loop

```
question ──► spec ──► targets ──► code + tests ──► tracked runs ──► review ──► report
   ▲                                                                             │
   └──────────────────────────── findings feed the next question ◄───────────────┘
```

| Step | Artifact | Gate to pass before moving on |
|---|---|---|
| 1. Question | one sentence in the commit or issue: what decision does this inform? | It names a decision (a spring rate, a shaft size), not just "explore". |
| 2. Spec | `docs/theory/NN_*.md`: equations with derivation, symbols and units, frame conventions, assumptions with bias direction, sources | Every equation is derived or cited; unverifiable sources are marked UNVERIFIED. |
| 3. Targets | `docs/theory/targets.py`: closed-form numbers computed **without** importing `cvtsim` | Targets are independent of the code they test. |
| 4. Code + tests | `cvtsim/<module>`, `tests/test_<module>.py` with the spec's target IDs | `pytest -q` green; tolerances exactly as specified. |
| 5. Runs | `runs/<stamp>-<name>/` via `cvtsim.runs.Run` | The validation set in §3 passes. |
| 6. Review | checklist in §5 | Each item ticked by someone other than the author (or an independent agent). |
| 7. Report | `REPORT.md`, `runs/INDEX.md` | Every number in the report traces to a run directory or a test ID. |

## 2. Runs are experiments

A run is to this project what a training run is to a model.

| ML training | cvtsim |
|---|---|
| hyperparameters | `config.yaml`: resolved parameters with source notes, plus overrides |
| code version | `meta.json`: git commit and dirty flag (a dirty run is not citable) |
| training curves | `series.csv.gz` and `dashboard.png` |
| eval metrics | `metrics.json` |
| validation set | convergence and energy checks (§3) |
| ablation / sweep | `group: sweep` runs, one change per run against a baseline |
| leaderboard | `runs/INDEX.md` (`python -m cvtsim.runs`) |

Rules:
* **Baseline first.** Every sweep includes an unmodified baseline in the same window, at the same commit.
* **Change one thing per run.** If two parameters move, it is two runs.
* **Never overwrite a run.** Runs are append-only; a wrong run is superseded, not edited.
* **Commit before you cite.** Results quoted in `REPORT.md` come from runs with `dirty: false`.

## 3. Validation: when a run may be believed

A run is valid only if all of these hold. They play the role a held-out set plays in training:
they catch a model that is fitting the numerics instead of the physics.

| Check | Criterion | Why |
|---|---|---|
| Energy closure | `max |ΔE − (W_in − W_diss)| / W_in < 1e-4` | Any missing force, wrong sign or unstable integration breaks the power balance. |
| Time-step convergence | outputs at `dt` and `dt/2` agree to 1e-3 | Results must not depend on the integrator. |
| Penalty convergence | trend reported for `k_n × {0.5, 2, 4}` | `k_n` is a physical belt property *and* a penalty; report sensitivity, don't hide it. |
| Closed-form anchor | at least one quantity checked against a hand calculation (e.g. torque capacity `2μR F_a / cos β` when slipping) | Confirms the mechanism, not just the numbers. |
| Placeholders listed | every `PLACEHOLDER` the result depends on is named next to it | Illustrative numbers must never be read as predictions. |

## 4. Numbers

* Every parameter in `cvtsim/io/dr26_default.yaml` carries a `source:`. Allowed tags: BRIEF,
  measured (with file or figure reference), DERIVED (with spec section), ASSUMED, UNVERIFIED,
  PLACEHOLDER.
* No number is invented. It is derived, cited, measured, or flagged.
* When data disagree (e.g. secondary stroke 24.48 mm vs 23.22 mm implied by measured radii),
  record both and the one used; do not average them away.

## 5. Review checklist

Run this on every change that touches physics.

- [ ] **Signs:** forces are "forces on" the named body; reactions are explicit.
- [ ] **Units:** SI inside the code; rpm only at I/O.
- [ ] **Frames:** world `+X` primary→secondary, `+Z` along shafts, rotation positive about `+Z`.
- [ ] **Energy:** new forces either derive from a potential the energy sum includes, or their power is in `W_diss`.
- [ ] **Tests:** no tolerance loosened, no target changed without a written derivation in the spec, no skip or xfail.
- [ ] **Spec drift:** the docs describe what the code does (formula, symbols, defaults).
- [ ] **Runs:** validation set (§3) passes; numbers in the report link to run directories.

If a test fails, the first question is whether the *spec* is wrong. Fixing a spec is allowed
and normal; it requires a derivation in the spec and a note in the commit. Loosening a
tolerance to make a failure disappear is never allowed.

## 6. Commits

* Small, single-purpose, green (`pytest -q`) at every commit.
* Message: imperative summary line, then what changed in the physics or the spec and why.
* Generated outputs are committed only from clean runs that the report cites.

## 7. Definition of done

A piece of work is done when the decision it was meant to inform can be made from the report
alone: the number, its uncertainty (convergence and sensitivity), the placeholders it depends
on, and a link to the run that produced it.
