# Viewer and plot data contract

The simulator, the viewer and the dashboards communicate only through this contract.

## Simulator → Python

`cvtsim.cvt.system.simulate(model, t_end, nout, n_frames=400)` returns `(res, snaps)`.

* `res`: `(n_out, len(CH))` array; channel names are in `CH` (09 §6 lists the additions).
* `snaps`: dict with
  * `t` `(F,)` frame times [s]
  * `belt` `(F, N, 2)` node positions [m], world frame (primary axis at origin, `+X` toward secondary)
  * `tension` `(F, N)` element tension [N] (element k joins nodes k and k+1)
  * `x_p, x_s, theta_e, theta_s, r_w` `(F,)` at the same times

## Python → viewer (`viewer.json` in the run directory)

Written by `cvtsim.viz.export.write_viewer(run_dir, res, CH, snaps, p, sh=None)`.

```json
{
  "meta":     {"run": "...", "commit": "...", "dt": 2.5e-7, "tune": "default"},
  "geometry": {"C": 0.232, "beta_deg": 12, "rp_min": 0.0315, "rp_mech_max": 0.09144,
               "rs_mech_min": 0.04242, "rs_mech_max": 0.09703, "belt_width": 0.03,
               "belt_thickness": 0.014, "gap": 0.002, "x_p_stop": 0.0178, "s_s": 0.02448,
               "n_weights": 3, "r_w0": 0.04445, "ramp": [[u, r_w], ...],
               "helix": [[x, h], ...], "r_cam": 0.0381,
               "spring_p": {"k": 13208, "F0": 683.04},
               "spring_s": {"k": 3353.66, "F0": 381.6, "k_t": 3.93, "theta_t0_deg": 30}},
  "frames":   {"t": [...], "belt": [[[x, y], ...], ...], "tension": [[...], ...],
               "x_p": [...], "x_s": [...], "theta_e": [...], "theta_s": [...], "r_w": [...]},
  "series":   {"t": [...], "<channel>": [...], ...},
  "shafts":   {"t": [...], "<bearing or stress channel>": [...]},
  "mechanism": {"spring_p_compression": [...], "spring_s_compression": [...],
                "torsion_twist": [...], "helix_twist": [...], "helix_angle_deg": [...],
                "ramp_slope": [...]},
  "events":   {"engagement": t, "peak_belt_tension": t, "shift_complete": t}
}
```

`mechanism` has one value per frame (same length as `frames.t`):

* Spring compressions are measured from free length, in m. Free lengths are not in the parameters,
  so draw each spring at an assumed free length minus its compression.
* `torsion_twist` is the torsion spring's total twist in rad: pretwist plus helix twist.
* `helix_twist` is the movable half's rotation relative to the fixed half, in rad.
* `helix_angle_deg` is the local helix angle under the roller.
* `ramp_slope` is `dr_w/du` at the current primary travel.

`events` holds the times in s that are worth bookmarking. A key is absent if the event did not
occur. `mechanism` and `events` are written only for full simulation output.

Floats are rounded to 5 significant figures. `series` is decimated to at most 2000 points.

## Viewer (`viewer/index.html`)

A static page: no build step, libraries from the cdnjs or jsdelivr CDN with pinned versions.
Open it with `?run=<path to viewer.json>`, or drop a file on it.

* **Transport:** play, pause, scrub, speed (0.1× to 4×), frame-accurate time readout.
* **Side view, belt plane:** both sheaves drawn at their current pitch radius, outer rims and hub.
  The belt polyline is coloured by tension with a perceptual gradient (viridis or turbo) and a
  colour bar. Sheave spokes rotate with `theta_e` and `theta_s`. Flyweight positions are drawn
  at radius `r_w`, at 120° spacing on the primary.
* **Cut sections, axial view:** for the primary, the fixed and movable halves (trapezoids at `±β`),
  the movable half displaced by `x_p`, the flyweight on its ramp, and the compression spring drawn
  at its current length. For the secondary, the halves displaced by `x_s`, the helix roller at
  `h(x)`, and the linear and torsion springs. Belt cross-section drawn in the groove at the
  current radius.
* **Synced plots:** engine rpm vs vehicle speed, ratio, sheave travel, belt tension, slip, and
  traction margin, each with a cursor at the current time; clicking a plot seeks to that time.
* **Readouts:** rpm, speed, ratio, clamp forces, belt temperature, bearing loads at the cursor.
