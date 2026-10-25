# Convergence of peak loads

Peak loads come from belt engagement and lock-up, in the first 2.5 s of a launch. This study
checks that they do not depend on three numerical settings:
* the time step;
* the friction smoothing speed `belt.v_eps`, which sets how sharply the belt switches between
  stick and slip;
* the output sampling rate.

To repeat the check at the current default, run `python examples/convergence_peaks.py`: five 2.5 s runs, about 30 min on 4 cores. The table below was recorded around the old default (0.05 m/s).

## Result

| Setting | Effect on peaks | Status |
|---|---|---|
| Time step 1.25e-7 s | halving it changes every peak by under 0.6 % | converged |
| Output at 4 kHz | within 0.4 % of 40 kHz for every load peak; the jackshaft Goodman factor reads 3 % high; the brief belt compression is under-read | adequate; use `--nout` for more |
| `v_eps` 0.05 m/s (old default) | halving it raised belt tension 5 % and PTO side load 4 %; lowered the jackshaft Goodman factor 3 % | not converged |
| `v_eps` 0.0125 m/s (new default) | halving it changes peaks by under 1.5 %; halving dt by under 0.6 % | converged |

Changing the default from 0.05 to 0.0125 m/s raises the reported peaks by roughly:
* belt tension 6 %;
* PTO side load 5 %;
* jackshaft stress 3 %;
* bearing load 1.4 %.

The jackshaft Goodman factor falls from 1.03 to 0.99. Runs recorded before this change used
0.05 m/s.

The model briefly puts the belt in compression (about −260 N) during engagement chatter, at
every setting. A real belt carries little compression, so treat that value with care.

## Data (first 2.5 s of the default launch, 40 kHz output; percentages against the first column)

| metric | base 40 kHz | base 4 kHz | dt/2 | v_eps x2 | v_eps x0.5 | v_eps x0.25 | v_eps x0.25, dt/2 |
|---|---|---|---|---|---|---|---|
| T_pto_peak | 253.7 | 252.7 (-0.4%) | 253.7 (-0.0%) | 255.6 (+0.8%) | 251.2 (-1.0%) | 256.1 (+1.0%) | 256.2 (+1.0%) |
| T_tight_max | 1425 | 1424 (-0.1%) | 1427 (+0.1%) | 1290 (-9.5%) | 1498 (+5.1%) | 1510 (+5.9%) | 1508 (+5.8%) |
| T_tight_min | -260.1 | -190.5 (+26.8%) | -259.9 (+0.1%) | -261.6 (-0.6%) | -261.2 (-0.4%) | -263.3 (-1.2%) | -263.3 (-1.2%) |
| hub_peak | 1325 | 1322 (-0.3%) | 1329 (+0.3%) | 1192 (-10.1%) | 1383 (+4.4%) | 1391 (+5.0%) | 1391 (+4.9%) |
| hub_peak_lp200 | 1184 | 1183 (-0.1%) | 1185 (+0.0%) | 1092 (-7.8%) | 1221 (+3.2%) | 1246 (+5.3%) | 1247 (+5.3%) |
| bearing_peak | 3238 | 3229 (-0.3%) | 3239 (+0.0%) | 2950 (-8.9%) | 3279 (+1.3%) | 3285 (+1.4%) | 3284 (+1.4%) |
| bearing_peak_lp200 | 2942 | 2942 (+0.0%) | 2942 (-0.0%) | 2779 (-5.5%) | 2943 (+0.0%) | 2900 (-1.4%) | 2899 (-1.5%) |
| vm_sec_nom | 165.4 | 164.9 (-0.3%) | 165.4 (+0.0%) | 154.7 (-6.5%) | 168.3 (+1.8%) | 170.3 (+3.0%) | 171.1 (+3.5%) |
| vm_prim_nom | 138.7 | 138.3 (-0.3%) | 138.7 (+0.0%) | 139.7 (+0.7%) | 137.1 (-1.1%) | 140.6 (+1.4%) | 140.7 (+1.5%) |
| n_good_prim | 1.434 | 1.439 (+0.4%) | 1.434 (-0.0%) | 1.43 (-0.3%) | 1.426 (-0.6%) | 1.406 (-2.0%) | 1.405 (-2.0%) |
| n_good_sec | 1.032 | 1.063 (+3.0%) | 1.032 (-0.0%) | 1.08 (+4.6%) | 1.002 (-3.0%) | 0.9877 (-4.3%) | 0.9869 (-4.4%) |
| t_bearing_peak | 1.916 | 1.916 (+0.0%) | 1.916 (-0.0%) | 1.894 (-1.1%) | 2.269 (+18.4%) | 2.315 (+20.8%) | 2.315 (+20.8%) |
| engage_rpm | 1914 | 1911 (-0.2%) | 1914 (+0.0%) | 1914 (+0.0%) | 1914 (-0.0%) | 1914 (-0.0%) | 1914 (+0.0%) |
| v_end | 8.79 | 8.79 (-0.0%) | 8.79 (-0.0%) | 8.787 (-0.0%) | 8.801 (+0.1%) | 8.812 (+0.2%) | 8.812 (+0.2%) |
| x_s_end | 9.074 | 9.073 (-0.0%) | 9.072 (-0.0%) | 9.323 (+2.8%) | 9.138 (+0.7%) | 9.172 (+1.1%) | 9.178 (+1.2%) |
