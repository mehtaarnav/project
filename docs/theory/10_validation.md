# 10 — Validation data

Measured data that simulator outputs are compared against. Each entry names its source, what
it measures, and how far it transfers to the DR26 car.

## SAE 2010-01-0861 (Martinez et al., Universidad Simón Bolívar, 2008 car)

Verified from the full text.

| Quantity | Value | Transfer to DR26 |
|---|---|---|
| 100 ft (30.48 m) acceleration, 12 kit configurations | 4.10 to 4.54 s (Table 4) | Plausibility band only. Different engine (B&S Intek 205432, 10 hp gross, governed 3800 rpm), unknown vehicle mass and final drive. |
| CVT ratio range (engine/secondary speed) | 4.15 low, 0.737 high (paper: 0.241 to 1.357 driven/drive) | CVTech hardware limit; DR26 uses an imposed 3.0 to 0.8 window. |
| CVT re-stabilization after a load step | 0.97 s (Fig. 7) | Order-of-magnitude check of shift dynamics (sheave masses, damping). |
| Peak output power through the CVT | about 7 to 8 hp at output speed 3000-4000 rpm (Figs. 8-11) | Efficiency check: ratio of output to engine power. |
| Heavier flyweights | earlier engagement, less power below about 2500 output rpm, less slip above 3500 | Qualitative trend the tune comparison must reproduce. |
| Stiffer primary spring | same direction as lighter flyweights | Qualitative trend. |
| Stiffer secondary spring | holds low ratio longer, earlier slip at high ratio | Qualitative trend. |

## How the simulator is checked against it

* `performance` metrics report the 100 ft time. A DR26 baseline outside roughly 3.5 to 5.5 s
  points to a modeling or parameter error, given the differences above.
* Tune comparisons (heavier vs lighter flyweights, spring swaps) must reproduce the
  qualitative trends in the table. A reversed trend is a defect.
