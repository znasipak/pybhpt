# `pybhpt.radial` performance

Timing of the homogeneous radial Teukolsky solutions. See the
[Performance overview](index) for the methodology and [reference machine](reference-machine).

`RadialTeukolsky(s, l, m, a, omega, r).solve(method)` computes the `In`/`Up` homogeneous
solutions on the supplied radial grid. Three things drive the cost, benchmarked
separately: the **solver method**, the **radial grid length**, and the **mode frequency**
`omega` (more oscillatory solutions cost more, differently for each method). All numbers
below use a representative mode `(s=-2, l=5, m=2, a=0.9)`.

## Method comparison

Median solve time by method, at `omega = 0.3` on a 256-point grid:

| method | median time |
|---|---|
| `HBL` | 0.53 ms |
| `AUTO` | 0.55 ms |
| `GSN` | 4.8 ms |
| `MST` | 649 ms |

`HBL` and `GSN` are numerical ODE integrators; `MST` is an analytic series solution in the
mode frequency. At this mode `MST` is roughly **1200× slower** than `HBL` — its series
requires many terms away from the low-frequency/near-superradiant regime where it is
typically preferred for accuracy. `AUTO` tracks `HBL` almost exactly here, confirming it
correctly avoids the expensive analytic method for a generic mode.

Full data: [`benchmarks/data/radial_method.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_method.csv).

## Grid-length scaling

Median solve time (method `AUTO`) versus the number of radial output points:

| nsamples | 64 | 256 | 1024 | 4096 |
|---|---|---|---|---|
| median time | 0.42 ms | 0.52 ms | 1.0 ms | 3.4 ms |

Cost grows sub-linearly at small grids (dominated by a fixed setup cost — boundary
condition solve, method dispatch) and approaches linear in `nsamples` at large grids,
where evaluating the solution on the output grid dominates.

Full data: [`benchmarks/data/radial_grid.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_grid.csv).

## Frequency sweep

Median solve time versus `omega`, for each method, on a fixed 256-point grid:

| method | omega=0.02 | 0.1 | 0.6 | 2.5 |
|---|---|---|---|---|
| `HBL` | 0.48 ms | 0.46 ms | 0.58 ms | 0.88 ms |
| `AUTO` | 0.50 ms | 0.49 ms | 0.60 ms | 0.93 ms |
| `GSN` | 5.1 ms | 4.7 ms | 5.5 ms | 15.7 ms |
| `MST` | 187 ms | 213 ms | 851 ms | 827 ms |

The ODE methods (`HBL`, `GSN`) grow modestly with `omega` — more oscillatory solutions
need more integration steps, roughly 2× over this range. `MST` grows much more sharply
(~4.5× from `omega=0.02` to `omega=0.6`) since higher frequency widens the series needed
for convergence; it is consistently the most expensive method across the whole sweep for
this mode.

Full data: [`benchmarks/data/radial_freq.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_freq.csv).
