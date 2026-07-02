# `pybhpt.swsh` performance

Timing of the spin-weighted spheroidal harmonic. See the [Performance overview](index) for
the methodology and [reference machine](reference-machine).

`SpinWeightedHarmonic(s, l, m, gamma, theta)` does two things at construction, which scale
with different parameters and are benchmarked separately:

1. a **spectral solve** for the eigenvalue and the spherical–spheroidal coupling
   coefficients, whose cost is driven by the spheroidicity `gamma = a·omega` (it widens
   the coupling band and the eigenproblem);
2. **evaluation on the `theta` grid** of `S`, `S'`, `S''`, whose cost is linear in the
   number of grid points on top of the fixed solve.

## Spectral solve

Median time to build the eigenvalue + coupling representation (measured on a tiny grid so
the solve dominates):

| gamma | median solve time |
|---|---|
| 0 (spherical) | 3–15 µs |
| 1 | 0.2–0.8 ms |
| 4 | 0.35–0.6 ms |
| 10 | 0.7–24 ms |

At `gamma = 0` the harmonic reduces to a spin-weighted *spherical* harmonic and no
eigenvalue solve is needed, so it is ~100× cheaper. For nonzero `gamma` the cost is set by
the coupling bandwidth rather than `l` alone, and grows sharply for large `gamma` (strong
spheroidicity, e.g. high-frequency modes of a near-extremal black hole), where the
truncated coupling matrix must be enlarged.

Full data: [`benchmarks/data/swsh_solve.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_solve.csv).

## On-grid evaluation

Median total construction time (spectral solve **plus** on-grid evaluation) versus grid
length, for three representative modes:

| mode `(s, l, m, gamma)` | ns = 64 | 512 | 4096 |
|---|---|---|---|
| `(-2, 2, 2, 2)` | 0.39 ms | 0.76 ms | 3.7 ms |
| `(-2, 8, 4, 4)` | 0.77 ms | 1.4 ms | 6.2 ms |
| `(2, 13, 7, 8)` | 9.5 ms | 10.6 ms | 21.7 ms |

The slope with `nsamples` is the on-grid evaluation cost; the intercept is the fixed
per-mode spectral-solve offset (visible as the near-constant time at small grids for the
high-`gamma` mode). Once solved, an instance can be re-evaluated at arbitrary angles via
`.eval(theta, deriv=…)` without repeating the solve.

Full data: [`benchmarks/data/swsh_grid.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_grid.csv).
