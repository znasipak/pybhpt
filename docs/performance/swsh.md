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
| 0 (spherical) | 4–15 µs |
| 1 | 0.2–1.3 ms |
| 4 | 0.35–0.87 ms |
| 10 | 0.74–23 ms |

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
| `(-2, 2, 2, 2)` | 0.40 ms | 0.74 ms | 3.7 ms |
| `(-2, 8, 4, 4)` | 0.74 ms | 1.3 ms | 6.1 ms |
| `(2, 13, 7, 8)` | 9.0 ms | 10.8 ms | 21.4 ms |

The slope with `nsamples` is the on-grid evaluation cost; the intercept is the fixed
per-mode spectral-solve offset (visible as the near-constant time at small grids for the
high-`gamma` mode). Once solved, an instance can be re-evaluated at arbitrary angles via
`.eval(theta, deriv=…)` without repeating the solve.

Full data: [`benchmarks/data/swsh_grid.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_grid.csv).

## Random parameter-space sample

The same treatment as the [radial random sample](radial), with one continuous axis instead
of two: the radial solver takes spin and frequency separately, while the harmonic sees only
their product `gamma = a·omega`. 600 random points — `l` uniform on [2, 50], one random `m`
per point in `[-l, l]`, `gamma` uniform on [0.02, 5] — at `s = -2` on a fixed 256-point
`theta` grid, timing the full construction (spectral solve plus on-grid evaluation). Every
point succeeded. Median 1.07 ms, p90 2.13 ms, p99 3.21 ms.

```{figure} ../_static/figures/swsh_random_l_vs_gamma.png
:alt: Construction time over mode number l and spheroidicity gamma
:width: 100%

Cost rises with both, `l` more steeply than `gamma` — but the hottest points do not sit in
the top-right corner where a smooth scaling would put them.
```

```{figure} ../_static/figures/swsh_random_l_vs_m.png
:alt: Construction time over mode number l and azimuthal number m
:width: 100%

The expensive modes cluster near `m = 0`, not at `|m| = l`. This is the opposite of the
radial solver, where large `|m|` is the costly end.
```

```{figure} ../_static/figures/swsh_random_m_vs_gamma.png
:alt: Construction time over azimuthal number m and spheroidicity gamma
:width: 100%

The `m ≈ 0` band darkens as `gamma` grows: low `|m|` at high spheroidicity is the
expensive combination.
```

Binning by `|m|/l` shows the effect cleanly, at fixed everything else:

| `abs(m)/l` | n | median | p90 |
|---|---|---|---|
| 0.0–0.2 | 116 | 1.49 ms | 3.01 ms |
| 0.2–0.5 | 157 | 1.38 ms | 2.23 ms |
| 0.5–0.8 | 184 | 0.98 ms | 1.72 ms |
| 0.8–1.0 | 143 | 0.84 ms | 1.28 ms |

A log-log fit (R² = 0.75) gives `log l` +0.39, `log gamma` +0.25, and `|m|/l` **−0.65** —
the only strongly negative coefficient anywhere in these benchmarks. Nearly axisymmetric
modes need the widest spherical–spheroidal coupling band, so they carry the largest
eigenproblem.

```{warning}
**Construction time is not a smooth function of `gamma`.** The slowest point in the sample
(`l = 15, m = -1, gamma = 4.8123044703`) takes ~29 ms, roughly 20× its neighbours, and it
reproduces on re-timing. Perturbing `gamma` in its tenth significant digit — to
`4.8123044751` — drops it back to 1.5 ms.

The **results are unaffected**: both values return the same eigenvalue, the same coupling
range, and harmonics agreeing to 2.8e-10, with `couplingstatus = 0` and no warning. The
extra time is wasted work, not a loss of accuracy.

The cost is entirely in the spectral solve (identical spike at a 4-point `theta` grid) and
comes from the truncation-growth loop in `swsh.cpp`: `nmax` is raised by 10 and the dense
eigenproblem re-solved until a convergence test passes. That test compares a coupling
coefficient across successive truncations against `SPECTRAL_COUPLING_CONVERGE_EPS = 1e-25`,
evaluated at the first index where the coefficient has already fallen to 1e-25 of the peak
— i.e. at the floating-point noise floor, where passing or failing is decided by the low
bits of the input. When it keeps failing, `nmax` climbs toward the limit of 600 and each
retry costs O(nmax³).

Incidence is low: 3 of 600 points here, and 1 of 800 in a separate random probe, exceed 10×
the median. Budget from the median and p90 columns and expect the occasional mode to cost
an order of magnitude more than its neighbours.
```

Full data: [`benchmarks/data/swsh_random.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_random.csv).
