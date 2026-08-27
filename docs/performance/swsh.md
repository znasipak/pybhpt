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
| 0 (spherical) | 3.6–15 µs |
| 1 | 62–112 µs |
| 4 | 73–250 µs |
| 10 | 0.23–0.45 ms |

At `gamma = 0` the harmonic reduces to a spin-weighted *spherical* harmonic and no
eigenvalue solve is needed, so it is ~50× cheaper. For nonzero `gamma` the cost is set by
the coupling bandwidth rather than `l` alone, and rises with `gamma` (strong
spheroidicity, e.g. high-frequency modes of a near-extremal black hole), where the
truncated coupling matrix must be enlarged. The growth is smooth: the widest spread within
a `gamma` column is a factor of ~3.4, across all `(s, l)` in the sweep.

Full data: [`benchmarks/data/swsh_solve.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_solve.csv).

## On-grid evaluation

Median total construction time (spectral solve **plus** on-grid evaluation) versus grid
length, for three representative modes:

| mode `(s, l, m, gamma)` | ns = 64 | 512 | 4096 |
|---|---|---|---|
| `(-2, 2, 2, 2)` | 0.18 ms | 0.56 ms | 3.8 ms |
| `(-2, 8, 4, 4)` | 0.37 ms | 0.99 ms | 6.4 ms |
| `(2, 13, 7, 8)` | 0.48 ms | 1.5 ms | 9.3 ms |

The slope with `nsamples` is the on-grid evaluation cost; the intercept is the fixed
per-mode spectral-solve offset, a few hundred µs even for the high-`gamma` mode, so from
`ns = 512` up the grid evaluation dominates. Once solved, an instance can be re-evaluated
at arbitrary angles via `.eval(theta, deriv=…)` without repeating the solve.

Full data: [`benchmarks/data/swsh_grid.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_grid.csv).

## Random parameter-space sample

The same treatment as the [radial random sample](radial), with one continuous axis instead
of two: the radial solver takes spin and frequency separately, while the harmonic sees only
their product `gamma = a·omega`. 600 random points — `l` uniform on [2, 50], one random `m`
per point in `[-l, l]`, `gamma` uniform on [0.02, 5] — at `s = -2` on a fixed 256-point
`theta` grid, timing the full construction (spectral solve plus on-grid evaluation). Every
point succeeded. Median 0.65 ms, p90 1.18 ms, p99 1.53 ms, slowest point 1.60 ms.

```{figure} ../_static/figures/swsh_random_l_vs_gamma.png
:alt: Construction time over mode number l and spheroidicity gamma
:width: 100%

Cost rises with both, `l` far more steeply than `gamma`; the hottest points are the
high-`l` end, at no particular `gamma`.
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
| 0.0–0.2 | 116 | 0.84 ms | 1.45 ms |
| 0.2–0.5 | 157 | 0.81 ms | 1.24 ms |
| 0.5–0.8 | 184 | 0.62 ms | 1.07 ms |
| 0.8–1.0 | 143 | 0.51 ms | 0.86 ms |

A log-log fit (R² = 0.69) gives `log l` +0.46, `log gamma` +0.15, and `|m|/l` **−0.46** —
the only strongly negative coefficient anywhere in these benchmarks. Nearly axisymmetric
modes need the widest spherical–spheroidal coupling band, so they carry the largest
eigenproblem.

```{note}
**Construction time used to be a discontinuous function of `gamma`.** Before the
truncation-test fix, the slowest point in this sample (`l = 15, m = -1, gamma = 4.8123044703`) took ~29 ms,
roughly 20× its neighbours, while perturbing `gamma` in its tenth significant digit — to
`4.8123044751` — dropped it back to 1.5 ms. The results were unaffected; the extra time was
wasted work.

The cause was the truncation-growth loop in `swsh.cpp`, which raises `nmax` by 10 and
re-solves the dense eigenproblem until a convergence test passes. That test took the
*relative* change of a single coupling coefficient across successive truncations, evaluated
at the first index where the coefficient had already fallen to 1e-25 of the peak — a ratio
of two round-off-level numbers, so passing or failing was decided by the low bits of the
input, and a failure sent `nmax` climbing toward its limit of 600 at O(nmax³) per retry.

The test now measures the largest change in the coupling vector against the peak
coefficient, requires the tail to have decayed below 1e-14 of the peak inside the current
basis, and tolerates 1e-14 (scaled by `gamma²` above `gamma = 1`). Two round-off-level tails
therefore compare as converged instead of at random. The same point now takes 0.42 ms with
an unchanged eigenvalue and coupling range, no sample point exceeds 10× the median, and the
solve is 3–20× faster across the whole `gamma` sweep. Eigenvalues and coupling coefficients
agree with the old truncation to ~1e-13 relative over `|gamma| <= 10`.
```

Full data: [`benchmarks/data/swsh_random.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/swsh_random.csv).
