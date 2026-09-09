# `pybhpt.teuk` performance

Timing of the inhomogeneous Teukolsky amplitude solve. See the
[Performance overview](index) for the methodology and [reference machine](reference-machine).

`TeukolskyMode(s, l, m, k, n, geo).solve(geo)` is the end-to-end amplitude computation:
build the spheroidal harmonic, solve the radial equation, and integrate the source over the
particle's trajectory. It is a single C++ call, so it is benchmarked two ways: the
**aggregate** cost users actually pay per mode, and a **per-stage breakdown** that isolates
how much of that cost is the spheroidal harmonic, the radial solve, and everything else
(mainly source integration).

## Aggregate cost

Median solve time (`s = -2`) at `a = 0.9`, across orbit class, mode number, and resolution:

| orbit class | mode `(l,m,k,n)` | nsamples = 64 | 512 | 4096 |
|---|---|---|---|---|
| circular–equatorial | (2,2,0,0) | 0.23 ms | 0.43 ms | 6.7 ms |
| eccentric–equatorial | (2,2,0,0) | 0.22 ms | 0.54 ms | 7.3 ms |
| spherical | (2,2,0,0) | 0.21 ms | 0.57 ms | 6.6 ms |
| generic | (2,2,0,0) | 0.30 ms | 3.4 ms | 29.2 ms |
| generic | (5,3,-2,3) | 0.35 ms | 3.4 ms | 30.4 ms |
| generic | (8,4,2,10) | 0.41 ms | 3.6 ms | 53.5 ms |
| generic | (13,7,-5,20) | 0.50 ms | 3.8 ms | 99.6 ms |
| generic | (20,10,5,30) | 0.55 ms | 4.1 ms | 99.3 ms |
| generic | (30,15,-8,50) | 0.78 ms | 4.2 ms | 192 ms |

`s = 0` and `s = +2` follow the same pattern at comparable cost to `s = -2` (see the full
data). The dominant driver is **orbit class**, not mode number: generic (eccentric +
inclined) orbits need a full 2D radial × polar source integral and are already ~8×
slower than the equatorial/spherical classes at `nsamples = 512`, growing to a much larger
gap at high resolution — see the stage breakdown below for why.

The last three rows are high-`(l,n)` spot checks (`SPOTCHECK_MODES` in `bench_teuk.py`),
added to bound the cost of the large-mode tail that appears in bigger sweeps (e.g. the
regression grid in `benchmarks/regression_sweep.py`, which goes up to `l=30`, `n=50`) but
isn't otherwise exercised by the curated modes above. Cost keeps climbing with `l`/`n`
rather than plateauing — `nsamples=4096` goes from 29 ms at `l=2` to 192 ms at `l=30, n=50`,
~6.6× over the tested range.

Full data (all spins, all resolutions): [`benchmarks/data/teuk_aggregate.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/teuk_aggregate.csv).

## Stage breakdown

`TeukolskyMode.solve()` does not expose its internal stages for direct timing — the
Python/cython API accepts precomputed `teuk`/`swsh` objects for reuse, but the current
binding does not actually consume them (it always calls one monolithic C++ solve). To
still get a genuine breakdown, each mode below is solved once to read off the *actual*
radial and polar grids and frequency it used, and then a standalone `SpinWeightedHarmonic`
and `RadialTeukolsky` are timed independently **on those same grids and parameters** — real
measurements, not a cheap proxy (see [`bench_teuk.py`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/bench_teuk.py)
for why that distinction matters). The remainder, `full − swsh − radial`, is a *derived*
estimate of the source-integration-plus-orchestration cost, not a direct measurement.

One representative mode per orbit class (`s = -2, l = 5, m = 3`), `nsamples = 512`:

| orbit class | full | swsh | radial | remainder (source integration) |
|---|---|---|---|---|
| circular–equatorial | 0.57 ms | 0.18 ms | 0.31 ms | 0.09 ms (**15%**) |
| eccentric–equatorial | 0.69 ms | 0.16 ms | 0.43 ms | 0.10 ms (**15%**) |
| spherical | 0.57 ms | 0.16 ms | 0.32 ms | 0.09 ms (**16%**) |
| generic | 3.35 ms | 0.16 ms | 0.43 ms | 2.76 ms (**82%**) |

For equatorial/spherical orbits, `swsh` and `radial` together account for ~85% of the solve
time — source integration is cheap because it only needs a 1D loop. Treat those three
remainder values as order-of-magnitude only: they are a small difference of two larger
measurements, and repeat runs put them anywhere between ~8% and ~20%. For generic orbits,
source integration dominates (~82% of the total), consistent with the
optimization work this project has focused on for the |s|=2 generic path: the source
integrand is evaluated on the full 2D radial × polar grid, and that 2D loop — not the
spheroidal harmonic or the radial solve — is where a generic-orbit mode actually spends
its time.

Full data: [`benchmarks/data/teuk_stages.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/teuk_stages.csv).

## Flux post-processing

`FluxMode(geo, teuk)` consumes an *already-solved* `TeukolskyMode`, so it is measured here
against the solve it depends on rather than on its own page. Median times, `s = -2`:

| orbit class | mode `(l,m,k,n)` | nsamples | solve | flux | flux share |
|---|---|---|---|---|---|
| circular–equatorial | (2,2,0,0) | 64 | 0.24 ms | 1 µs | 0.61% |
| circular–equatorial | (2,2,0,0) | 4096 | 6.5 ms | 4 µs | 0.06% |
| spherical | (2,2,0,0) | 64 | 0.23 ms | 1 µs | 0.55% |
| generic | (5,3,-2,3) | 512 | 3.6 ms | 1 µs | 0.04% |
| generic | (8,4,2,10) | 4096 | 54.5 ms | 4 µs | 0.01% |

Across every orbit class, mode, and resolution tested, the flux step is **at most 0.61%**
of the mode solve, and its share *falls* with resolution — the flux is a fixed algebraic
combination of the mode amplitudes, while the solve grows with `nsamples`. This is why
[`pybhpt.flux` performance](flux) reports the mode-solve cost instead of a flux-specific
timing table.

Full data: [`benchmarks/data/teuk_flux.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/teuk_flux.csv).
