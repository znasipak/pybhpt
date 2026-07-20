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
| circular–equatorial | (2,2,0,0) | 0.29 ms | 0.50 ms | 6.8 ms |
| eccentric–equatorial | (2,2,0,0) | 0.38 ms | 0.69 ms | 8.3 ms |
| spherical | (2,2,0,0) | 0.28 ms | 0.53 ms | 8.0 ms |
| generic | (2,2,0,0) | 0.42 ms | 3.4 ms | 29.7 ms |
| generic | (5,3,-2,3) | 0.40 ms | 3.4 ms | 30.1 ms |
| generic | (8,4,2,10) | 0.57 ms | 3.6 ms | 53.1 ms |

`s = 0` and `s = +2` follow the same pattern at comparable cost to `s = -2` (see the full
data). The dominant driver is **orbit class**, not mode number: generic (eccentric +
inclined) orbits need a full 2D radial × polar source integral and are already ~5-7×
slower than the equatorial/spherical classes at `nsamples = 512`, growing to a much larger
gap at high resolution — see the stage breakdown below for why.

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
| circular–equatorial | 0.57 ms | 0.22 ms | 0.33 ms | 0.015 ms (**3%**) |
| eccentric–equatorial | 0.67 ms | 0.23 ms | 0.41 ms | 0.034 ms (**5%**) |
| spherical | 0.55 ms | 0.20 ms | 0.32 ms | 0.029 ms (**5%**) |
| generic | 3.40 ms | 0.22 ms | 0.40 ms | 2.77 ms (**82%**) |

For equatorial/spherical orbits, `swsh` and `radial` together account for essentially all
of the solve time — source integration is cheap because it only needs a 1D loop. For
generic orbits, source integration dominates (~82% of the total), consistent with the
optimization work this project has focused on for the |s|=2 generic path: the source
integrand is evaluated on the full 2D radial × polar grid, and that 2D loop — not the
spheroidal harmonic or the radial solve — is where a generic-orbit mode actually spends
its time.

Full data: [`benchmarks/data/teuk_stages.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/teuk_stages.csv).
