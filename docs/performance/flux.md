# `pybhpt.flux` performance

Timing of the gravitational-wave flux computation. See the [Performance overview](index) for the methodology and
[reference machine](reference-machine).

`FluxMode(geo, teuk)` turns an **already-solved** `TeukolskyMode` into the energy, angular
momentum, and Carter-constant fluxes at infinity and on the horizon. It runs no solver of
its own: it reads the mode amplitudes and forms a fixed algebraic combination of them. So
there is no useful flux-specific timing table — the cost of a flux is the cost of the
Teukolsky modes that feed it.

## Flux is a rounding error on the mode solve

Measured against the solve it consumes, at the same mode (`s = -2`, `a = 0.9`):

| orbit class | mode `(l,m,k,n)` | nsamples | solve | flux | flux share |
|---|---|---|---|---|---|
| circular–equatorial | (2,2,0,0) | 64 | 0.30 ms | 2 µs | 0.52% |
| circular–equatorial | (2,2,0,0) | 4096 | 7.8 ms | 5 µs | 0.06% |
| spherical | (2,2,0,0) | 64 | 0.30 ms | 1 µs | 0.43% |
| generic | (5,3,-2,3) | 512 | 3.6 ms | 2 µs | 0.04% |
| generic | (8,4,2,10) | 4096 | 53.5 ms | 11 µs | 0.02% |

The largest flux share anywhere in the swept grid is **0.52%**, at the cheapest mode on the
coarsest grid; it drops to a few hundredths of a percent for the expensive modes, because
the flux step is essentially constant in `nsamples` while the solve is not.

Full data: [`benchmarks/data/teuk_flux.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/teuk_flux.csv),
produced by [`bench_teuk.py`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/bench_teuk.py)
alongside the mode-solve benchmarks.

## Budgeting a flux calculation

Plan from the [`pybhpt.teuk` numbers](teuk) instead, and treat the flux itself as free:

```
wall time ≈ (number of (l, m, k, n) modes summed) × (per-mode solve cost)
```

with the per-mode cost read off the teuk aggregate table for your orbit class and
resolution. Two consequences from that table:

- **Orbit class dominates.** A generic (eccentric + inclined) orbit costs ~7× a
  circular–equatorial one at `nsamples = 512`, because its source integral is a full 2D
  radial × polar loop.
- **The mode count, not the flux algebra, sets the bill.** Cost per mode keeps climbing
  with `l` and `n` rather than plateauing, so how deep the mode sum is truncated matters
  far more than anything on this page.

The one-off `KerrGeodesic` construction ([`pybhpt.geo` performance](geo)) is shared by
every mode in the sum and is negligible against it, except for near-polar orbits where the
geodesic solve itself costs ~10 ms.
