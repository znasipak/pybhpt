# `pybhpt.geo` performance

Timing of `KerrGeodesic` construction across orbit geometry and grid resolution. See the
[Performance overview](index) for the methodology and [reference machine](reference-machine).

`KerrGeodesic(a, p, e, x, nsamples)` is solved once per orbit and then reused by every
mode, so its cost is amortized over an entire `(l, m, k, n)` sweep. The cost is set mainly
by the orbit geometry and, at large grids, grows roughly linearly in `nsamples`.

Representative median times (`a = 0.9`):

| orbit class | nsamples = 64 | 512 | 4096 |
|---|---|---|---|
| circular–equatorial | 7 µs | 9 µs | 28 µs |
| eccentric–equatorial | 0.03–0.04 ms | 0.09–0.13 ms | 0.6–0.9 ms |
| spherical (inclined circular) | 0.2–11 ms | 0.3–12 ms | 1–15 ms |
| generic (eccentric + inclined) | 0.3–12 ms | 0.4–12 ms | 1.7–17 ms |

- Circular–equatorial orbits are essentially free — the frequencies and trajectory are
  closed-form, so there is no grid-dependent solve.
- Eccentric–equatorial orbits scale cleanly with `nsamples`.
- Inclined orbits (spherical and generic) carry a large fixed cost that is nearly
  independent of `nsamples`; the wide range reflects the inclination `x` (near-polar
  orbits are the expensive end). Because this cost is paid once per orbit and shared by
  all modes, it is rarely the bottleneck of a production run.

```{figure} ../_static/figures/geo_hist.png
:alt: Distribution of KerrGeodesic construction times by orbit class
:width: 100%

Construction time across every sampled `(a, p, e, x, nsamples)`, split by orbit class.
The classes separate by orders of magnitude, and the two inclined classes are visibly
bimodal — the upper cluster is the near-polar sampling.
```

The bimodality is inclination, not resolution: at `x = 0.5` the spherical and generic
orbits sit at 0.29 and 0.40 ms median, while at `x = 0.1` both jump to ~11.7 ms, a factor
of ~30 with everything else held equal. That split is what widens the ranges in the table
above, so read those as two populations rather than a spread.

Full data: [`benchmarks/data/geodesic.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/geodesic.csv).
