# `pybhpt.geo` performance

Timing of `KerrGeodesic` construction across orbit geometry and grid resolution. See the
[Performance overview](index) for the methodology and [reference machine](reference-machine).

`KerrGeodesic(a, p, e, x, nsamples)` is solved once per orbit and then reused by every
mode, so its cost is amortized over an entire `(l, m, k, n)` sweep. The cost is set mainly
by the orbit geometry and, at large grids, grows roughly linearly in `nsamples`.

Representative median times (`a = 0.9`):

| orbit class | nsamples = 64 | 512 | 4096 |
|---|---|---|---|
| circular–equatorial | 6 µs | 9 µs | 26 µs |
| eccentric–equatorial | 0.03–0.04 ms | 0.09–0.13 ms | 0.6–0.8 ms |
| spherical (inclined circular) | 0.2–11 ms | 0.3–11 ms | 1–14 ms |
| generic (eccentric + inclined) | 0.2–11 ms | 0.4–12 ms | 1.6–15 ms |

- Circular–equatorial orbits are essentially free — the frequencies and trajectory are
  closed-form, so there is no grid-dependent solve.
- Eccentric–equatorial orbits scale cleanly with `nsamples`.
- Inclined orbits (spherical and generic) carry a large fixed cost that is nearly
  independent of `nsamples`; the wide range reflects the inclination `x` (near-polar
  orbits are the expensive end). Because this cost is paid once per orbit and shared by
  all modes, it is rarely the bottleneck of a production run.

Full data: [`benchmarks/data/geodesic.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/geodesic.csv).
