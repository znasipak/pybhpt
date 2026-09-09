# `pybhpt.radial` performance

Timing of the homogeneous radial Teukolsky solutions. See the
[Performance overview](index) for the methodology and [reference machine](reference-machine).

`RadialTeukolsky(s, l, m, a, omega, r).solve(method)` computes the `In`/`Up` homogeneous
solutions on the supplied radial grid. Four things drive the cost, benchmarked
separately: the **solver method**, the **radial grid length**, the **mode frequency**
`omega` (more oscillatory solutions cost more, differently for each method), and the
**mode numbers and spin** `(l, m, a)`. The first three use a representative mode
`(s=-2, l=5, m=2, a=0.9)`; the last sweeps `l`, `m`, and `a` directly.

## Method comparison

Median solve time by method, at `omega = 0.3` on a 256-point grid:

| method | median time |
|---|---|
| `HBL` | 0.54 ms |
| `AUTO` | 0.54 ms |
| `GSN` | 4.9 ms |

`HBL` and `GSN` are numerical ODE integrators, and `AUTO` dispatches between them. `GSN`
is ~10× the cost of `HBL` at this mode; `AUTO` tracks `HBL` almost exactly, confirming it
picks the cheap integrator for a generic mode.

```{note}
The analytic `MST` series solver is not benchmarked here. It was removed from `AUTO`'s
dispatch in `6d214d5` as unreliable at high `l`, so `AUTO` will never select it. It
remains callable explicitly via `solve("MST")`, but it is not part of the supported
solve path these numbers describe.
```

Full data: [`benchmarks/data/radial_method.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_method.csv).

## Grid-length scaling

Median solve time (method `AUTO`) versus the number of radial output points:

| nsamples | 64 | 256 | 1024 | 4096 |
|---|---|---|---|---|
| median time | 0.43 ms | 0.54 ms | 1.1 ms | 3.5 ms |

Cost grows sub-linearly at small grids (dominated by a fixed setup cost — boundary
condition solve, method dispatch) and approaches linear in `nsamples` at large grids,
where evaluating the solution on the output grid dominates.

Full data: [`benchmarks/data/radial_grid.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_grid.csv).

## Frequency sweep

Median solve time versus `omega`, for each method, on a fixed 256-point grid:

| method | omega=0.02 | 0.1 | 0.6 | 2.5 |
|---|---|---|---|---|
| `HBL` | 0.49 ms | 0.45 ms | 0.60 ms | 0.87 ms |
| `AUTO` | 0.49 ms | 0.42 ms | 0.59 ms | 0.92 ms |
| `GSN` | 5.1 ms | 4.9 ms | 5.7 ms | 16.0 ms |

Both integrators grow with `omega` — more oscillatory solutions need more integration
steps — but not at the same rate: `HBL` costs 1.8× more at `omega = 2.5` than at `0.02`,
while `GSN` costs 3.1× more, so the gap between them widens with frequency. `AUTO` tracks
`HBL` across the whole sweep.

Full data: [`benchmarks/data/radial_freq.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_freq.csv).

## Mode number and spin

Method `AUTO`, `omega = 0.3`, 256-point grid — so the mode `(l, m)` and the background spin
`a` are the only things varying. Each cell is the median over five equally spaced `m` in
`[-l, l]`:

| l | a = 0 | a = 0.3 | a = 0.7 | a = 0.9 | a = 0.99 |
|---|---|---|---|---|---|
| 2 | 0.38 ms | 0.42 ms | 0.43 ms | 0.46 ms | 0.48 ms |
| 5 | 0.43 ms | 0.48 ms | 0.49 ms | 0.54 ms | 0.61 ms |
| 10 | 0.48 ms | 0.55 ms | 0.57 ms | 0.64 ms | 0.95 ms |
| 15 | 0.57 ms | 0.64 ms | 0.69 ms | 0.77 ms | 1.32 ms |
| 30 | 0.96 ms | 1.06 ms | 1.16 ms | 1.37 ms | 2.35 ms |
| 50 | 1.87 ms | 2.05 ms | 2.16 ms | 2.47 ms | 3.94 ms |

All 150 `(l, m, a)` combinations solved successfully — nothing failed, even at
`l = 50, a = 0.99`.

The two parameters compound rather than acting independently. Growth from `l=2` to `l=50`
is **4.9×** at `a = 0`, but **8.1×** at `a = 0.99`; equivalently, going from `a = 0` to
`a = 0.99` costs 1.3× at `l = 2` but 2.1× at `l = 50`. Spin is nearly free for slowly
rotating black holes and only becomes a real cost driver in the near-extremal, high-`l`
corner.

`m` matters mostly through its sign at high spin. At `a = 0.9, l = 30` the spread across
`m` is mild (1.22 ms at `m = 0` to 2.01 ms at `m = -30`), but at `a = 0.99, l = 15` the
retrograde end costs 4.3× the `m = 0` mode (3.33 ms vs 0.77 ms). Prograde modes stay
cheaper than retrograde ones of the same `|m|` throughout.

```{note}
**One reproducible spike.** `l = 50, m = -25, a = 0.99` takes ~18 ms — over 3× its
immediate neighbours (`m = -24`: 5.8 ms, `m = -26`: 4.6 ms), breaking the otherwise smooth
trend. It is not a method-dispatch artifact: timing the integrators directly gives
`HBL` = 18.4 ms and `AUTO` = 17.9 ms at that mode versus 5.5 ms for both at `m = -24`, so
the extra cost is inside the `HBL` integrator itself, not in `AUTO` falling back to
something slower. Left in the data rather than smoothed away.
```

Full data: [`benchmarks/data/radial_mode.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_mode.csv).

## Random parameter-space sample

The sweeps above each hold two of the three cost drivers fixed, so none of them shows how
`l`, `a` and `omega` interact. This one samples all of them at once: 600 random points with
`l` uniform on [2, 50], one random `m` per point in `[-l, l]`, `a` log-uniform in `1 - a`
over [0.01, 1], and `omega` uniform on [0.02, 5]. Method `AUTO`, 256-point grid. Every
point solved — nothing failed anywhere in that box. Times run 0.42–4.81 ms, median 1.61 ms.

Spin is plotted as `1 - a` on a log axis, near-extremal to the right: the cost structure
bunches against `a = 1`, and a linear spin axis compresses the entire `a = 0.9` to `0.99`
range — where most of the variation lives — into the last tenth of the plot.

```{figure} ../_static/figures/radial_random_l_vs_omega.png
:alt: Solve time over mode number l and frequency omega
:width: 100%

Cost climbs along both axes with no visible interaction between them; the cheap corner is
low `l`, low `omega`.
```

```{figure} ../_static/figures/radial_random_l_vs_a.png
:alt: Solve time over mode number l and black hole spin
:width: 100%

The gradient is almost entirely vertical: `l` dominates, and spin only sharpens the tail in
the near-extremal column on the right.
```

```{figure} ../_static/figures/radial_random_omega_vs_a.png
:alt: Solve time over frequency omega and black hole spin
:width: 100%

With `l` averaged out, the remaining structure is a mild climb in `omega` and an elevated
near-extremal edge.
```

A log-log fit of `log(time)` against the four parameters (R² = 0.80) ranks them:

| term | exponent |
|---|---|
| `log l` | +0.40 |
| `log omega` | +0.20 |
| `−log(1 − a)` | +0.09 |
| `m/l` | −0.02 |

So cost scales roughly as `l^0.4 ω^0.2`, with `l` the leading driver by about 2× over
frequency and 4× over spin. Rank correlations tell the same story: `l` +0.75, `omega`
+0.51, `1 - a` −0.29.

The `m/l` coefficient is near zero, but raw `|m|` correlates at +0.54 — that is `|m|`
riding along with `l`, not an independent effect, consistent with the `l` versus `l - |m|`
comparison above. The sign asymmetry does survive and is spin-dependent: near-extremal
(`1 - a < 0.05`) retrograde modes run 2.04 ms against 1.83 ms prograde, narrowing to
1.37 vs 1.27 ms at slow spin.

```{note}
These exponents describe the box sampled here, not a universal scaling law. The earlier
draft of this sample drew `omega` log-uniformly, which concentrated points at low frequency
where the `omega` cost is flat, and reported `log omega` at +0.07 instead of +0.20 — the
same physics weighted differently. Read them as a guide to relative cost within this range.
```

Full data: [`benchmarks/data/radial_random.csv`](https://github.com/znasipak/pybhpt/blob/main/benchmarks/data/radial_random.csv).
