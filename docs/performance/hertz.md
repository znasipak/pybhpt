# `pybhpt.hertz` performance

Timing of the Hertz-potential solve. See the [Performance overview](index) for the methodology and
[reference machine](reference-machine).

```{note}
No timing table on this page. `HertzMode(teuk, gauge).solve()` is a transform of an
*already-solved* `TeukolskyMode`, so what a user waits for is the mode solve itself —
see [`pybhpt.teuk` performance](teuk). The benchmark harness covers `geo`, `radial`,
`swsh`, and `teuk` (with flux); `hertz` is deliberately not measured separately.
```
