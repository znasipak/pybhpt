# `pybhpt.redshift` performance

Timing of the Detweiler redshift invariant. See the [Performance overview](index) for the methodology and
[reference machine](reference-machine).

```{note}
No timing table on this page. `RedshiftCoefficients` carries a real per-call cost of its
own (unlike `flux`/`hertz`, it is not post-processing of a solved mode), but it is
outside the current benchmark scope — the harness covers `geo`, `radial`, `swsh`, and
`teuk` (with flux). Treat the [`pybhpt.teuk` numbers](teuk) as a lower bound on any
pipeline that also computes the redshift invariant.
```
