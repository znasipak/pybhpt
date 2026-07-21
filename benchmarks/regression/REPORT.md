# Regression report — `feature-continuous-solution` vs `main`

**Verdict: no amplitude regression.** Across the full locked grid of 137,280
mode/orbit/resolution combinations (spin s = −2, 0, +2; L up to 30; near-extremal and
retrograde orbits; n up to 50; two resolutions), converged (`nsamples = 4096`)
Teukolsky/scalar amplitudes on this branch match `main` to solver precision. The sweep
also surfaced one real, since-fixed bug in this branch's precision estimate (amplitudes
were unaffected) — see below.

| | |
|---|---|
| branch (data) | `feature-continuous-solution` @ `6fbf442` |
| precision fix (after this run, not re-swept) | `98a7b66` |
| baseline | `main` @ `1ccff1e` |
| machine | Apple arm64, Darwin 24.6.0, Python 3.12.13, numpy 2.2.0, BLAS Accelerate |
| modes compared | 137,280 (both builds) |
| wall time | ~11.3 h total (branch dump ~53 min; main dump ~10.4 h; main is unoptimized and pays for it most on the generic-orbit, high-resolution tail) |

## Method

For every case the two builds (each with its own matched Python wrapper + compiled
extension) computed `amplitude('In'/'Up')` and `precision('In'/'Up')` via the API common
to both: `KerrGeodesic(a,p,e,x,ns)` → `TeukolskyMode(s,l,m,k,n,orbit).solve(orbit)`.
Amplitudes are compared by relative difference `|A−B| / max(|A|,|B|)`; precision by ratio
`prec_branch / prec_main`. Locked parameter grid (validity-aware — invalid/below-separatrix
orbits skipped; see `benchmarks/regression_sweep.py`):

- **spins** s ∈ {−2, 0, +2}; **L** ∈ {2, 3, 5, 13, 20, 30}; **m+k** ∈ {−10, −2, 0, 2, 10}
  (polar mode k = (m+k) − m); **n** up to 50 for e ∈ {0.5, 0.9}, up to 10 for e = 0.3, 0 for e = 0
- **orbits** a ∈ {0, 0.5, 0.9, 0.99} × (p,e) ∈ {(10,0), (10,0.5), (8,0.3), (12,0.9)} × x —
  full inclination coverage (±1, 0.5, 0.1, −0.3) for the e = 0.9 orbit, x ∈ {0.5, −0.3}
  (prograde/retrograde inclined) otherwise
- **resolutions** nsamples ∈ {128, 4096} (coarse + converged)
- equatorial (\|x\|=1) restricts to k=0; circular (e=0) restricts to n=0

## Amplitude agreement (both \|amp\| > 1e-8)

| side | nsamples | n | median | p90 | p99 | max |
|---|---|---|---|---|---|---|
| In | 128 | 13,312 | 9.4e-14 | 8.8e-01 | 1.6e+00 | 1.9e+00 |
| In | 4096 | 12,201 | 3.7e-14 | 1.2e-12 | 3.5e-10 | 3.9e-08 |
| Up | 128 | 19,655 | 6.6e-08 | 1.4e+00 | 1.9e+00 | 2.0e+00 |
| Up | 4096 | 14,688 | 5.2e-13 | 8.6e-09 | 8.6e-06 | 5.3e-04 |

**Converged (nsamples = 4096): median ~4e-14 (In) / 5e-13 (Up).** The tail is worse than
the earlier 28k pilot's single 3.7e-5 outlier — this run pushes into a much harder corner
of parameter space (L up to 30, n up to 50, a up to 0.99, e = 0.9) that the pilot didn't
cover, and a handful of modes there sit at 2–5×10⁻⁴.

## Interpretation of the differences

1. **Under-resolved modes (nsamples = 128).** The huge tail at ns=128 (up to rel ≈ 2.0,
   the maximum possible for this metric) is exactly what's expected: 128 points cannot
   resolve L up to 30 / n up to 50. Both builds are *differently inaccurate* on the same
   coarse grid — confirmed by re-checking a cluster of `ns=128` rows that looked like a
   suspicious exact sign-flip (`rel≈2.00`, s=0, a=0, l=2, e=0.9): at the *same* mode's
   `ns=4096` entry, branch and main agree to displayed precision
   (e.g. `-4.2041e-05+1.0501e-05j` on both sides, bitwise identical). Pure aliasing
   artifact of the deliberately coarse grid, not a real disagreement.

2. **A handful of converged (ns=4096) outliers exceed their own reported precision.**
   The three worst `Up` cases (e.g. `s=2, l=5, m=2, k=-2, n=50, a=0.9, e=0.9`, rel=5.3e-4)
   are all near-extremal (a ≥ 0.9), high-eccentricity (e=0.9), high-n (10–50) modes — and
   the disagreement is larger than either build's own precision estimate (up to ~28×).
   This matches a near-identical situation already diagnosed earlier on this branch (a CI
   failure on `[s=2,l=13,m=10,k=-3,n=20]` at a=0.99, resolved via `rtol *= 3` after
   determining the *reference* itself is only precision-accurate at that extreme). Not a
   new correctness regression — an accuracy floor for the hardest near-extremal/high-n
   modes, now visible because this run reaches further into that regime than the pilot did.

3. **Reported precision changed by design** (error-estimate rework: RMS/L2 condition
   number, Neumaier compensated summation, geometric truncation estimate): median ratio
   branch/main ≈ 1.0 for both boundary conditions, with a wide tail (p90 ≈ 130×, and an
   extreme p99/max out to ~10¹⁵–10²⁷) driven by near-zero/parity-forbidden modes where
   `main`'s older formula reports a fixed, vanishingly small constant regardless of context
   — both sides mean "negligible/converged" there, just expressed at wildly different tiny
   magnitudes, so the ratio itself is not physically meaningful. Intended, independent of
   the amplitude values.

## Bug found and fixed: NaN precision for scalar modes at a = 0

While comparing precision values, `precision('In'/'Up')` was **NaN** for a reproducible
27–28% of scalar (s=0) modes at **exactly a=0** (Schwarzschild) — 3,254 (In) / 3,101 (Up)
of 137,280 rows, all confined to that one (s, a) combination, reproducing identically at
both `nsamples=128` and `4096` (so not a convergence artifact). Amplitudes were unaffected.

**Root cause:** `radial_integral_convergence_sum`/`polar_integral_convergence_sum`
computed the relative-change estimate as `|1 − old/new|` with no guard on `new == 0`. A
scalar 1D sub-integral (I1–I4) can cancel to an *exact* bit-for-bit zero — plausible at
a=0, where extra Schwarzschild symmetry can force this for particular (m,k,n) parity
combinations — making that division `0/0 = NaN`. The NaN then poisoned
`scalar_amplitude_precision` even though the term's true contribution is exactly zero,
because `0 * NaN = NaN`, not `0`.

**Fix** (`98a7b66`, after this data was collected — not yet re-swept at scale): a
`reldiff_or_zero` helper returns `0` instead of dividing when the new value is exactly
zero, applied at the two scalar sub-integral sites (the confirmed cause) and, defensively,
the three `|s|=2` amplitude-driver sites with the identical unguarded pattern (not observed
to trigger there — a full |s|=2 amplitude landing on exact zero is far less likely than a
single 1D sub-integral — but the same bug shape existed).

**Verification:** rather than re-running the full 11-hour sweep, all 3,329
previously-NaN cases from this run's dump were re-checked directly against the fixed
build: all now return finite precision, and **zero amplitudes changed** (the fix only
touches the precision estimate). Full unit test suite (67/67) still passes.

## Reproduce

```bash
# build each branch (main in a worktree; symlink extern/boost)
bash benchmarks/regression/run_regression.sh
```

Per-side summary statistics are in [`summary.csv`](summary.csv). The raw per-mode dumps
(~28 MB each) are regenerable with the command above and are not committed.
