# Regression report — `feature-continuous-solution` vs `main`

> **Status: 28k-mode pilot run. Rerun pending.** The results below are from the initial
> 28,512-case sweep. The grid in `benchmarks/regression_sweep.py` has since been expanded
> to ~137,280 cases (retrograde x, e=0.9, L up to 30, m+k ∈ {−10,−2,0,2,10}, n up to 50);
> this report will be regenerated from that run via `benchmarks/regression/run_regression.sh`.

**Verdict: no amplitude regression.** Across 28,512 mode/orbit/resolution combinations
(spin s = −2, 0, +2; broad a, p, e, x; three resolutions), converged Teukolsky/scalar
amplitudes on this branch match `main` to solver precision. The only differences are
(1) grossly under-resolved modes, which are sampling-dependent on both branches and
converge to agreement at adequate resolution, and (2) *reported precision*, which this
branch deliberately changed (error-estimate rework). Neither is a correctness regression.

| | |
|---|---|
| branch | `feature-continuous-solution` @ `64fde6b` |
| baseline | `main` @ `1ccff1e` |
| machine | Apple arm64, Darwin 24.6.0, Python 3.12.13, numpy 2.2.0, BLAS Accelerate |
| modes compared | 28,512 (0 NaN either branch) |

## Method

For every case the two builds (each with its own matched Python wrapper + compiled
extension) computed `amplitude('In'/'Up')` and `precision('In'/'Up')` via the API common
to both: `KerrGeodesic(a,p,e,x,ns)` → `TeukolskyMode(s,l,m,k,n,orbit).solve(orbit)`.
Amplitudes are compared by relative difference `|A−B| / max(|A|,|B|)`; precision by ratio
`prec_branch / prec_main`. Parameter grid (validity-aware — invalid/below-separatrix
orbits skipped):

- **spins** s ∈ {−2, 0, +2}; **L** ∈ {2, 3, 5, 8}; **m** a spread in [−L, L]; **k** ∈ {−2, 0, 2}; **n** ∈ {−3, 0, 3, 10}
- **orbits** a ∈ {0, 0.5, 0.9, 0.99} × (p,e) ∈ {(10,0), (10,0.5), (8,0.3)} × x ∈ {1, 0.5, 0.1}
- **resolutions** nsamples ∈ {2⁶, 2⁹, 2¹¹}
- equatorial (x=1) restricts to k=0; circular (e=0) restricts to n=0.

## Amplitude agreement (both \|amp\| > 1e-8)

| side | nsamples | n | median | p90 | p99 | max |
|---|---|---|---|---|---|---|
| In | 64 | 4725 | 1.3e-14 | 1.0e-11 | 2.0e-7 | 8.5e-6 |
| In | 512 | 4725 | 1.2e-14 | 1.9e-13 | 6.3e-11 | 2.0e-8 |
| In | 2048 | 4725 | 1.7e-14 | 2.5e-13 | 7.4e-11 | 2.1e-8 |
| Up | 64 | 6505 | 6.9e-14 | 2.0e-7 | 7.7e-2 | **1.7** |
| Up | 512 | 6505 | 5.0e-14 | 5.2e-11 | 8.2e-8 | 3.5e-5 |
| Up | 2048 | 6505 | 5.0e-14 | 5.1e-11 | 8.3e-8 | 3.7e-5 |

**Converged (nsamples ≥ 512): median ~1e-14 (In) / 5e-14 (Up), 99th ~1e-10 / 1e-7.**

## Interpretation of the differences

1. **Under-resolved modes (nsamples = 64).** The `Up` tail at ns=64 (up to rel ≈ 1.7) is
   entirely high radial-harmonic modes (n = 10) on a 64-point grid — well below the ~2|n|+2
   samples needed. Both branches are *differently inaccurate* here because the geodesic
   sampling changed (the branch defaults to the Darwin phase parametrization). At ns ≥ 512
   these vanish (max rel-diff 3.5e-5). Not a regression.

2. **One near-extremal converged outlier.** The single ns=2048 mode above 1e-5 is
   `s=2, l=3, m=−3, n=10, a=0.99` at rel = 3.7e-5 — within its own reported precision
   (~3e-5), i.e. at the accuracy floor of a near-extremal high-n mode, not a code error.

3. **Reported precision changed by design.** The branch's error-estimate rework
   (RMS/L2 condition number, Neumaier compensated summation, geometric truncation
   estimate) changes `precision(...)`: median ratio branch/main ≈ 1.2 (In) / 1.6 (Up)
   — branch is mostly slightly more conservative — with a wide two-sided tail (the largest
   ratios are parity-forbidden / near-zero modes where `main` reports machine epsilon).
   This is intended and independent of the amplitude values.

*(Near-zero modes — parity-forbidden or physically negligible amplitudes ≲ 1e-15 — show
large relative differences purely from the sign of round-off and are excluded above by the
\|amp\| > 1e-8 cut; both branches return effective zeros there.)*

## Reproduce

```bash
# build each branch (main in a worktree; symlink extern/boost)
python benchmarks/regression_sweep.py --module <branch-so> --pybhpt <branch-repo> --dump branch.npz
python benchmarks/regression_sweep.py --module <main-so>   --pybhpt <main-repo>   --dump main.npz
python benchmarks/regression_sweep.py --compare branch.npz main.npz
```

Per-side / per-resolution summary statistics are in [`summary.csv`](summary.csv). The raw
per-mode dumps (~5 MB each) are regenerable with the commands above and are not committed.
