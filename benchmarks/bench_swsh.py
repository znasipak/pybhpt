"""Benchmark: spin-weighted spheroidal harmonic, split into two costs.

The grid-based SpinWeightedHarmonic does two things at construction: (A) a spectral
solve for the eigenvalue + spherical-spheroidal coupling coefficients, whose cost grows
with the spheroidal mode l and the spheroidicity gamma = a*omega (both widen the coupling
band), and (B) evaluation of S, S', S'' on the supplied theta grid, whose cost is linear
in the number of grid points on top of the fixed solve. The two scale with different
parameters, so they are measured separately:

  * swsh_solve.csv -- vary (s, l, gamma) on a tiny theta grid; dominated by the solve.
  * swsh_grid.csv  -- fix representative modes, vary nsamples; the slope over nsamples is
                      the on-grid eval cost, above the fixed per-mode solve offset.
  * swsh_random.csv -- a random sample of the (l, m, gamma) parameter space, behind the
                      Performance page's scatter plots, mirroring radial_random.csv. The
                      radial solver takes spin and frequency separately; the harmonic sees
                      only their product, so this sweep has one continuous axis (gamma)
                      where the radial one has two. l is uniform, m uniform in [-l, l],
                      gamma uniform, and the theta grid is held fixed so the sample varies
                      the mode alone. Seeded, so the same points come back on a re-run.

Emits benchmarks/data/swsh_solve.csv, swsh_grid.csv, and swsh_random.csv.
"""
import numpy as np

from .timing import bench, save
from .grids import FIELD_SPINS, RESOLUTIONS

SWSH_L = (2, 3, 5, 8, 13, 20)          # spheroidal mode number
SWSH_GAMMA = (0.0, 1.0, 4.0, 10.0)     # spheroidicity a*omega (0 -> spherical)
SOLVE_NS = 4                           # tiny grid so the solve dominates axis A

# representative (s, l, m, gamma) for the on-grid eval scaling (axis B)
GRID_MODES = [(-2, 2, 2, 2.0), (-2, 8, 4, 4.0), (2, 13, 7, 8.0)]

# random sample (axis C) of the joint (l, m, gamma) space
RANDOM_N = 600
RANDOM_SEED = 20260827
RANDOM_L_RANGE = (2, 50)      # uniform integer, matching radial_random.csv
RANDOM_GAMMA = (0.02, 5.0)    # uniform spheroidicity a*omega
RANDOM_NS = 256               # fixed theta grid


def _theta(n):
    # interior grid, away from the poles where the |s| harmonics vanish
    return np.linspace(0.05, np.pi - 0.05, n)


def run(outdir, quick=False):
    from pybhpt.swsh import SpinWeightedHarmonic

    ls = (2, 5, 13) if quick else SWSH_L
    gammas = (0.0, 4.0) if quick else SWSH_GAMMA
    resolutions = (2**6, 2**9, 2**12) if quick else RESOLUTIONS

    # --- axis A: spectral solve (eigenvalue + coupling) ---
    th_solve = _theta(SOLVE_NS)
    rowsA = []
    for s in FIELD_SPINS:
        for l in ls:
            if l < max(abs(s), 1):
                continue
            m = l
            for g in gammas:
                t = bench(lambda s=s, l=l, m=m, g=g:
                          SpinWeightedHarmonic(s, l, m, g, th_solve))
                rowsA.append({"s": s, "l": l, "m": m, "gamma": g,
                              "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("swsh_solve", rowsA, outdir)

    # --- axis B: on-grid eval scaling with nsamples ---
    modes = GRID_MODES[:2] if quick else GRID_MODES
    rowsB = []
    for (s, l, m, g) in modes:
        for ns in resolutions:
            th = _theta(ns)
            t = bench(lambda s=s, l=l, m=m, g=g, th=th:
                      SpinWeightedHarmonic(s, l, m, g, th))
            rowsB.append({"s": s, "l": l, "m": m, "gamma": g, "nsamples": ns,
                          "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("swsh_grid", rowsB, outdir)

    # --- axis C: random sample of (l, m, gamma) at a fixed theta grid ---
    rng = np.random.default_rng(RANDOM_SEED)
    n_samples = 60 if quick else RANDOM_N
    th_rand = _theta(RANDOM_NS)
    rowsC = []
    for _ in range(n_samples):
        l = int(rng.integers(RANDOM_L_RANGE[0], RANDOM_L_RANGE[1] + 1))
        m = int(rng.integers(-l, l + 1))
        g = float(rng.uniform(*RANDOM_GAMMA))
        row = {"s": -2, "l": l, "m": m, "gamma": g, "nsamples": RANDOM_NS, "status": "ok"}
        try:
            t = bench(lambda l=l, m=m, g=g: SpinWeightedHarmonic(-2, l, m, g, th_rand))
            row.update(min=t["min"], median=t["median"], p90=t["p90"])
        except Exception as exc:
            row.update(status=f"error: {type(exc).__name__}", min="", median="", p90="")
        rowsC.append(row)
    return save("swsh_random", rowsC, outdir)
