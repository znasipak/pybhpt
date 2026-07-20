"""Benchmark: homogeneous radial Teukolsky solutions (RadialTeukolsky.solve).

Three independent axes, each isolating one cost driver:

  * radial_method.csv -- solve time by SolutionMethod (AUTO/MST/HBL/GSN) at a fixed
                          representative mode and radial grid. MST is analytic (series in
                          the frequency); HBL/GSN/TEUK are numerical ODE integrators. AUTO
                          picks whichever is expected to be cheapest/most stable.
  * radial_grid.csv    -- solve time vs radial grid length (nsamples), method AUTO, at a
                          fixed mode. Exposes the per-point cost of evaluating solutions
                          on the output grid on top of any fixed setup cost.
  * radial_freq.csv    -- solve time vs mode frequency omega, for each method, at a fixed
                          grid length. Higher |omega| means more oscillatory solutions:
                          MST needs more series terms and the ODE integrators take more
                          steps, so this exposes how each method scales with frequency.

Emits benchmarks/data/radial_method.csv, radial_grid.csv, radial_freq.csv.
"""
import numpy as np

from .timing import bench, save
from .grids import RESOLUTIONS

METHODS = ("AUTO", "MST", "HBL", "GSN")

# representative mode for the method-comparison and grid-scaling axes
REP_S, REP_L, REP_M, REP_A, REP_OMEGA = -2, 5, 2, 0.9, 0.3

OMEGAS = (0.02, 0.05, 0.1, 0.3, 0.6, 1.2, 2.5)
FREQ_NS = 256  # fixed grid length for the frequency sweep


def _r_grid(a, nsamples, rmax=50.0):
    # stay safely outside the horizon r_+ = 1 + sqrt(1 - a^2)
    rplus = 1.0 + np.sqrt(1.0 - a * a)
    rmin = rplus * 1.05
    return np.linspace(rmin, rmax, nsamples)


def run(outdir, quick=False):
    from pybhpt.radial import RadialTeukolsky

    methods = ("AUTO", "HBL") if quick else METHODS
    resolutions = (2**6, 2**9, 2**12) if quick else RESOLUTIONS
    omegas = (0.05, 0.6, 2.5) if quick else OMEGAS

    # --- axis A: method comparison at a fixed mode/grid ---
    r = _r_grid(REP_A, 256)
    rowsA = []
    for method in methods:
        t = bench(lambda method=method:
                  RadialTeukolsky(REP_S, REP_L, REP_M, REP_A, REP_OMEGA, r).solve(method))
        rowsA.append({"method": method, "s": REP_S, "l": REP_L, "m": REP_M,
                      "a": REP_A, "omega": REP_OMEGA, "nsamples": r.shape[0],
                      "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("radial_method", rowsA, outdir)

    # --- axis B: grid-length scaling, method AUTO ---
    rowsB = []
    for ns in resolutions:
        r = _r_grid(REP_A, ns)
        t = bench(lambda r=r:
                  RadialTeukolsky(REP_S, REP_L, REP_M, REP_A, REP_OMEGA, r).solve("AUTO"))
        rowsB.append({"nsamples": ns, "s": REP_S, "l": REP_L, "m": REP_M,
                      "a": REP_A, "omega": REP_OMEGA,
                      "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("radial_grid", rowsB, outdir)

    # --- axis C: frequency sweep, per method, fixed grid length ---
    r = _r_grid(REP_A, FREQ_NS)
    rowsC = []
    for method in methods:
        for omega in omegas:
            t = bench(lambda method=method, omega=omega:
                      RadialTeukolsky(REP_S, REP_L, REP_M, REP_A, omega, r).solve(method))
            rowsC.append({"method": method, "omega": omega, "s": REP_S, "l": REP_L,
                          "m": REP_M, "a": REP_A, "nsamples": r.shape[0],
                          "min": t["min"], "median": t["median"], "p90": t["p90"]})
    return save("radial_freq", rowsC, outdir)
