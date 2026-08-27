"""Benchmark: homogeneous radial Teukolsky solutions (RadialTeukolsky.solve).

Five axes, each isolating one cost driver:

  * radial_method.csv -- solve time by SolutionMethod (AUTO/HBL/GSN) at a fixed
                          representative mode and radial grid. HBL/GSN are numerical ODE
                          integrators; AUTO picks whichever is expected to be
                          cheapest/most stable. The analytic MST series is deliberately
                          NOT benchmarked: it was removed from AUTO's dispatch in 6d214d5
                          (unreliable at high l), so AUTO can never select it, and it is
                          not a method the docs steer users toward.
  * radial_grid.csv    -- solve time vs radial grid length (nsamples), method AUTO, at a
                          fixed mode. Exposes the per-point cost of evaluating solutions
                          on the output grid on top of any fixed setup cost.
  * radial_freq.csv    -- solve time vs mode frequency omega, for each method, at a fixed
                          grid length. Higher |omega| means more oscillatory solutions and
                          the ODE integrators take more steps, so this exposes how each
                          method scales with frequency.

  * radial_mode.csv    -- solve time over mode number and black hole spin: l x m x a, at
                          fixed omega and grid length, method AUTO. Isolates how cost
                          scales with the mode itself (the angular/radial structure the
                          solver must resolve) rather than with the output grid or the
                          frequency. Rows that fail to solve are recorded with status
                          "error" and no timing, rather than dropped.

  * radial_random.csv  -- a random sample of the (l, m, a, omega) parameter space, behind
                          the Performance page's scatter plots. The other axes hold two of
                          the three drivers fixed, so none of them shows how l, a and omega
                          interact; this one samples all of them at once. Spin is drawn
                          log-uniformly in (1 - a), since the cost structure bunches
                          against a = 1 and a linear spin axis would hide it; omega is
                          drawn uniformly and plotted on a linear axis. Seeded, so the
                          same points come back on a re-run.

Emits benchmarks/data/radial_method.csv, radial_grid.csv, radial_freq.csv,
radial_mode.csv, radial_random.csv.
"""
import numpy as np

from .timing import bench, save
from .grids import RESOLUTIONS

METHODS = ("AUTO", "HBL", "GSN")

# representative mode for the method-comparison and grid-scaling axes
REP_S, REP_L, REP_M, REP_A, REP_OMEGA = -2, 5, 2, 0.9, 0.3

OMEGAS = (0.02, 0.05, 0.1, 0.3, 0.6, 1.2, 2.5)
FREQ_NS = 256  # fixed grid length for the frequency sweep

# mode/spin sweep (axis D): omega and grid length are held fixed so the only thing varying
# is the mode and the background spin.
MODE_SPINS = (0.0, 0.3, 0.7, 0.9, 0.99)
MODE_LS = (2, 5, 10, 15, 30, 50)
MODE_NS = 256
MODE_M_COUNT = 5   # up to 5 equally spaced m per l, spanning -l..l


# random sample (axis E) of the joint (l, m, a, omega) space
RANDOM_N = 600
RANDOM_SEED = 20260827
RANDOM_L_RANGE = (2, 50)          # uniform integer
RANDOM_ONE_MINUS_A = (0.01, 1.0)  # log-uniform in (1 - a), i.e. a in [0, 0.99]
RANDOM_OMEGA = (0.02, 5.0)        # uniform
RANDOM_NS = 256


def _m_values(l, count=MODE_M_COUNT):
    """Up to `count` equally spaced m in [-l, l] (fewer when 2l+1 < count)."""
    if 2 * l + 1 <= count:
        return tuple(range(-l, l + 1))
    return tuple(sorted({int(round(v)) for v in np.linspace(-l, l, count)}))


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
    save("radial_freq", rowsC, outdir)

    # --- axis D: mode number x spin, fixed omega and grid length, method AUTO ---
    spins = MODE_SPINS[::2] if quick else MODE_SPINS
    ls = (2, 10, 30) if quick else MODE_LS
    rowsD = []
    for a in spins:
        r = _r_grid(a, MODE_NS)
        for l in ls:
            for m in _m_values(l):
                row = {"a": a, "s": REP_S, "l": l, "m": m, "omega": REP_OMEGA,
                       "nsamples": MODE_NS, "status": "ok"}
                try:
                    t = bench(lambda a=a, l=l, m=m, r=r:
                              RadialTeukolsky(REP_S, l, m, a, REP_OMEGA, r).solve("AUTO"))
                    row.update(min=t["min"], median=t["median"], p90=t["p90"])
                except Exception as exc:
                    row.update(status=f"error: {type(exc).__name__}",
                               min="", median="", p90="")
                rowsD.append(row)
    save("radial_mode", rowsD, outdir)

    # --- axis E: random sample of (l, m, a, omega), method AUTO ---
    # One random m per point (uniform in [-l, l]) rather than a median over several: the
    # sample is of the mode space itself, so each row is a real mode, not an aggregate.
    rng = np.random.default_rng(RANDOM_SEED)
    n_samples = 60 if quick else RANDOM_N
    rowsE = []
    for _ in range(n_samples):
        l = int(rng.integers(RANDOM_L_RANGE[0], RANDOM_L_RANGE[1] + 1))
        m = int(rng.integers(-l, l + 1))
        one_minus_a = float(np.exp(rng.uniform(np.log(RANDOM_ONE_MINUS_A[0]),
                                               np.log(RANDOM_ONE_MINUS_A[1]))))
        a = 1.0 - one_minus_a
        omega = float(rng.uniform(*RANDOM_OMEGA))
        r = _r_grid(a, RANDOM_NS)
        row = {"l": l, "m": m, "a": a, "one_minus_a": one_minus_a, "omega": omega,
               "s": REP_S, "nsamples": RANDOM_NS, "status": "ok"}
        try:
            tt = bench(lambda l=l, m=m, a=a, omega=omega, r=r:
                       RadialTeukolsky(REP_S, l, m, a, omega, r).solve("AUTO"))
            row.update(min=tt["min"], median=tt["median"], p90=tt["p90"])
        except Exception as exc:
            row.update(status=f"error: {type(exc).__name__}", min="", median="", p90="")
        rowsE.append(row)
    return save("radial_random", rowsE, outdir)
