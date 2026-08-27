"""Benchmark: full Teukolsky amplitude solve (TeukolskyMode.solve).

Two axes:

  * teuk_aggregate.csv -- TeukolskyMode.solve() time across orbit class, mode number, and
                          spin, at a few resolutions. This is the end-to-end cost users
                          actually pay per (s, l, m, k, n) mode. Includes a few high-l/n
                          spot checks (SPOTCHECK_MODES) beyond the curated MODES list, to
                          bound the cost of the large-(l,n) tail that shows up in larger
                          sweeps (e.g. benchmarks/regression_sweep.py, l up to 30/n up to
                          50) but isn't otherwise exercised here.
  * teuk_stages.csv    -- a per-stage breakdown at one representative mode per orbit
                          class. TeukolskyMode.solve() is a single monolithic C++ call
                          (the Python/cython API accepts precomputed teuk/swsh objects for
                          reuse, but the current binding does not actually consume them --
                          checked in cython/teukolsky_wrap.pyx), so there is no in-process
                          hook to time its internal stages directly. Instead this solves
                          the mode once to learn the *actual* radial/polar grids and
                          frequency it used, then independently times a standalone
                          SpinWeightedHarmonic and RadialTeukolsky built on those same
                          grids/parameters -- real measurements, not a cheap proxy (see the
                          note in the Performance docs about why an earlier eigenvalue-only
                          proxy was misleading). The remainder (full minus swsh minus
                          radial) is reported in the docs as a derived estimate of the
                          source-integration + orchestration cost, not a direct measurement.

  * teuk_flux.csv      -- FluxMode cost measured against the TeukolskyMode.solve it
                          consumes, at the same mode. FluxMode(geo, teuk) is pure
                          post-processing of an already-solved mode, so the flux page in
                          the Performance docs reports the mode-solve cost rather than a
                          flux-specific timing; this dataset is what justifies that.

Emits benchmarks/data/teuk_aggregate.csv, benchmarks/data/teuk_stages.csv, and
benchmarks/data/teuk_flux.csv.
"""
import warnings

from .timing import bench, save
from .grids import FIELD_SPINS

warnings.filterwarnings("ignore")  # non-convergence warnings are expected at low nsamples

# (orbit_class, a, p, e, x)
ORBITS = [
    ("circular-equatorial", 0.9, 10.0, 0.0, 1.0),
    ("eccentric-equatorial", 0.9, 10.0, 0.5, 1.0),
    ("spherical", 0.9, 10.0, 0.0, 0.5),
    ("generic", 0.9, 7.0, 0.6, 0.1),
]

# curated (s, l, m, k, n); invalid/redundant combos per orbit class are filtered below
MODES = [
    (-2, 2, 2, 0, 0),
    (-2, 5, 3, -2, 3),
    (-2, 8, 4, 2, 10),
    (0, 2, 2, 0, 0),
    (0, 5, 3, -2, 3),
    (2, 2, 2, 0, 0),
    (2, 5, 3, -2, 3),
    (2, 8, 4, 2, 10),
]

# High-l/high-n spot checks, beyond anything above -- added to bound the cost of the
# large-(l,n) tail of the 137k-case regression grid (benchmarks/regression_sweep.py sweeps
# l up to 30 and n up to 50), which the modes above don't reach. Only the "generic" orbit
# keeps these (nonzero k/n are filtered out for the equatorial/circular classes by
# _mode_valid), which is also the class that dominates that grid's case count.
SPOTCHECK_MODES = [
    (-2, 13, 7, -5, 20),
    (-2, 20, 10, 5, 30),
    (-2, 30, 15, -8, 50),
]

# Flux modes: s=-2 only (FluxMode is built from the s=-2 Teukolsky amplitudes), spanning
# the low-l/high-l ends of MODES so the flux:solve ratio is bracketed rather than sampled
# at one point.
FLUX_MODES = [
    (-2, 2, 2, 0, 0),
    (-2, 5, 3, -2, 3),
    (-2, 8, 4, 2, 10),
]

TEUK_RESOLUTIONS = (2**6, 2**9, 2**12)

STAGE_NS = 512
# one representative mode per orbit class (s=-2, l=5, m=3), matching each class's k/n rule
STAGE_MODES = {
    "circular-equatorial": (-2, 5, 3, 0, 0),
    "eccentric-equatorial": (-2, 5, 3, 0, 3),
    "spherical": (-2, 5, 3, -2, 0),
    "generic": (-2, 5, 3, -2, 3),
}


def _mode_valid(s, l, m, k, n, e, x):
    if l < max(abs(s), abs(m)):
        return False
    if abs(x) == 1.0 and k != 0:   # equatorial: no polar structure
        return False
    if e == 0.0 and n != 0:        # circular: no radial harmonics
        return False
    return True


def run(outdir, quick=False):
    from pybhpt.geo import KerrGeodesic
    from pybhpt.teuk import TeukolskyMode
    from pybhpt.swsh import SpinWeightedHarmonic
    from pybhpt.radial import RadialTeukolsky
    from pybhpt.flux import FluxMode

    orbits = ORBITS[:2] if quick else ORBITS
    modes = MODES[:3] if quick else MODES + SPOTCHECK_MODES
    resolutions = (2**6, 2**10) if quick else TEUK_RESOLUTIONS

    # --- axis A: aggregate solve cost across the parameter space ---
    rowsA = []
    for cls, a, p, e, x in orbits:
        for ns in resolutions:
            orbit = KerrGeodesic(a, p, e, x, ns)
            for (s, l, m, k, n) in modes:
                if not _mode_valid(s, l, m, k, n, e, x):
                    continue
                t = bench(lambda s=s, l=l, m=m, k=k, n=n, orbit=orbit, ns=ns:
                          TeukolskyMode(s, l, m, k, n, orbit).solve(orbit, nsamples=ns))
                rowsA.append({"orbit_class": cls, "a": a, "p": p, "e": e, "x": x,
                              "s": s, "l": l, "m": m, "k": k, "n": n, "nsamples": ns,
                              "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("teuk_aggregate", rowsA, outdir)

    # --- axis B: per-stage breakdown, one mode per orbit class ---
    rowsB = []
    for cls, a, p, e, x in orbits:
        s, l, m, k, n = STAGE_MODES[cls]
        orbit = KerrGeodesic(a, p, e, x, STAGE_NS)

        probe = TeukolskyMode(s, l, m, k, n, orbit)
        probe.solve(orbit, nsamples=STAGE_NS)
        r = probe.radialpoints
        th = probe.polarpoints
        omega = probe.frequency
        gamma = a * omega

        t_full = bench(lambda s=s, l=l, m=m, k=k, n=n, orbit=orbit:
                       TeukolskyMode(s, l, m, k, n, orbit).solve(orbit, nsamples=STAGE_NS))
        t_swsh = bench(lambda s=s, l=l, m=m, gamma=gamma, th=th:
                       SpinWeightedHarmonic(s, l, m, gamma, th))
        t_radial = bench(lambda s=s, l=l, m=m, a=a, omega=omega, r=r:
                         RadialTeukolsky(s, l, m, a, omega, r).solve("AUTO"))

        for stage, t in (("full", t_full), ("swsh", t_swsh), ("radial", t_radial)):
            rowsB.append({"orbit_class": cls, "s": s, "l": l, "m": m, "k": k, "n": n,
                          "nsamples_r": len(r), "nsamples_th": len(th), "stage": stage,
                          "min": t["min"], "median": t["median"], "p90": t["p90"]})
    save("teuk_stages", rowsB, outdir)

    # --- axis C: flux post-processing vs the mode solve it consumes ---
    # FluxMode reads an already-solved TeukolskyMode, so its cost is reported next to that
    # solve; the flux docs page cites the ratio instead of carrying its own timing table.
    rowsC = []
    for cls, a, p, e, x in orbits:
        for ns in resolutions:
            orbit = KerrGeodesic(a, p, e, x, ns)
            for (s, l, m, k, n) in FLUX_MODES:
                if not _mode_valid(s, l, m, k, n, e, x):
                    continue
                solved = TeukolskyMode(s, l, m, k, n, orbit)
                solved.solve(orbit, nsamples=ns)
                t_solve = bench(lambda s=s, l=l, m=m, k=k, n=n, orbit=orbit, ns=ns:
                                TeukolskyMode(s, l, m, k, n, orbit).solve(orbit, nsamples=ns))
                t_flux = bench(lambda orbit=orbit, solved=solved: FluxMode(orbit, solved))
                for stage, tt in (("solve", t_solve), ("flux", t_flux)):
                    rowsC.append({"orbit_class": cls, "a": a, "p": p, "e": e, "x": x,
                                  "s": s, "l": l, "m": m, "k": k, "n": n, "nsamples": ns,
                                  "stage": stage, "min": tt["min"], "median": tt["median"],
                                  "p90": tt["p90"]})
    return save("teuk_flux", rowsC, outdir)
