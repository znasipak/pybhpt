"""Benchmark: KerrGeodesic construction+solve across (a, p, e, x) and resolution.

Cost is driven mainly by nsamples (grid) and orbit complexity (eccentric / inclined).
The geodesic is solved once per orbit and reused by every mode, so this is an amortised
per-orbit cost. Emits benchmarks/data/geodesic.csv (one row per valid orbit x resolution).
"""
from .timing import bench, save
from .grids import orbit_params, orbit_class, make_orbit, RESOLUTIONS


def run(outdir, quick=False):
    from pybhpt.geo import KerrGeodesic
    resolutions = (2**6, 2**9, 2**12) if quick else RESOLUTIONS
    rows = []
    for (a, p, e, x) in orbit_params():
        for ns in resolutions:
            if make_orbit(KerrGeodesic, a, p, e, x, ns) is None:
                continue   # below separatrix / invalid
            t = bench(lambda a=a, p=p, e=e, x=x, ns=ns: KerrGeodesic(a, p, e, x, ns))
            rows.append({"a": a, "p": p, "e": e, "x": x,
                         "orbit_class": orbit_class(e, x), "nsamples": ns,
                         "min": t["min"], "median": t["median"], "p90": t["p90"]})
    return save("geodesic", rows, outdir)
