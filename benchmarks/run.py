"""Benchmark dispatcher. Pins BLAS/OpenMP to a single thread for reproducible
single-core timings BEFORE importing numpy/pybhpt (multi-thread scaling is documented
separately). Writes benchmarks/data/<name>.csv + metadata.json.
"""
import argparse
import os

# must precede numpy / LAPACK load
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

DEFAULT_OUT = os.path.join(os.path.dirname(__file__), "data")

# name -> module providing run(outdir, quick)
BENCHMARKS = {
    "geodesic": "benchmarks.bench_geodesic",
    "swsh":    "benchmarks.bench_swsh",       # spectral solve + on-grid eval
    # "radial":  "benchmarks.bench_radial",   # solve vs grid length x method (TODO)
    # "teuk":    "benchmarks.bench_teuk",      # aggregate mode solve + stage breakdown (TODO)
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="?", default="all", choices=["all", *BENCHMARKS])
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--quick", action="store_true", help="fast subset for smoke tests")
    args = ap.parse_args()
    import importlib
    names = list(BENCHMARKS) if args.which == "all" else [args.which]
    for name in names:
        mod = importlib.import_module(BENCHMARKS[name])
        path = mod.run(args.out, quick=args.quick)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
