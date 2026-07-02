"""pybhpt performance benchmarks.

Run locally on a quiet reference machine to (re)generate the static data under
benchmarks/data/, which the documentation renders. Timing is never done at doc-build
or in CI.

    python -m benchmarks.run --all            # full grids -> benchmarks/data/*.csv
    python -m benchmarks.run geodesic --quick  # fast subset (smoke)
"""
