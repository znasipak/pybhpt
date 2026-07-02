"""Timing utilities for the pybhpt performance benchmarks.

Reports the min, median, and p90 of several timed repetitions (after warmup), and
records the platform/build so committed results are interpretable. Timing is meant to
be run locally on a quiet reference machine and the resulting data committed under
benchmarks/data/ -- the documentation renders those static files; nothing is timed at
doc-build or in CI.
"""
import json
import platform
import subprocess
import time

import numpy as np


def bench(fn, reps=7, warmup=2):
    """Time `fn` (no args). Returns {min, median, p90} in seconds.

    warmup runs are discarded (first-call effects: lazy tables, page faults, caches).
    median is the headline number; p90 exposes the tail; min is the best case.
    """
    for _ in range(warmup):
        fn()
    ts = np.empty(reps)
    for i in range(reps):
        t0 = time.perf_counter()
        fn()
        ts[i] = time.perf_counter() - t0
    return {"min": float(ts.min()),
            "median": float(np.median(ts)),
            "p90": float(np.percentile(ts, 90))}


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True,
            stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def _blas_backend():
    """Best-effort name of the BLAS/LAPACK numpy (and, by proxy, the build) links."""
    try:
        cfg = np.show_config(mode="dicts")  # numpy >= 1.25
        for _, v in cfg.get("Build Dependencies", {}).items():
            if isinstance(v, dict) and "name" in v:
                return v["name"]
    except Exception:
        pass
    return "unknown"


def platform_info():
    """Machine/build stamp recorded alongside every benchmark data file."""
    import pybhpt  # noqa: F401  (record that the import works / its location)
    return {
        "cpu": platform.processor() or platform.machine(),
        "machine": platform.machine(),
        "system": f"{platform.system()} {platform.release()}",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "blas": _blas_backend(),
        "pybhpt_commit": _git_sha(),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
    }


def save(name, rows, outdir, meta=None):
    """Write a benchmark table to <outdir>/<name>.csv and merge platform metadata.

    rows: list of dicts (parameter columns + min/median/p90 seconds).
    """
    import csv
    import os
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}.csv")
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    mpath = os.path.join(outdir, "metadata.json")
    allmeta = {}
    if os.path.exists(mpath):
        allmeta = json.load(open(mpath))
    allmeta[name] = {"platform": platform_info(), **(meta or {})}
    json.dump(allmeta, open(mpath, "w"), indent=2)
    return path
