#!/usr/bin/env bash
# Regression sweep: this branch vs main, over the locked ~137k-case grid defined in
# benchmarks/regression_sweep.py. Runs LOCALLY on a quiet machine (main is slow; expect
# several hours to overnight). Nothing here runs in CI.
#
# Usage:
#   bash benchmarks/regression/run_regression.sh
# Optional env:
#   MAIN_REF   git ref for the baseline build         (default: origin/main)
#   ENV_BIN    dir of the python/cmake toolchain       (default: current PATH)
#   WORKTREE   path for the temporary main checkout    (default: /tmp/pybhpt-main)
#   OUTDIR     where dumps + summary land              (default: benchmarks/regression/data)
#   QUICK_LIMIT  cap cases for a smoke run (unset = full grid)
set -euo pipefail

BRANCH_REPO="$(cd "$(dirname "$0")/../.." && pwd)"
MAIN_REF="${MAIN_REF:-origin/main}"
WORKTREE="${WORKTREE:-/tmp/pybhpt-main}"
OUTDIR="${OUTDIR:-$BRANCH_REPO/benchmarks/regression/data}"
BRANCH_SO="$(ls "$BRANCH_REPO"/build_opt/cybhpt_full*.so 2>/dev/null | head -1)"
LIMIT_ARG=""; [ -n "${QUICK_LIMIT:-}" ] && LIMIT_ARG="--limit $QUICK_LIMIT"
mkdir -p "$OUTDIR"

echo "== 1. build branch (expects an existing build_opt) =="
[ -n "$BRANCH_SO" ] || { echo "no build_opt .so; run: cmake --build build_opt -j8"; exit 1; }

echo "== 2. build main worktree at $MAIN_REF =="
git -C "$BRANCH_REPO" fetch origin main --quiet || true
rm -rf "$WORKTREE"; git -C "$BRANCH_REPO" worktree prune
git -C "$BRANCH_REPO" worktree add "$WORKTREE" "$MAIN_REF"
# boost is a submodule with pre-built headers in the branch checkout: symlink it in
rm -rf "$WORKTREE/extern/boost"
ln -s "$BRANCH_REPO/extern/boost" "$WORKTREE/extern/boost"
( cd "$WORKTREE" && cmake -S . -B build_main -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$(command -v python3)" >/dev/null && cmake --build build_main -j 8 )
MAIN_SO="$(ls "$WORKTREE"/build_main/cybhpt_full*.so | head -1)"

echo "== 3. dump both builds (each with its own matched pybhpt python) =="
SWEEP="$BRANCH_REPO/benchmarks/regression_sweep.py"
python "$SWEEP" --module "$BRANCH_SO" --pybhpt "$BRANCH_REPO" --dump "$OUTDIR/branch.npz" $LIMIT_ARG
python "$SWEEP" --module "$MAIN_SO"   --pybhpt "$WORKTREE"    --dump "$OUTDIR/main.npz"   $LIMIT_ARG

echo "== 4. compare -> $OUTDIR/compare.txt =="
python "$SWEEP" --compare "$OUTDIR/branch.npz" "$OUTDIR/main.npz" | tee "$OUTDIR/compare.txt"

echo "Done. Raw dumps (~large) in $OUTDIR; regenerate REPORT.md from them."
