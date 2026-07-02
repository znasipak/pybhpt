"""Large-scale regression sweep: dump amplitudes + precisions for many
(s, l, m, k, n) x (a, p, e, x) x resolution, then compare two builds
(this branch vs main).

    # build under test loaded explicitly (editable install hijacks cybhpt_full)
    python dev/regression_sweep.py --module <branch>/cybhpt_full*.so --dump branch.npz
    python dev/regression_sweep.py --module <main>/cybhpt_full*.so   --dump main.npz
    python dev/regression_sweep.py --compare branch.npz main.npz

Uses only the API common to both builds: KerrGeodesic(a,p,e,x,ns),
TeukolskyMode(s,l,m,k,n,orbit).solve(orbit), amplitude('In'/'Up'), precision('In'/'Up').
"""
import argparse, importlib.util, sys, warnings, itertools
import numpy as np
warnings.filterwarnings("ignore")

# ------------------------------- parameter grid -------------------------------
# Locked grid (~137k cases). Retrograde orbits (x<0), high eccentricity (e=0.9),
# L up to 30, m+k swept in {-10,-2,0,2,10}, n up to 50 for e in {0.5, 0.9}.
SPINS = (-2, 0, 2)
LS = (2, 3, 5, 13, 20, 30)
MK = (-10, -2, 0, 2, 10)                 # sweep m+k; polar harmonic k = (m+k) - m
A_SPINS = (0.0, 0.5, 0.9, 0.99)
PE = ((10.0, 0.0), (10.0, 0.5), (8.0, 0.3), (12.0, 0.9))
RESOLUTIONS = (2**7, 2**12)              # >=128 so n=50 is representable; coarse + converged


def x_values(p, e):
    # full inclination coverage only for the high-eccentricity orbit; else prograde-
    # and retrograde-inclined only (equatorial x=+-1 kept just for e=0.9).
    return (1.0, 0.5, 0.1, -0.3, -1.0) if e == 0.9 else (0.5, -0.3)


def m_values(l):
    # coarser m sampling for large l to bound the combinatorics
    return sorted(set([-l, 0, l // 2, l])) if l >= 13 else sorted(set([-l, -1, 0, 1, l // 2, l]))


def n_values(e):
    if e == 0.0:
        return (0,)
    if e <= 0.3:
        return (-3, 0, 3, 10)
    return (-3, 0, 3, 10, 20, 50)        # high radial harmonics for e in {0.5, 0.9}


def iter_cases(limit=None):
    count = 0
    for a in A_SPINS:
        for (p, e) in PE:
            for x in x_values(p, e):
                equatorial = (abs(x) == 1.0)
                for ns in RESOLUTIONS:
                    for s in SPINS:
                        for l in LS:
                            if l < max(abs(s), 1):
                                continue
                            for m in m_values(l):
                                if abs(m) > l:
                                    continue
                                for t in MK:
                                    # equatorial: only k=0 contributes (keep only t == m)
                                    if equatorial and t != m:
                                        continue
                                    k = 0 if equatorial else (t - m)
                                    for n in n_values(e):
                                        if e == 0.0 and n != 0:   # circular: only n=0
                                            continue
                                        yield (s, l, m, k, n, a, p, e, x, ns)
                                        count += 1
                                        if limit and count >= limit:
                                            return


def run_dump(module, out, limit, pybhpt_path=None):
    spec = importlib.util.spec_from_file_location("cybhpt_full", module)
    cyb = importlib.util.module_from_spec(spec); sys.modules["cybhpt_full"] = cyb
    spec.loader.exec_module(cyb)
    if pybhpt_path:   # use the matching repo's pure-python pybhpt wrapper (main vs branch)
        # drop the editable meta-path finder (it redirects pybhpt to the installed branch,
        # overriding sys.path); cybhpt_full is already registered above.
        sys.meta_path = [f for f in sys.meta_path
                         if type(f).__name__ != "ScikitBuildRedirectingFinder"]
        sys.path.insert(0, pybhpt_path)
    from pybhpt.geo import KerrGeodesic
    from pybhpt.teuk import TeukolskyMode
    print("cybhpt:", cyb.__file__, "| pybhpt:", sys.modules["pybhpt"].__file__)

    keys, ampIn, ampUp, precIn, precUp = [], [], [], [], []
    orbit_cache = {}
    done = 0
    for (s, l, m, k, n, a, p, e, x, ns) in iter_cases(limit):
        key = f"{s}:{l}:{m}:{k}:{n}:{a}:{p}:{e}:{x}:{ns}"
        try:
            ok = (a, p, e, x, ns)
            orbit = orbit_cache.get(ok)
            if orbit is None:
                orbit = KerrGeodesic(a, p, e, x, ns)
                orbit_cache[ok] = orbit
            md = TeukolskyMode(s, l, m, k, n, orbit)
            md.solve(orbit)
            ai, au = md.amplitude("In"), md.amplitude("Up")
            pi, pu = md.precision("In"), md.precision("Up")
        except Exception:
            ai = au = complex(np.nan); pi = pu = np.nan
        keys.append(key); ampIn.append(ai); ampUp.append(au); precIn.append(pi); precUp.append(pu)
        done += 1
        if done % 2000 == 0:
            print(f"  {done} modes...", flush=True)
    np.savez(out, keys=np.array(keys), ampIn=np.array(ampIn, dtype=complex),
             ampUp=np.array(ampUp, dtype=complex),
             precIn=np.array(precIn), precUp=np.array(precUp))
    print(f"wrote {out}: {done} modes  ({cyb.__file__})")


def _reldiff(x, y):
    x = np.asarray(x); y = np.asarray(y)
    d = np.abs(x - y)
    scale = np.maximum(np.abs(x), np.abs(y))
    out = np.where(scale > 0, d / np.where(scale == 0, 1, scale), 0.0)
    return out


def run_compare(fa, fb):
    A = np.load(fa, allow_pickle=True); B = np.load(fb, allow_pickle=True)
    ka = {k: i for i, k in enumerate(A["keys"])}
    kb = {k: i for i, k in enumerate(B["keys"])}
    common = [k for k in A["keys"] if k in kb]
    print(f"{len(A['keys'])} vs {len(B['keys'])} rows; {len(common)} common keys\n")
    ia = np.array([ka[k] for k in common]); ib = np.array([kb[k] for k in common])

    for side, amp, prec in (("In", "ampIn", "precIn"), ("Up", "ampUp", "precUp")):
        aA, aB = A[amp][ia], B[amp][ib]
        pA, pB = A[prec][ia], B[prec][ib]
        finite = np.isfinite(aA) & np.isfinite(aB) & np.isfinite(np.abs(aA)) & np.isfinite(np.abs(aB))
        # amplitude relative difference where at least one side is non-negligible
        rel = _reldiff(aA, aB)
        big = finite & (np.maximum(np.abs(aA), np.abs(aB)) > 1e-30)
        r = rel[big]
        print(f"=== {side} amplitude (branch=A vs main=B), {big.sum()} nonzero finite ===")
        if r.size:
            for q in (50, 90, 99, 100):
                print(f"   {q:3d}th pct rel-diff: {np.percentile(r, q):.2e}")
            # worst offenders
            order = np.argsort(-r)[:10]
            ck = np.array(common)[big]
            print("   worst:")
            for j in order[:8]:
                print(f"     {ck[j]:38s} rel={r[j]:.2e}  A={aA[big][j]:.3e} B={aB[big][j]:.3e}"
                      f"  precA={pA[big][j]:.1e} precB={pB[big][j]:.1e}")
        # precision comparison
        pf = finite & np.isfinite(pA) & np.isfinite(pB) & (pB > 0)
        ratio = (pA[pf] / pB[pf])
        if ratio.size:
            print(f"   precision ratio A/B: median={np.median(ratio):.2f} "
                  f"10th={np.percentile(ratio,10):.2f} 90th={np.percentile(ratio,90):.2f}")
        else:
            print("   precision ratio A/B: (no comparable rows)")
        print()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--module")
    ap.add_argument("--dump")
    ap.add_argument("--compare", nargs=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pybhpt", default=None, help="repo dir whose pybhpt/ python to use")
    args = ap.parse_args()
    if args.compare:
        run_compare(*args.compare)
    else:
        assert args.module and args.dump
        run_dump(args.module, args.dump, args.limit, args.pybhpt)
