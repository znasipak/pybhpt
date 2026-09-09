"""Validity-aware parameter grids for the benchmarks.

Orbit validity (p above the separatrix for the given a, e, x) is enforced by attempting
construction and skipping anything that raises -- pybhpt does not expose a separatrix
helper, and this is robust to the exact bound. Mode validity uses l >= max(|s|, |m|).
"""

# Orbit sampling: spins a, and (p, e, x) covering circular/eccentric x equatorial/inclined.
ORBIT_SPINS = (0.0, 0.5, 0.9, 0.99)
ORBIT_PE = ((10.0, 0.0), (10.0, 0.5), (8.0, 0.3))   # (p, e): circular, eccentric, moderate
ORBIT_X = (1.0, 0.5, 0.1)                            # equatorial, inclined, near-polar

RESOLUTIONS = (2**6, 2**7, 2**8, 2**9, 2**10, 2**11, 2**12)

FIELD_SPINS = (-2, 0, 2)


def orbit_params():
    """Yield (a, p, e, x). Caller must guard construction (validity by try/except)."""
    for a in ORBIT_SPINS:
        for (p, e) in ORBIT_PE:
            for x in ORBIT_X:
                yield (a, p, e, x)


def orbit_class(e, x):
    if e == 0.0 and abs(x) == 1.0:
        return "circular-equatorial"
    if e == 0.0:
        return "spherical"
    if abs(x) == 1.0:
        return "eccentric-equatorial"
    return "generic"


def make_orbit(KerrGeodesic, a, p, e, x, nsamples):
    """Construct an orbit, or None if the parameters are invalid (below separatrix)."""
    try:
        return KerrGeodesic(a, p, e, x, nsamples)
    except Exception:
        return None


def mode_params(s, l_values=(2, 3, 5, 8, 13), k_values=(-2, 0, 2), n_values=(-3, 0, 3, 10)):
    """Yield valid (s, l, m, k, n): l >= max(|s|, |m|), a spread of m."""
    for l in l_values:
        lmin = max(abs(s), 1)
        if l < lmin:
            continue
        for m in sorted(set([-l, -1, 0, 1, l // 2, l])):
            if abs(m) > l:
                continue
            for k in k_values:
                for n in n_values:
                    yield (s, l, m, k, n)
