"""Generate the static PNG figures embedded in the Performance docs pages.

Reads the committed CSVs under benchmarks/data/ and writes PNGs straight into the docs
tree (default docs/_static/figures/), which is where the Performance pages reference them
from -- one committed copy, no duplicate under benchmarks/. Like the timing harness, this
runs locally and the resulting images are committed; nothing is generated at
documentation-build time or in CI.

    python -m benchmarks.plots                 # -> docs/_static/figures/*.png
    python -m benchmarks.plots --out <dir>

Three styles are produced per dataset (scatter vs resolution, histogram of the time
distribution, and a bar/strip category comparison) so the most legible one can be chosen
per component.
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MS = 1e3  # seconds -> milliseconds

# consistent colors for the orbit classes / spins
_CLASS_COLOR = {
    "circular-equatorial": "#4c72b0",
    "eccentric-equatorial": "#dd8452",
    "spherical": "#55a868",
    "generic": "#c44e52",
}


def _load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (ValueError, TypeError):
                pass
    return rows


def _save(fig, outdir, name):
    fig.tight_layout()
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ------------------------------- KerrGeodesic -------------------------------
def geo_scatter(rows, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for cls, color in _CLASS_COLOR.items():
        sub = [r for r in rows if r["orbit_class"] == cls]
        if not sub:
            continue
        x = [r["nsamples"] for r in sub]
        y = [r["median"] * MS for r in sub]
        ax.scatter(x, y, s=22, alpha=0.7, color=color, label=cls)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("nsamples")
    ax.set_ylabel("median time [ms]")
    ax.set_title("KerrGeodesic: construction time vs resolution")
    ax.legend(fontsize=8, framealpha=0.9)
    ax.grid(True, which="both", alpha=0.25)
    return _save(fig, outdir, "geo_scatter.png")


def geo_hist(rows, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    allt = np.array([r["median"] * MS for r in rows])
    bins = np.logspace(np.log10(allt.min()), np.log10(allt.max()), 30)
    for cls, color in _CLASS_COLOR.items():
        t = [r["median"] * MS for r in rows if r["orbit_class"] == cls]
        if not t:
            continue
        ax.hist(t, bins=bins, alpha=0.6, color=color, label=cls)
    ax.set_xscale("log")
    ax.set_xlabel("median time [ms]")
    ax.set_ylabel("count (cases across a, p, e, x, nsamples)")
    ax.set_title("KerrGeodesic: distribution of construction times")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    return _save(fig, outdir, "geo_hist.png")


def geo_bar(rows, outdir):
    ns = max(r["nsamples"] for r in rows)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    classes = [c for c in _CLASS_COLOR if any(
        r["orbit_class"] == c and r["nsamples"] == ns for r in rows)]
    meds, p90s, colors = [], [], []
    for cls in classes:
        sub = [r for r in rows if r["orbit_class"] == cls and r["nsamples"] == ns]
        meds.append(np.median([r["median"] * MS for r in sub]))
        p90s.append(np.max([r["p90"] * MS for r in sub]))
        colors.append(_CLASS_COLOR[cls])
    xpos = np.arange(len(classes))
    err = np.array(p90s) - np.array(meds)
    ax.bar(xpos, meds, yerr=[np.zeros_like(err), err], color=colors,
           capsize=4, alpha=0.85)
    ax.set_yscale("log")
    ax.set_xticks(xpos)
    ax.set_xticklabels(classes, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("time [ms]")
    ax.set_title(f"KerrGeodesic: median & p90 whisker at nsamples={int(ns)}")
    ax.grid(True, axis="y", alpha=0.25)
    return _save(fig, outdir, "geo_bar.png")


# ------------------------------ SWSH: on-grid ------------------------------
def _mode_label(r):
    return f"(s={int(r['s'])}, l={int(r['l'])}, m={int(r['m'])}, γ={r['gamma']:g})"


def swsh_grid_scatter(rows, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    modes = sorted({(r["s"], r["l"], r["m"], r["gamma"]) for r in rows})
    cmap = plt.get_cmap("viridis", len(modes))
    for i, mode in enumerate(modes):
        sub = [r for r in rows if (r["s"], r["l"], r["m"], r["gamma"]) == mode]
        sub.sort(key=lambda r: r["nsamples"])
        ax.plot([r["nsamples"] for r in sub], [r["median"] * MS for r in sub],
                "o-", color=cmap(i), ms=4, label=_mode_label(sub[0]))
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("nsamples")
    ax.set_ylabel("median time [ms]")
    ax.set_title("SpinWeightedHarmonic: construction time vs grid length")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    return _save(fig, outdir, "swsh_grid_scatter.png")


# ------------------------------ SWSH: solve ------------------------------
def swsh_solve_hist(rows, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    gammas = sorted({r["gamma"] for r in rows})
    allt = np.array([r["median"] * MS for r in rows])
    bins = np.logspace(np.log10(allt.min()), np.log10(allt.max()), 24)
    cmap = plt.get_cmap("plasma", len(gammas))
    for i, g in enumerate(gammas):
        t = [r["median"] * MS for r in rows if r["gamma"] == g]
        ax.hist(t, bins=bins, alpha=0.6, color=cmap(i), label=f"γ={g:g}")
    ax.set_xscale("log")
    ax.set_xlabel("median solve time [ms]")
    ax.set_ylabel("count (cases across s, l)")
    ax.set_title("SpinWeightedHarmonic: spectral-solve time distribution")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    return _save(fig, outdir, "swsh_solve_hist.png")


def swsh_solve_bar(rows, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    gammas = sorted({r["gamma"] for r in rows})
    meds, p90s = [], []
    for g in gammas:
        sub = [r for r in rows if r["gamma"] == g]
        meds.append(np.median([r["median"] * MS for r in sub]))
        p90s.append(np.max([r["p90"] * MS for r in sub]))
    xpos = np.arange(len(gammas))
    cmap = plt.get_cmap("plasma", len(gammas))
    err = np.array(p90s) - np.array(meds)
    ax.bar(xpos, meds, yerr=[np.zeros_like(err), err],
           color=[cmap(i) for i in range(len(gammas))], capsize=4, alpha=0.85)
    ax.set_yscale("log")
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"γ={g:g}" for g in gammas])
    ax.set_ylabel("solve time [ms]")
    ax.set_title("SpinWeightedHarmonic: spectral-solve median & p90")
    ax.grid(True, axis="y", alpha=0.25)
    return _save(fig, outdir, "swsh_solve_bar.png")


# --------------------- RadialTeukolsky: l x a x omega cube ---------------------
def _random_scatter(rows, outdir, xkey, ykey, xlabel, ylabel, title, name,
                    xlog=True, ylog=False, invert_x=False):
    """Scatter of a random (l, m, a, omega) sample, colored by median solve time.

    Spin is plotted as 1 - a on a log axis: the interesting structure is bunched against
    a = 1, and a linear spin axis compresses everything from a = 0.9 to 0.99 -- most of
    the cost variation -- into the last tenth of the plot.
    """
    x = np.array([r[xkey] for r in rows])
    y = np.array([r[ykey] for r in rows])
    c = np.array([r["median"] * MS for r in rows])

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    sc = ax.scatter(x, y, c=c, s=26, cmap="viridis",
                    norm=matplotlib.colors.LogNorm(vmin=c.min(), vmax=c.max()),
                    linewidths=0)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    if invert_x:
        ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_axisbelow(True)
    ax.grid(True, which="both", alpha=0.2)
    cb = fig.colorbar(sc, ax=ax, label="median solve time [ms]")
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    return _save(fig, outdir, name)


def radial_random_l_vs_a(rows, outdir):
    return _random_scatter(
        rows, outdir, "one_minus_a", "l", "1 − a  (near-extremal to the right)",
        "mode number l", "RadialTeukolsky: l vs spin",
        "radial_random_l_vs_a.png", invert_x=True)


def radial_random_l_vs_omega(rows, outdir):
    return _random_scatter(
        rows, outdir, "omega", "l", "frequency omega", "mode number l",
        "RadialTeukolsky: l vs omega", "radial_random_l_vs_omega.png", xlog=False)


def radial_random_omega_vs_a(rows, outdir):
    return _random_scatter(
        rows, outdir, "one_minus_a", "omega", "1 − a  (near-extremal to the right)",
        "frequency omega", "RadialTeukolsky: omega vs spin",
        "radial_random_omega_vs_a.png", ylog=False, invert_x=True)


# ------------------- SWSH: random (l, m, gamma) sample -------------------
# Same treatment as the radial random sample, with one continuous axis instead of two:
# the radial solver takes spin and frequency separately, the harmonic sees only their
# product gamma = a*omega.
def swsh_random_l_vs_gamma(rows, outdir):
    return _random_scatter(
        rows, outdir, "gamma", "l", "spheroidicity gamma", "mode number l",
        "SpinWeightedHarmonic: l vs spheroidicity",
        "swsh_random_l_vs_gamma.png", xlog=False)


def swsh_random_m_vs_gamma(rows, outdir):
    return _random_scatter(
        rows, outdir, "gamma", "m", "spheroidicity gamma", "azimuthal number m",
        "SpinWeightedHarmonic: m vs spheroidicity",
        "swsh_random_m_vs_gamma.png", xlog=False)


def swsh_random_l_vs_m(rows, outdir):
    return _random_scatter(
        rows, outdir, "m", "l", "azimuthal number m", "mode number l",
        "SpinWeightedHarmonic: l vs m", "swsh_random_l_vs_m.png", xlog=False)


def run(datadir, outdir):
    os.makedirs(outdir, exist_ok=True)
    geo = _load(os.path.join(datadir, "geodesic.csv"))
    grid = _load(os.path.join(datadir, "swsh_grid.csv"))
    solve = _load(os.path.join(datadir, "swsh_solve.csv"))
    rand = [r for r in _load(os.path.join(datadir, "radial_random.csv"))
            if r.get("status") == "ok"]
    swsh_rand = [r for r in _load(os.path.join(datadir, "swsh_random.csv"))
                 if r.get("status") == "ok"]
    paths = [
        geo_scatter(geo, outdir), geo_hist(geo, outdir), geo_bar(geo, outdir),
        swsh_grid_scatter(grid, outdir),
        swsh_solve_hist(solve, outdir), swsh_solve_bar(solve, outdir),
        radial_random_l_vs_a(rand, outdir),
        radial_random_l_vs_omega(rand, outdir),
        radial_random_omega_vs_a(rand, outdir),
        swsh_random_l_vs_gamma(swsh_rand, outdir),
        swsh_random_m_vs_gamma(swsh_rand, outdir),
        swsh_random_l_vs_m(swsh_rand, outdir),
    ]
    for p in paths:
        print("wrote", p)
    return paths


if __name__ == "__main__":
    here = os.path.dirname(__file__)
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(here, "data"))
    ap.add_argument("--out", default=os.path.join(here, os.pardir, "docs",
                                                  "_static", "figures"))
    args = ap.parse_args()
    run(args.data, args.out)
