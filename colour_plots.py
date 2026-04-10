"""
plot_colour_curves.py
---------------------
Reads one or more SN photometry .txt files (as produced by e.g. the
Las Cumbres / AAVSO pipeline) and plots colour curves of the form
(filter_X - filter_Y) vs MJD, with propagated photometric uncertainties.

File format expected
--------------------
Line 1: <SN_name>  <redshift>  <RA>  <Dec>
Subsequent blocks:
  "filter <name>"  — starts a new filter block
  <MJD>  <mag>  <mag_err>  — one observation per line

Usage
-----
    python plot_colour_curves.py <file1.txt> [file2.txt ...]

Colour pairs to plot are defined in COLOUR_PAIRS below.
"""

import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path


# ── Configuration ─────────────────────────────────────────────────────────────

# Define which colour curves to plot.
# Each tuple is (filter_blue, filter_red), plotted as blue - red.
# The script will only plot a pair if both filters are present in the data.
COLOUR_PAIRS = [
    ("Bs", "Vs"),   # B - V
    ("Vs", "Rs"),   # V - R
    ("Rs", "Is"),   # R - I
    ("Bs", "Rs"),   # B - R
    ("Vs", "Is"),   # V - I
]

# Sigma-clipping threshold for outlier rejection (set to None to disable)
SIGMA_CLIP = None

# Maximum allowed photometric uncertainty; points above this are excluded
MAX_ERR = 0.7 # mag


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_file(filepath):
    """
    Parse a photometry text file.

    Returns
    -------
    meta : dict
        Keys: 'name', 'redshift', 'ra', 'dec'
    data : dict
        Keyed by filter name; each value is a dict with arrays
        'mjd', 'mag', 'mag_err'.
    """
    meta = {}
    data = {}
    current_filter = None

    with open(filepath, "r") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    if not lines:
        raise ValueError(f"Empty file: {filepath}")

    # First non-empty line is the header
    header = lines[0].split()
    if len(header) < 4:
        raise ValueError(
            f"Header line must have at least 4 fields "
            f"(name z RA Dec), got: {lines[0]!r}"
        )
    meta["name"] = header[0]
    meta["redshift"] = float(header[1])
    meta["ra"] = float(header[2])
    meta["dec"] = float(header[3])

    mjd_buf, mag_buf, err_buf = [], [], []

    def _flush():
        if current_filter is not None and mjd_buf:
            data[current_filter] = {
                "mjd":     np.array(mjd_buf,  dtype=float),
                "mag":     np.array(mag_buf,  dtype=float),
                "mag_err": np.array(err_buf,  dtype=float),
            }

    for line in lines[1:]:
        parts = line.split()
        if parts[0].lower() == "filter":
            _flush()
            current_filter = parts[1]
            mjd_buf, mag_buf, err_buf = [], [], []
            data.setdefault(current_filter, None)
        elif len(parts) == 3:
            try:
                mjd_buf.append(float(parts[0]))
                mag_buf.append(float(parts[1]))
                err_buf.append(float(parts[2]))
            except ValueError:
                pass  # skip malformed lines silently

    _flush()

    # Remove placeholder None entries (filters declared but empty)
    data = {k: v for k, v in data.items() if v is not None}

    return meta, data


# ── Colour construction ───────────────────────────────────────────────────────

def interpolate_to_common_mjd(d1, d2, tol=0.3):
    """
    Match observations in two filter datasets that fall within `tol` days
    of each other, averaging if multiple matches exist.

    Returns matched arrays: mjd_mid, colour, colour_err
    """
    mjd1, mag1, err1 = d1["mjd"], d1["mag"], d1["mag_err"]
    mjd2, mag2, err2 = d2["mjd"], d2["mag"], d2["mag_err"]

    mjd_out, col_out, col_err_out = [], [], []

    for i, t1 in enumerate(mjd1):
        dt = np.abs(mjd2 - t1)
        match = np.where(dt <= tol)[0]
        if match.size == 0:
            continue
        # Use the closest match if multiple
        j = match[np.argmin(dt[match])]
        t_mid = 0.5 * (t1 + mjd2[j])
        colour = mag1[i] - mag2[j]
        colour_err = np.sqrt(err1[i]**2 + err2[j]**2)
        mjd_out.append(t_mid)
        col_out.append(colour)
        col_err_out.append(colour_err)

    return (
        np.array(mjd_out),
        np.array(col_out),
        np.array(col_err_out),
    )


def sigma_clip(mjd, col, col_err, sigma=3.0):
    """Remove points deviating more than `sigma` from the median colour."""
    if len(col) < 4:
        return mjd, col, col_err
    med = np.median(col)
    std = np.std(col)
    mask = np.abs(col - med) < sigma * std
    return mjd[mask], col[mask], col_err[mask]


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_colour_curves(meta, data, colour_pairs, outfile=None):
    """
    Plot one panel per colour pair on a single figure.
    """
    # Work out which pairs are actually available
    available_pairs = [
        (f1, f2) for f1, f2 in colour_pairs
        if f1 in data and f2 in data
    ]

    if not available_pairs:
        print(
            f"[{meta['name']}] No matching colour pairs found in data. "
            f"Available filters: {list(data.keys())}"
        )
        return

    n_panels = len(available_pairs)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(8, 3.0 * n_panels),
        sharex=True,
        constrained_layout=True,
    )
    if n_panels == 1:
        axes = [axes]

    # Colour-map for multi-file use (here one object but extensible)
    colours = cm.tab10(np.linspace(0, 0.9, max(n_panels, 1)))

    fig.suptitle(
        f"{meta['name']}   "
        f"($z = {meta['redshift']:.4f}$,  "
        f"RA = {meta['ra']:.4f}°,  "
        f"Dec = {meta['dec']:.4f}°)",
        fontsize=12,
    )

    for ax, (f1, f2), colour in zip(axes, available_pairs, colours):
        d1 = data[f1]
        d2 = data[f2]

        # Apply per-filter quality cut on error
        def quality_cut(d):
            mask = d["mag_err"] <= MAX_ERR
            return {k: v[mask] for k, v in d.items()}

        d1q = quality_cut(d1)
        d2q = quality_cut(d2)

        mjd, col, col_err = interpolate_to_common_mjd(d1q, d2q)

        if mjd.size == 0:
            ax.set_visible(False)
            continue

        if SIGMA_CLIP is not None:
            mjd, col, col_err = sigma_clip(mjd, col, col_err, sigma=SIGMA_CLIP)

        # Plot
        ax.errorbar(
            mjd, col, yerr=col_err,
            fmt="o",
            color=colour,
            markersize=5,
            capsize=3,
            elinewidth=1.0,
            linewidth=0,
            label=f"{f1}−{f2}",
            zorder=3,
        )

        # Light connecting line to show temporal evolution
        sort_idx = np.argsort(mjd)
        ax.plot(
            mjd[sort_idx], col[sort_idx],
            color=colour, alpha=0.35, linewidth=1.0, zorder=2,
        )

        ax.set_ylabel(f"{f1} − {f2}  (mag)", fontsize=10)
        #ax.invert_yaxis()
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.tick_params(axis="both", labelsize=9)

    axes[-1].set_xlabel("MJD", fontsize=10)

    if outfile:
        fig.savefig(outfile, dpi=150, bbox_inches="tight")
        print(f"Saved: {outfile}")
    else:
        plt.show()

    plt.close(fig)


# ── Entry point ───────────────────────────────────────────────────────────────

def main(filepaths):
    for fp in filepaths:
        fp = Path(fp)
        if not fp.exists():
            print(f"File not found: {fp}")
            continue
        print(f"Processing: {fp}")
        meta, data = parse_file(fp)
        print(
            f"  {meta['name']}  z={meta['redshift']}  "
            f"Filters found: {list(data.keys())}"
        )
        outfile = fp.with_name(fp.stem + "_colour_curves.pdf")
        plot_colour_curves(meta, data, COLOUR_PAIRS)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_colour_curves.py <file1.txt> [file2.txt ...]")
        sys.exit(1)
    main(sys.argv[1:])