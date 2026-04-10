"""
SN Ia 3D lightcurve surface plot — Burns et al. (2011) style, magnitude space
==============================================================================
Reproduces Figure 4 of Burns et al. (2011, AJ, 141, 19).

Axes: (1) t - t_max   (2) dm15(B)   (3) Apparent magnitude
      z-axis is inverted so bright epochs sit at the top, matching
      the visual convention of Burns et al. (2011).

Input:  - one Excel (.xlsx) file per SN, each containing all bands
        - a separate lookup file (xlsx or csv) with SN parameters
Output: a single figure with one subplot for the chosen band

Usage
-----
    python snia_3d_surface_plot_xlsx.py

Only edit the CONFIGURATION block below.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from scipy.interpolate import SmoothBivariateSpline

# =============================================================================
# CONFIGURATION — edit this block only
# =============================================================================

DATA_DIR     = "../Photometry_final"        # Directory containing the per-SN Excel files
FILE_PATTERN = "*.csv"   # Glob pattern to locate your lightcurve files

# -- Band to plot -------------------------------------------------------------
BAND = "R"                # Single filter label, exactly as it appears in BAND_COL

# -- Column names inside each lightcurve Excel file ---------------------------
MJD_COL  = "MJD"         # Modified Julian Date of each observation
MAG_COL  = "magnitude"   # Apparent magnitude
BAND_COL = "Filter"      # Photometric band label
ERR_COL  = "mag_error_up"       # Photometric uncertainty; set to None to omit errorbars

# -- dm15 / t_max lookup file -------------------------------------------------
# Layout expected:
#   - Columns = SN names (e.g. "SN2005cf", "SN2006X", ...)
#   - Rows    = parameters, identified by a label column on the left
#   - The row where the label column reads DM15_ROW_LABEL contains dm15
#   - The row where the label column reads TMAX_ROW_LABEL contains t_max (MJD)
#
LOOKUP_PATH      = "snpy_txt_clean/snoopy_results_EBV_fixed.xlsx"  # path to the lookup file (.xlsx or .csv)
LOOKUP_LABEL_COL = "Parameter"           # name of the left-hand label column
DM15_ROW_LABEL   = "dm15"               # value in LOOKUP_LABEL_COL for the dm15 row
TMAX_ROW_LABEL   = "T_max"               # value in LOOKUP_LABEL_COL for the t_max row

# -- SN name prefix/suffix in the lookup file ---------------------------------
# The lookup column headers include a prefix not present in the filenames.
# e.g. lookup header "SN2005cf"  <->  filename "2005cf.xlsx"
# Set LOOKUP_SN_PREFIX = "SN" and LOOKUP_SN_SUFFIX = ""
# If there is no prefix or suffix, set them to "".
LOOKUP_SN_PREFIX = "SN"
LOOKUP_SN_SUFFIX = ""

# -- Phase range to display ---------------------------------------------------
PHASE_MIN, PHASE_MAX = -10, 70          # days relative to t_max

# -- dm15 axis limits ---------------------------------------------------------
DM15_MIN, DM15_MAX = 0.6, 1.9          # mag; widen if your sample demands it

# -- Smoothing surface --------------------------------------------------------
N_GRID_PHASE = 60
N_GRID_DM15  = 30
# Smoothing factor for SmoothBivariateSpline (s = factor * n_points).
# Increase → smoother surface; decrease → tighter fit to data.
SPLINE_SMOOTHING_FACTOR = 0.5

# -- Figure output ------------------------------------------------------------
OUTPUT_PATH    = "snia_3d_surface.pdf"  # extension sets the format
FIG_WIDTH      = 7.0    # inches
FIG_HEIGHT     = 5.5    # inches

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def load_lookup(lookup_path, label_col, dm15_row, tmax_row, prefix, suffix):
    """
    Load the parameter lookup file and return two dicts keyed by filename stem:
        dm15_map  : { "2005cf": 1.05, ... }
        tmax_map  : { "2005cf": 53549.2, ... }

    The lookup file has SN names as column headers (e.g. "SN2005cf").
    The prefix/suffix are stripped to obtain the filename stem.
    """
    if lookup_path.endswith(".csv"):
        ldf = pd.read_csv(lookup_path, index_col=label_col)
    else:
        ldf = pd.read_excel(lookup_path, index_col=label_col)

    if dm15_row not in ldf.index:
        raise KeyError(
            f"Row '{dm15_row}' not found in lookup file. "
            f"Available rows: {list(ldf.index)}"
        )
    if tmax_row not in ldf.index:
        raise KeyError(
            f"Row '{tmax_row}' not found in lookup file. "
            f"Available rows: {list(ldf.index)}"
        )

    def strip_name(col):
        """Remove configured prefix and suffix from a column header."""
        name = str(col)
        if prefix and name.startswith(prefix):
            name = name[len(prefix):]
        if suffix and name.endswith(suffix):
            name = name[:-len(suffix)]
        return name

    dm15_map = {strip_name(col): float(ldf.loc[dm15_row, col])
                for col in ldf.columns}
    tmax_map = {strip_name(col): float(ldf.loc[tmax_row, col])
                for col in ldf.columns}

    return dm15_map, tmax_map


def load_all_data(data_dir, file_pattern, dm15_map, tmax_map, cfg):
    """
    Load all Excel lightcurve files and return a concatenated DataFrame with
    columns: phase, mag, dm15, err, band, sn_name.
    phase is computed as MJD - t_max; magnitudes are kept as-is.
    Only the band specified in cfg['band'] is retained.
    """
    filepaths = sorted(glob.glob(os.path.join(data_dir, file_pattern)))
    if not filepaths:
        raise FileNotFoundError(
            f"No files matching '{file_pattern}' found in '{data_dir}'.\n"
            "Check DATA_DIR and FILE_PATTERN in the CONFIGURATION block."
        )
    print(f"Found {len(filepaths)} lightcurve file(s).\n")

    frames = []
    for fp in filepaths:
        sn_name = os.path.splitext(os.path.basename(fp))[0]
        try:
            if sn_name not in dm15_map:
                print(f"  SKIPPED: {os.path.basename(fp)} — "
                      f"'{sn_name}' not found in lookup file.")
                continue
            if sn_name not in tmax_map:
                print(f"  SKIPPED: {os.path.basename(fp)} — "
                      f"t_max for '{sn_name}' not found in lookup file.")
                continue

            dm15_val = dm15_map[sn_name]
            tmax_val = tmax_map[sn_name]

            df = pd.read_csv(fp)
            _validate_columns(df, fp, cfg)

            df = df.rename(columns={
                cfg["mjd_col"]:  "mjd",
                cfg["mag_col"]:  "mag",
                cfg["band_col"]: "band",
            })

            if cfg["err_col"] and cfg["err_col"] in df.columns:
                df = df.rename(columns={cfg["err_col"]: "err"})
            else:
                df["err"] = np.nan

            # Retain only the requested band
            df = df[df["band"] == cfg["band"]].copy()
            if df.empty:
                print(f"  SKIPPED: {os.path.basename(fp)} — "
                      f"no data for band '{cfg['band']}'.")
                continue

            df["phase"]   = df["mjd"] - tmax_val
            df["dm15"]    = dm15_val
            df["sn_name"] = sn_name

            df = df[(df["phase"] >= cfg["phase_min"]) &
                    (df["phase"] <= cfg["phase_max"])].copy()

            if df.empty:
                print(f"  EMPTY  : {os.path.basename(fp)} — "
                      "no data in phase range.")
            else:
                frames.append(df)
                print(f"  Loaded : {os.path.basename(fp)}"
                      f"  ({len(df)} points, "
                      f"dm15 = {dm15_val:.3f}, "
                      f"t_max = {tmax_val:.2f})")

        except Exception as exc:
            print(f"  SKIPPED: {os.path.basename(fp)} — {exc}")

    if not frames:
        raise RuntimeError("No data loaded. Check the CONFIGURATION block.")

    combined = pd.concat(frames, ignore_index=True)
    print(f"\nLoaded {combined['sn_name'].nunique()} SNe Ia, "
          f"{len(combined)} photometric points total.\n")
    return combined


def _validate_columns(df, filepath, cfg):
    required = [cfg["mjd_col"], cfg["mag_col"], cfg["band_col"]]
    missing  = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"Column(s) not found: {missing}. "
            f"Available: {list(df.columns)}"
        )


def fit_surface(phase, dm15, mag, n_phase, n_dm15, smoothing_factor):
    """Fit a smooth bivariate cubic spline over (phase, dm15) -> mag."""
    mask = np.isfinite(phase) & np.isfinite(dm15) & np.isfinite(mag)
    phase, dm15, mag = phase[mask], dm15[mask], mag[mask]

    s      = smoothing_factor * len(phase)
    spline = SmoothBivariateSpline(phase, dm15, mag, kx=3, ky=3, s=s)

    p_grid = np.linspace(phase.min(), phase.max(), n_phase)
    d_grid = np.linspace(dm15.min(),  dm15.max(),  n_dm15)
    P, D   = np.meshgrid(p_grid, d_grid)
    Z      = spline(p_grid, d_grid).T   # shape (n_dm15, n_phase)
    return P, D, Z


def make_figure(combined, band, cfg):
    """
    Create and save the 3D figure for the chosen band.
    z-axis uses raw apparent magnitudes, inverted so bright is at the top.
    """
    phase = combined["phase"].values.astype(float)
    dm15  = combined["dm15"].values.astype(float)
    mag   = combined["mag"].values.astype(float)
    err   = combined["err"].values.astype(float)

    fig = plt.figure(figsize=(cfg["fig_width"], cfg["fig_height"]))
    ax  = fig.add_subplot(111, projection="3d")

    # -- Smooth surface -------------------------------------------------------
    if len(phase) >= 16:
        try:
            P, D, Z = fit_surface(
                phase, dm15, mag,
                cfg["n_grid_phase"], cfg["n_grid_dm15"],
                cfg["spline_smoothing_factor"]
            )
            ax.plot_surface(P, D, Z,
                            rstride=1, cstride=1,
                            color="silver", alpha=0.55,
                            linewidth=0, antialiased=True, shade=True)
        except Exception as exc:
            print(f"  Warning: surface fit failed: {exc}")
    else:
        print(f"  Warning: only {len(phase)} points — "
              "surface fit skipped (need ≥ 16).")

    # -- Errorbars ------------------------------------------------------------
    if np.any(np.isfinite(err)):
        for x, y, z, e in zip(phase, dm15, mag, err):
            if np.isfinite(e):
                ax.plot([x, x], [y, y], [z - e, z + e],
                        color="black", linewidth=0.4, alpha=0.4)

    # -- Data scatter ---------------------------------------------------------
    ax.scatter(phase, dm15, mag,
               s=6,c="black" , depthshade=True, alpha=0.7, linewidths=0)

    # -- Labels & appearance --------------------------------------------------
    ax.set_xlabel(r"$t - t_{\mathrm{max}}$ (days)", labelpad=8,  fontsize=10)
    ax.set_ylabel(r"$\Delta m_{15}(B)$",             labelpad=8,  fontsize=10)
    ax.set_zlabel("Apparent magnitude",              labelpad=10, fontsize=10)
    ax.set_title(f"{band} band", fontsize=11, pad=4)

    ax.set_xlim(cfg["phase_min"], cfg["phase_max"])
    ax.set_ylim(cfg["dm15_min"],  cfg["dm15_max"])

    # Invert z: bright (small mag) at top, faint (large mag) at bottom
    mag_bright = np.nanmin(mag) - 0.3
    mag_faint  = np.nanmax(mag) + 0.3
    ax.set_zlim(mag_faint, mag_bright)

    ax.tick_params(axis="both", labelsize=8)
    ax.view_init(elev=20, azim=-60)     # viewpoint matching Burns+2011

    plt.tight_layout()
    fig.savefig(cfg["output_path"], dpi=200, bbox_inches="tight")
    print(f"Figure saved to: {cfg['output_path']}")
    plt.show()


def main():
    cfg = dict(
        band                    = BAND,
        mjd_col                 = MJD_COL,
        mag_col                 = MAG_COL,
        band_col                = BAND_COL,
        err_col                 = ERR_COL,
        phase_min               = PHASE_MIN,
        phase_max               = PHASE_MAX,
        dm15_min                = DM15_MIN,
        dm15_max                = DM15_MAX,
        n_grid_phase            = N_GRID_PHASE,
        n_grid_dm15             = N_GRID_DM15,
        spline_smoothing_factor = SPLINE_SMOOTHING_FACTOR,
        output_path             = OUTPUT_PATH,
        fig_width               = FIG_WIDTH,
        fig_height              = FIG_HEIGHT,
    )

    # -- Load dm15 and t_max from the lookup file -----------------------------
    print(f"Loading parameter lookup from: {LOOKUP_PATH}\n")
    dm15_map, tmax_map = load_lookup(
        LOOKUP_PATH, LOOKUP_LABEL_COL,
        DM15_ROW_LABEL, TMAX_ROW_LABEL,
        LOOKUP_SN_PREFIX, LOOKUP_SN_SUFFIX
    )
    print(f"Found {len(dm15_map)} SNe in lookup file.\n")

    # -- Load lightcurve data -------------------------------------------------
    combined = load_all_data(DATA_DIR, FILE_PATTERN, dm15_map, tmax_map, cfg)

    # -- Check the requested band is present ----------------------------------
    available_bands = combined["band"].unique()
    if BAND not in available_bands:
        raise RuntimeError(
            f"Band '{BAND}' not found in the loaded data. "
            f"Available bands: {sorted(available_bands)}"
        )

    print(f"Plotting band '{BAND}' for "
          f"{combined['sn_name'].nunique()} SNe Ia...\n")
    make_figure(combined, BAND, cfg)


if __name__ == "__main__":
    main()