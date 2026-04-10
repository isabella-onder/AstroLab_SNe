"""
SN Ia 3D lightcurve surface plot — Burns et al. (2011) style, magnitude space
==============================================================================
Reproduces Figure 4 of Burns et al. (2011, AJ, 141, 19).

Axes: (1) t - t_max   (2) dm15(B)   (3) Absolute magnitude M = mag - DM
      z-axis is inverted so bright epochs sit at the top, matching
      the visual convention of Burns et al. (2011).

Input:  - one CSV (.csv) file per SN, each containing all bands
        - a separate lookup file (.xlsx) with SN parameters including DM
Output: a single figure with one subplot for the chosen band

Surface fitting strategy
------------------------
Each SN's apparent magnitudes are converted to absolute magnitudes using the
distance modulus DM from the lookup file: M = mag - DM. This places all SNe
on a common absolute scale, removing the distance spread that would otherwise
cause the spline to produce unphysical wiggles. Both the surface and the
scatter are plotted in absolute magnitude, each SN in a distinct colour.

Usage
-----
    python snia_3d_surface_plot.py

Only edit the CONFIGURATION block below.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from scipy.interpolate import SmoothBivariateSpline

# =============================================================================
# CONFIGURATION — edit this block only
# =============================================================================

DATA_DIR     = "../Photometry_final_Ia"
FILE_PATTERN = "*.csv"

# -- Band to plot -------------------------------------------------------------
BAND = "I"                # Single filter label, exactly as it appears in BAND_COL

# -- Column names inside each lightcurve CSV file -----------------------------
MJD_COL  = "MJD"
MAG_COL  = "magnitude"
BAND_COL = "Filter"
ERR_COL  = "mag_error_up"   # Set to None to omit errorbars

# -- dm15 / t_max lookup file (.xlsx) -----------------------------------------
# Layout:
#   - Columns = SN names (e.g. "SN2005cf")
#   - Rows    = parameters, with a left-hand label column
LOOKUP_PATH      = "snpy_txt_clean/snoopy_results_EBV_fixed.xlsx"
LOOKUP_LABEL_COL = "Parameter"
DM15_ROW_LABEL   = "dm15"
TMAX_ROW_LABEL   = "T_max"
DM_ROW_LABEL     = "DM"       # row label for the distance modulus

# Prefix/suffix in lookup column headers not present in filenames:
# e.g. lookup "SN2005cf" <-> filename "2005cf.csv"
LOOKUP_SN_PREFIX = "SN"
LOOKUP_SN_SUFFIX = ""

# -- Phase range to display ---------------------------------------------------
PHASE_MIN, PHASE_MAX = -10, 70          # days relative to t_max

# -- dm15 axis limits ---------------------------------------------------------
DM15_MIN, DM15_MAX = 0.79, 1.3

# -- Surface smoothing --------------------------------------------------------
N_GRID_PHASE = 60
N_GRID_DM15  = 30
# s = SPLINE_SMOOTHING_FACTOR * n_points.
# Start at 1.0 and increase if the surface is still too wiggly.
SPLINE_SMOOTHING_FACTOR = 5.0

# -- Colour map for individual SNe --------------------------------------------
# Any matplotlib qualitative colormap works: "tab20", "tab10", "Set1", etc.
SN_COLORMAP = "tab20"

# -- Figure output ------------------------------------------------------------
OUTPUT_PATH = "snia_3d_surface.pdf"
FIG_WIDTH   = 8.0
FIG_HEIGHT  = 6.0

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def load_lookup(lookup_path, label_col, dm15_row, tmax_row, dm_row, prefix, suffix):
    """Return dm15_map, tmax_map, and dm_map dicts keyed by filename stem."""
    if lookup_path.endswith(".csv"):
        ldf = pd.read_csv(lookup_path, index_col=label_col)
    else:
        ldf = pd.read_excel(lookup_path, index_col=label_col)

    for row in (dm15_row, tmax_row, dm_row):
        if row not in ldf.index:
            raise KeyError(
                f"Row '{row}' not found in lookup file. "
                f"Available rows: {list(ldf.index)}"
            )

    def strip_name(col):
        name = str(col)
        if prefix and name.startswith(prefix):
            name = name[len(prefix):]
        if suffix and name.endswith(suffix):
            name = name[:-len(suffix)]
        return name

    dm15_map = {strip_name(c): float(ldf.loc[dm15_row, c]) for c in ldf.columns}
    tmax_map = {strip_name(c): float(ldf.loc[tmax_row, c]) for c in ldf.columns}
    dm_map   = {strip_name(c): float(ldf.loc[dm_row,   c]) for c in ldf.columns}
    return dm15_map, tmax_map, dm_map


def load_all_data(data_dir, file_pattern, dm15_map, tmax_map, dm_map, cfg):
    """
    Load all CSV lightcurve files. Returns a DataFrame with columns:
        phase, mag, abs_mag, dm15, err, band, sn_name
    abs_mag = mag - DM  (absolute magnitude, used for both surface and scatter).
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
                      f"'{sn_name}' not in lookup file.")
                continue

            dm15_val = dm15_map[sn_name]
            tmax_val = tmax_map[sn_name]
            dm_val   = dm_map[sn_name]

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
                print(f"  EMPTY  : {os.path.basename(fp)} — no data in phase range.")
                continue

            # Absolute magnitude: removes distance modulus spread
            df["abs_mag"] = df["mag"] - dm_val

            frames.append(df)
            print(f"  Loaded : {os.path.basename(fp)}"
                  f"  ({len(df)} pts, dm15={dm15_val:.3f}, "
                  f"DM={dm_val:.3f}, t_max={tmax_val:.2f})")

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
        raise KeyError(f"Column(s) not found: {missing}. Available: {list(df.columns)}")


def fit_surface(phase, dm15, rel_mag, n_phase, n_dm15, smoothing_factor):
    """
    Fit a smooth bivariate cubic spline over (phase, dm15) -> rel_mag.
    Using rel_mag (peak-normalised) prevents distance modulus spread from
    producing unphysical wiggles in the surface.
    """
    mask = np.isfinite(phase) & np.isfinite(dm15) & np.isfinite(rel_mag)
    phase, dm15, rel_mag = phase[mask], dm15[mask], rel_mag[mask]

    s      = smoothing_factor * len(phase)
    spline = SmoothBivariateSpline(phase, dm15, rel_mag, kx=3, ky=3, s=s)

    p_grid = np.linspace(phase.min(), phase.max(), n_phase)
    d_grid = np.linspace(dm15.min(),  dm15.max(),  n_dm15)
    P, D   = np.meshgrid(p_grid, d_grid)
    Z      = spline(p_grid, d_grid).T
    return P, D, Z


def make_figure(combined, band, cfg):
    """
    3D figure: both surface and scatter plotted on abs_mag = mag - DM,
    so all SNe are on a common absolute magnitude scale.
    z-axis is inverted so bright (more negative M) sits at the top.
    """
    sn_names   = sorted(combined["sn_name"].unique())
    n_sne      = len(sn_names)
    cmap       = cm.get_cmap(cfg["sn_colormap"], n_sne)
    colour_map = {sn: cmap(i) for i, sn in enumerate(sn_names)}

    fig = plt.figure(figsize=(cfg["fig_width"], cfg["fig_height"]))
    ax  = fig.add_subplot(111, projection="3d")

    # -- Smooth surface (fitted on abs_mag) -----------------------------------
    phase_all   = combined["phase"].values.astype(float)
    dm15_all    = combined["dm15"].values.astype(float)
    abs_mag_all = combined["abs_mag"].values.astype(float)

    if len(phase_all) >= 16:
        try:
            P, D, Z = fit_surface(
                phase_all, dm15_all, abs_mag_all,
                cfg["n_grid_phase"], cfg["n_grid_dm15"],
                cfg["spline_smoothing_factor"]
            )
            ax.plot_surface(P, D, Z,
                            rstride=1, cstride=1,
                            color="silver", alpha=0.45,
                            linewidth=0, antialiased=True, shade=True)
        except Exception as exc:
            print(f"  Warning: surface fit failed: {exc}")
    else:
        print(f"  Warning: only {len(phase_all)} points — surface fit skipped.")

    # -- Per-SN coloured scatter (absolute magnitudes) ------------------------
    for sn, grp in combined.groupby("sn_name"):
        phase   = grp["phase"].values.astype(float)
        dm15    = grp["dm15"].values.astype(float)
        abs_mag = grp["abs_mag"].values.astype(float)
        err     = grp["err"].values.astype(float)
        col     = colour_map[sn]

        if np.any(np.isfinite(err)):
            for x, y, z, e in zip(phase, dm15, abs_mag, err):
                if np.isfinite(e):
                    ax.plot([x, x], [y, y], [z - e, z + e],
                            color=col, linewidth=0.5, alpha=0.5)

        ax.scatter(phase, dm15, abs_mag,
                   s=8, color=col, depthshade=True,
                   alpha=0.85, linewidths=0, label=sn)

    # -- Legend ---------------------------------------------------------------
    if n_sne <= 30:
        ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1),
                  fontsize=6, framealpha=0.7,
                  markerscale=1.5, title="SN name", title_fontsize=7)

    # -- Labels & appearance --------------------------------------------------
    ax.set_xlabel(r"$t - t_{\mathrm{max}}$ (days)", labelpad=8,  fontsize=10)
    ax.set_ylabel(r"$\Delta m_{15}(B)$",             labelpad=8,  fontsize=10)
    ax.set_zlabel(r"$M$ (absolute magnitude)",       labelpad=10, fontsize=10)
    ax.set_title(f"{band} band", fontsize=11, pad=4)

    ax.set_xlim(cfg["phase_min"], cfg["phase_max"])
    ax.set_ylim(cfg["dm15_min"],  cfg["dm15_max"])

    # Invert z: bright (more negative M) at top
    all_abs    = combined["abs_mag"].values
    z_top      = np.nanmin(all_abs) - 0.3
    z_bottom   = np.nanmax(all_abs) + 0.3
    ax.set_zlim(z_bottom, z_top)

    ax.tick_params(axis="both", labelsize=8)
    ax.view_init(elev=20, azim=-60)

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
        sn_colormap             = SN_COLORMAP,
        output_path             = OUTPUT_PATH,
        fig_width               = FIG_WIDTH,
        fig_height              = FIG_HEIGHT,
    )

    print(f"Loading parameter lookup from: {LOOKUP_PATH}\n")
    dm15_map, tmax_map, dm_map = load_lookup(
        LOOKUP_PATH, LOOKUP_LABEL_COL,
        DM15_ROW_LABEL, TMAX_ROW_LABEL, DM_ROW_LABEL,
        LOOKUP_SN_PREFIX, LOOKUP_SN_SUFFIX
    )
    print(f"Found {len(dm15_map)} SNe in lookup file.\n")

    combined = load_all_data(DATA_DIR, FILE_PATTERN, dm15_map, tmax_map, dm_map, cfg)

    available_bands = combined["band"].unique()
    if BAND not in available_bands:
        raise RuntimeError(
            f"Band '{BAND}' not found in the loaded data. "
            f"Available bands: {sorted(available_bands)}"
        )

    print(f"Plotting band '{BAND}' for {combined['sn_name'].nunique()} SNe Ia...\n")
    make_figure(combined, BAND, cfg)


if __name__ == "__main__":
    main()
