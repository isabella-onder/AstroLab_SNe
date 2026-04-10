"""
SN Ia 3D lightcurve surface plot — Burns et al. (2011) style, magnitude space
==============================================================================
Reads per-band .dat files of the form:
    SN20xxaa_lc_<band>_data.dat   — observed photometry
    SN20xxaa_lc_<band>_model.dat  — SNooPy model curve
Both are plotted on the same 3D axes, distinguished by marker style.
Apparent magnitudes are converted to absolute magnitudes: M = m - DM.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import SmoothBivariateSpline

# =============================================================================
# CONFIGURATIONe
# =============================================================================

DATA_DIRS = ['.']

BAND = "Rs"

# dm15 / t_max / DM lookup (.xlsx) -------------------------------------------
LOOKUP_PATH      = "../snpy_txt_clean/snoopy_results_EBV_fixed.xlsx"
LOOKUP_LABEL_COL = "Parameter"
DM15_ROW_LABEL   = "dm15"
TMAX_ROW_LABEL   = "T_max"
DM_ROW_LABEL     = "DM"

LOOKUP_SN_PREFIX = "SN"
LOOKUP_SN_SUFFIX = ""

# Phase range to display (days relative to t_max) ----------------------------
PHASE_MIN, PHASE_MAX = -10, 70

# Surface smoothing (fitted to model points) ----------------------------------
N_GRID_PHASE            = 60
N_GRID_DM15             = 30
SPLINE_SMOOTHING_FACTOR = 5.0

# Appearance ------------------------------------------------------------------
SN_COLORMAP = "tab20"

# Marker styles for data vs model
DATA_MARKER  = 'o'
MODEL_MARKER = 's'      # square
DATA_SIZE    = 18
MODEL_SIZE   = 8

# Figure output ---------------------------------------------------------------
OUTPUT_PATH = "snia_3d_surface_merged.pdf"
FIG_WIDTH   = 9.0
FIG_HEIGHT  = 6.5

# Padding around data extent for axis limits ----------------------------------
PAD_PHASE = 2.0
PAD_DM15  = 0.03
PAD_MAG   = 0.2

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def parse_dat_file(filepath):
    """
    Read a SNooPy .dat file (data or model).
    Returns arrays (t, mag, mag_err).
    """
    rows = []
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.replace(',', ' ').split()
            rows.append([float(x) for x in parts[:3]])
    if not rows:
        raise ValueError(f"No data rows found in {filepath}")
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def load_lookup(lookup_path, label_col, dm15_row, tmax_row, dm_row, prefix, suffix):
    """Return dm15_map, tmax_map, dm_map dicts keyed by filename stem."""
    if lookup_path.endswith('.csv'):
        ldf = pd.read_csv(lookup_path, index_col=label_col)
    else:
        ldf = pd.read_excel(lookup_path, index_col=label_col)

    for row in (dm15_row, tmax_row, dm_row):
        if row not in ldf.index:
            raise KeyError(
                f"Row '{row}' not found in lookup. "
                f"Available: {list(ldf.index)}"
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


def load_files(data_dirs, band, suffix, dm15_map, tmax_map, dm_map, kind):
    """
    Scan each directory for files matching SN*_lc_<band>_<suffix>.dat.
    Returns a DataFrame with columns:
        phase, abs_mag, mag_err, dm15, sn_name, kind
    `kind` is a string label: 'data' or 'model'.
    """
    pattern = f"*_lc_{band}_{suffix}.dat"
    frames  = []

    for data_dir in data_dirs:
        filepaths = sorted(glob.glob(os.path.join(data_dir, pattern)))
        if not filepaths:
            print(f"  [{kind}] No files matching '{pattern}' in '{data_dir}'.")
            continue

        for fp in filepaths:
            basename = os.path.basename(fp)
            sn_full  = basename.split('_lc_')[0]
            sn_name  = sn_full
            if LOOKUP_SN_PREFIX and sn_name.startswith(LOOKUP_SN_PREFIX):
                sn_name = sn_name[len(LOOKUP_SN_PREFIX):]
            if LOOKUP_SN_SUFFIX and sn_name.endswith(LOOKUP_SN_SUFFIX):
                sn_name = sn_name[:-len(LOOKUP_SN_SUFFIX)]

            if sn_name not in dm15_map:
                print(f"  [{kind}] SKIPPED: {basename} — '{sn_name}' not in lookup.")
                continue

            try:
                t, mag, mag_err = parse_dat_file(fp)
            except Exception as exc:
                print(f"  [{kind}] SKIPPED: {basename} — {exc}")
                continue

            tmax    = tmax_map[sn_name]
            dm15    = dm15_map[sn_name]
            mu      = dm_map[sn_name]
            phase   = t - tmax
            abs_mag = mag - mu

            mask = (phase >= PHASE_MIN) & (phase <= PHASE_MAX)
            if not np.any(mask):
                print(f"  [{kind}] EMPTY  : {basename} — no data in phase range.")
                continue

            df = pd.DataFrame({
                'phase':   phase[mask],
                'abs_mag': abs_mag[mask],
                'mag_err': mag_err[mask],
                'dm15':    dm15,
                'sn_name': sn_name,
                'kind':    kind,
            })
            frames.append(df)
            print(f"  [{kind}] Loaded : {basename}  "
                  f"({mask.sum()} pts, dm15={dm15:.3f}, DM={mu:.3f})")

    if not frames:
        return pd.DataFrame(
            columns=['phase', 'abs_mag', 'mag_err', 'dm15', 'sn_name', 'kind']
        )
    return pd.concat(frames, ignore_index=True)


def load_all_data(data_dirs, band, dm15_map, tmax_map, dm_map):
    """Load both data and model files; return a single combined DataFrame."""
    df_data  = load_files(data_dirs, band, 'data',  dm15_map, tmax_map, dm_map, 'data')
    df_model = load_files(data_dirs, band, 'model', dm15_map, tmax_map, dm_map, 'model')

    combined = pd.concat([df_data, df_model], ignore_index=True)
    if combined.empty:
        raise RuntimeError("No data loaded. Check DATA_DIRS, BAND, and LOOKUP_PATH.")

    n_sne = combined['sn_name'].nunique()
    print(f"\nLoaded {n_sne} unique SNe, {len(combined)} points total "
          f"({len(df_data)} data, {len(df_model)} model).\n")
    return combined


def fit_surface(phase, dm15, abs_mag, n_phase, n_dm15, smoothing_factor):
    mask = np.isfinite(phase) & np.isfinite(dm15) & np.isfinite(abs_mag)
    phase, dm15, abs_mag = phase[mask], dm15[mask], abs_mag[mask]
    s      = smoothing_factor * len(phase)
    spline = SmoothBivariateSpline(phase, dm15, abs_mag, kx=3, ky=3, s=s)
    p_grid = np.linspace(phase.min(), phase.max(), n_phase)
    d_grid = np.linspace(dm15.min(),  dm15.max(),  n_dm15)
    P, D   = np.meshgrid(p_grid, d_grid)
    Z      = spline(p_grid, d_grid).T
    return P, D, Z


def make_figure(combined, band):
    sn_names   = sorted(combined['sn_name'].unique())
    n_sne      = len(sn_names)
    cmap       = cm.get_cmap(SN_COLORMAP, n_sne)
    colour_map = {sn: cmap(i) for i, sn in enumerate(sn_names)}

    # Axis limits from combined data -----------------------------------------
    p_min    = combined['phase'].min()   - PAD_PHASE
    p_max    = combined['phase'].max()   + PAD_PHASE
    d_min    = combined['dm15'].min()    - PAD_DM15
    d_max    = combined['dm15'].max()    + PAD_DM15
    z_top    = combined['abs_mag'].min() - PAD_MAG
    z_bottom = combined['abs_mag'].max() + PAD_MAG

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    ax  = fig.add_subplot(111, projection='3d')

    # --- Optional smooth surface fitted to model points ---------------------
    df_model = combined[combined['kind'] == 'model']
    if len(df_model) >= 16:
        try:
            P, D, Z = fit_surface(
                df_model['phase'].values.astype(float),
                df_model['dm15'].values.astype(float),
                df_model['abs_mag'].values.astype(float),
                N_GRID_PHASE, N_GRID_DM15, SPLINE_SMOOTHING_FACTOR,
            )
            # Uncomment to enable the surface:
            # ax.plot_surface(P, D, Z, rstride=1, cstride=1,
            #                 color='silver', alpha=0.35,
            #                 linewidth=0, antialiased=True, shade=True)
        except Exception as exc:
            print(f"  Warning: surface fit failed — {exc}")

    # --- Per-SN scatter: data points (circles) and model (squares) ----------
    legend_handles = []

    for sn, grp in combined.groupby('sn_name'):
        col = colour_map[sn]

        grp_data  = grp[grp['kind'] == 'data']
        grp_model = grp[grp['kind'] == 'model']

        # Observed data — larger circles with black outline
        if not grp_data.empty:
            sc = ax.scatter(
                grp_data['phase'], grp_data['dm15'], grp_data['abs_mag'],
                s=DATA_SIZE, color=col, marker=DATA_MARKER,
                depthshade=False, alpha=1.0,
                linewidths=0.4, edgecolors='k',
                zorder=10,
            )
            # Errorbars for data
            for _, row in grp_data.iterrows():
                e = row['mag_err']
                if np.isfinite(e) and e > 0:
                    ax.plot(
                        [row['phase'],   row['phase']],
                        [row['dm15'],    row['dm15']],
                        [row['abs_mag'] - e, row['abs_mag'] + e],
                        color=col, linewidth=0.8, alpha=0.7,
                    )

        # Model curve — smaller squares, no outline, connected by line
        if not grp_model.empty:
            grp_model_sorted = grp_model.sort_values('phase')
            ax.scatter(
                grp_model_sorted['phase'],
                grp_model_sorted['dm15'],
                grp_model_sorted['abs_mag'],
                s=MODEL_SIZE, color=col, marker=MODEL_MARKER,
                depthshade=False, alpha=0.7,
                linewidths=0, zorder=9,
            )
            # Connect model points with a thin line for each SN
            ax.plot(
                grp_model_sorted['phase'].values,
                grp_model_sorted['dm15'].values,
                grp_model_sorted['abs_mag'].values,
                color=col, linewidth=0.8, alpha=0.5, zorder=8,
            )

        # One legend entry per SN (use data scatter handle if available)
        if not grp_data.empty or not grp_model.empty:
            legend_handles.append(
                plt.Line2D([0], [0], marker=DATA_MARKER, color='w',
                           markerfacecolor=col, markeredgecolor='k',
                           markeredgewidth=0.4, markersize=5, label=sn)
            )

    # --- Proxy handles for data vs model style legend -----------------------
    style_handles = [
        plt.Line2D([0], [0], marker=DATA_MARKER, color='w',
                   markerfacecolor='grey', markeredgecolor='k',
                   markeredgewidth=0.4, markersize=6, label='Observed data'),
        plt.Line2D([0], [0], marker=MODEL_MARKER, color='grey',
                   markersize=5, markerfacecolor='grey',
                   linewidth=0.8, label='SNooPy model'),
        #plt.Line2D([0], [0], marker='D', color='w',
                  # markerfacecolor='grey', markeredgecolor='k',
                  # markeredgewidth=0.5, markersize=6, label='Peak magnitude'),
    ]

    # --- Peak brightness points in 3D space ---------------------------------
    # For each SN, find the epoch of peak brightness (minimum absolute mag)
    # and plot it at its true (phase, dm15, abs_mag) coordinates with
    # uncertainties on the magnitude axis.

    for sn, grp in combined.groupby('sn_name'):
        col = colour_map[sn]

        idx_max  = grp['abs_mag'].idxmin()
        mag_peak = grp.loc[idx_max, 'abs_mag']
        err_peak = grp.loc[idx_max, 'mag_err']

        ax.scatter(p_max, d_max, mag_peak,
                   s=30, color=col, depthshade=False,
                   alpha=1.0, 
                   marker='o')

        #if np.isfinite(err_peak) and err_peak > 0:
        #    ax.plot([p_max, p_max],
        #            [d_max, d_max],
        #            [mag_peak - err_peak, mag_peak + err_peak],
        #            color=col, linewidth=1.2, alpha=0.9, zorder=11)

    # Legend: SN names on the right, style legend inside
    if n_sne <= 30:
        leg1 = ax.legend(
            handles=legend_handles,
            loc='center left', #bbox_to_anchor=(1.05, 1),
            fontsize=6, framealpha=0,
            markerscale=1.2, title='SN', #title_fontsize=7,
        )
        ax.add_artist(leg1)

    ax.legend(
        handles=style_handles,
        loc='upper left', fontsize=8, framealpha=0,
    )

    # --- Labels & limits ----------------------------------------------------
    ax.set_xlabel(r'$t - t_{\mathrm{max}}$ (days)', labelpad=8,  fontsize=10)
    ax.set_ylabel(r'$\Delta m_{15}(B)$',             labelpad=8,  fontsize=10)
    ax.set_zlabel(r'$M$ (absolute magnitude)',        labelpad=10, fontsize=10)
    ax.set_title(f'{band} band — data & model', fontsize=11, pad=4)

    ax.set_xlim(p_min, p_max)
    ax.set_ylim(d_min, d_max)
    ax.set_zlim(z_bottom, z_top)   # inverted: bright at top

    ax.tick_params(axis='both', labelsize=8)
    ax.view_init(elev=20, azim=-35)

    plt.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=200, bbox_inches='tight')
    print(f"Figure saved to: {OUTPUT_PATH}")
    plt.show()


def main():
    print(f"Loading parameter lookup from: {LOOKUP_PATH}\n")
    dm15_map, tmax_map, dm_map = load_lookup(
        LOOKUP_PATH, LOOKUP_LABEL_COL,
        DM15_ROW_LABEL, TMAX_ROW_LABEL, DM_ROW_LABEL,
        LOOKUP_SN_PREFIX, LOOKUP_SN_SUFFIX,
    )
    print(f"Found {len(dm15_map)} SNe in lookup.\n")

    combined = load_all_data(DATA_DIRS, BAND, dm15_map, tmax_map, dm_map)

    print(f"Plotting band '{BAND}'...\n")
    make_figure(combined, BAND)


if __name__ == '__main__':
    main()