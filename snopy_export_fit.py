"""
SNooPy fit curve exporter
=========================
Reads all .snpy files in SNPY_DIR, evaluates the EBV_model2 fitted lightcurve
for each band at regular phase points, and writes everything to a single CSV
(snopy_fit_curves.csv) for use by snia_3d_surface_plot.py.

Run this script once from the terminal:
    python snopy_export_fit.py

It does NOT need to be run inside iPython.
"""

import os
import glob
import pickle
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================

SNPY_DIR   = "snpy_txt_clean"   # directory containing your .snpy files
SNPY_GLOB  = "SN*.snpy"         # glob pattern — adjust if filenames differ
OUTPUT_CSV = "snopy_fit_curves.csv"

PHASE_MIN = -10    # days relative to t_max
PHASE_MAX =  70
N_POINTS  =  200   # sample points per band per SN

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def load_snpy_file(fp):
    """Load a .snpy file using pickle (works regardless of SNooPy version)."""
    with open(fp, "rb") as fh:
        s = pickle.load(fh)
    return s


def get_distmod(s):
    """Extract distance modulus from a fitted SNooPy object."""
    for attr in ["mu", "DM", "dm"]:
        if hasattr(s.model, attr):
            val = getattr(s.model, attr)
            if val is not None and np.isfinite(float(val)):
                return float(val)
    raise RuntimeError(
        f"Could not find distance modulus on s.model for {s.name}. "
        "Tried: mu, DM, dm."
    )


def eval_model_band(s, band, mjd_arr):
    """
    Evaluate the EBV_model2 model for a given band at MJD values mjd_arr.
    Returns array of model magnitudes (NaN where model is undefined).
    EBV_model2 stores the template via s.model.M0, evaluated through
    s.model.get_template() or directly via s.model.__call__.
    """
    try:
        # Primary method: model is callable per band
        mags, _ = s.model(band, mjd_arr)
        return np.array(mags, dtype=float)
    except Exception:
        pass

    try:
        # Fallback: use the lc object's model attribute if fitted
        lc = s.data[band]
        mags = lc.model(mjd_arr)
        return np.array(mags, dtype=float)
    except Exception:
        pass

    return np.full(len(mjd_arr), np.nan)


def main():
    filepaths = sorted(glob.glob(os.path.join(SNPY_DIR, SNPY_GLOB)))
    if not filepaths:
        raise FileNotFoundError(
            f"No files matching '{SNPY_GLOB}' found in '{SNPY_DIR}'."
        )
    print(f"Found {len(filepaths)} .snpy file(s).\n")

    phases  = np.linspace(PHASE_MIN, PHASE_MAX, N_POINTS)
    all_rows = []

    for fp in filepaths:
        fname   = os.path.basename(fp)
        sn_name = os.path.splitext(fname)[0]   # e.g. "SN2005cf"
        stem    = sn_name[2:] if sn_name.startswith("SN") else sn_name

        try:
            s = load_snpy_file(fp)

            dm15 = float(s.model.dm15)
            tmax = float(s.model.Tmax)
            dm   = get_distmod(s)

            mjd_arr = phases + tmax

            # Bands present in this SN's data
            bands = list(s.data.keys())
            sn_rows = 0

            for band in bands:
                mags = eval_model_band(s, band, mjd_arr)

                # Only keep phases where the model is defined
                valid = np.isfinite(mags)
                if not np.any(valid):
                    print(f"    Band {band}: no valid model points, skipping.")
                    continue

                for ph, mg in zip(phases[valid], mags[valid]):
                    all_rows.append({
                        "sn_name": stem,      # strip "SN" to match CSV stems
                        "band":    band,
                        "phase":   ph,
                        "mag":     mg,
                        "abs_mag": mg - dm,
                        "dm15":    dm15,
                        "tmax":    tmax,
                        "DM":      dm,
                    })
                sn_rows += np.sum(valid)

            print(f"  OK : {fname}  "
                  f"(dm15={dm15:.3f}, DM={dm:.3f}, {len(bands)} bands, "
                  f"{sn_rows} points)")

        except Exception as exc:
            print(f"  SKIPPED: {fname} — {exc}")

    if not all_rows:
        raise RuntimeError("No fit curve data exported. Check your .snpy files.")

    df = pd.DataFrame(all_rows)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nExported {len(df)} rows for {df['sn_name'].nunique()} SNe "
          f"across {df['band'].nunique()} bands → {OUTPUT_CSV}")
    print(f"Bands present: {sorted(df['band'].unique())}")


if __name__ == "__main__":
    main()
