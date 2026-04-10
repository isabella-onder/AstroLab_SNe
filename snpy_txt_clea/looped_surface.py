# =============================================================================
# SNooPy fit curve exporter
# =============================================================================
# Paste this into your iPython session AFTER fitting each supernova with
# EBV_model2. It samples the model lightcurve at regular phase points and
# appends a row per (phase, band) to a master CSV file.
#
# Assumes your SNooPy supernova object is called `s` and you have already run:
#     s.choose_model('EBV_model2')
#     s.fit(...)
#
# Usage — paste once per SN, or wrap in your existing loop:
#
#     exec(open('snopy_export_fit.py').read())
#
# =============================================================================

import numpy as np
import pandas as pd
import os

# -- Configuration ------------------------------------------------------------
PHASE_MIN   = -10      # days relative to t_max
PHASE_MAX   =  70
N_POINTS    =  200     # number of sample points along the curve
OUTPUT_CSV  = "snopy_fit_curves.csv"   # appended each time; safe to re-run
# -----------------------------------------------------------------------------

sn_name = s.name          # taken from the SNooPy object
dm15    = s.model.dm15    # SNooPy EBV_model2 fitted dm15
tmax    = s.model.Tmax    # fitted t_max (MJD)

SN_list = ['2025acsn', '2019np', '2020ue', '2025ahkt',
            '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc',
            '2026kc_subbed', '2025acsn_subbed']

EBV_host_list = [0.1586433260393873, 0.020058351568198397, 0.028446389496717725, 0.07221006564551423, 0.38439095550692925, 0.025893508388037927, 0.016776075857038657, 0.04412837345003647, 0.03646973012399708, 0.05397520058351568, 0.3931436907366886, 0.03391684901531729, 0.01276440554339898, 0.03391684901531729, 0.1586433260393873]   # or e.g. [0.1, 0.2, None, 0.05, ...]
the_list = [SN+'.txt' for SN in SN_list]
print(the_list)
all_results = {}
fout   = open('failures.log', 'w')

# ── Main loop (CSP Example 2 pattern) ────────────────────────────────────────
for file in the_list:

    try:
        s = get_sn(file)
        s.choose_model('EBV_model2', stype='dm15')

        if EBV_host_list is not None and EBV_host_list[idx] is not None:
            ebv = EBV_host_list[idx]
            print(f'  Fitting {s.name} with EBVhost = {ebv} (fixed prior)')
            s.fit(EBVhost=ebv)
        else:
            print(f'  Fitting {s.name} with EBVhost free')
            s.fit()
    except:
        fout.write('Failed on object %s\n' % s.name)
        all_results[s.name] = {p: None for p in PARAMETERS}
        
    # DM from SNooPy: mu = distance modulus
    try:
        dm = s.model.DM
    except AttributeError:
        # fallback if mu not directly on model
        dm = s.get_distmod()[0]

    phases = np.linspace(PHASE_MIN, PHASE_MAX, N_POINTS)

    rows = []
    for band in s.model.bands:
        try:
            # EBV_model2: model.eval(band, phases) returns model magnitudes
            mags = s.model.eval(band, phases + tmax)
            for ph, mg in zip(phases, mags):
                rows.append({
                    "sn_name": sn_name,
                    "band":    band,
                    "phase":   ph,
                    "mag":     mg,
                    "dm15":    dm15,
                    "tmax":    tmax,
                    "DM":      dm,
                })
        except Exception as e:
            print(f"  Could not evaluate band {band} for {sn_name}: {e}")

    if rows:
        df_new = pd.DataFrame(rows)
        write_header = not os.path.exists(OUTPUT_CSV)
        df_new.to_csv(OUTPUT_CSV, mode="a", header=write_header, index=False)
        print(f"Exported {len(rows)} points for {sn_name} "
            f"(dm15={dm15:.3f}, DM={dm:.3f}) → {OUTPUT_CSV}")
    else:
        print(f"No points exported for {sn_name}.")
