import pandas as pd
import numpy as np

SN_list = ['2019np','2020ue','2025ahsa','2026acd','2026atb','2026gm']
SN_list_subbed= ['2026kc_subbed']

SN_list_nonIa = ['2025aebt','2025afwg', '2026dix']
SNe = ['2025advo']
for SN in SNe:

    # --- Load Excel file ---
    df = pd.read_csv(f"../Photometry_final/{SN}.csv")

    # --- Rename columns (adjust if names differ slightly) ---
    df = df.rename(columns={
        'MJD': 'time',
        'magnitude': 'mag',
        'mag_err_upper': 'magerr',
        'Filter': 'filter',
        'zpt_mean': 'zp'
    })

    # --- Map Johnson–Cousins filters to sncosmo names ---
    filter_map = {
        'I': 'besselli',
        'V': 'bessellv',
        'B': 'bessellb',
        'R': 'bessellr'
    }

    df['band'] = df['filter'].map(filter_map)

    # --- Check for unmapped filters (important sanity check) ---
    if df['band'].isnull().any():
        print("WARNING: Some filters were not recognized:")
        print(df[df['band'].isnull()]['filter'].unique())

    # --- Convert magnitude to flux ---
    df['flux'] = 10**(-0.4 * (df['mag'] - df['zp']))

    # --- Convert magnitude error to flux error ---
    df['fluxerr'] = df['flux'] * (np.log(10) / 2.5) * df['magerr']

    # --- Set zeropoint system ---
    df['zpsys'] = 'vega'

    # --- Select required columns ---
    out = df[['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']]

    # --- Sort by time (optional but good practice) ---
    out = out.sort_values(by='time')

    # --- Save to sncosmo-readable file ---
    out.to_csv(f"sncosmo_nonia/{SN}_input.dat", sep=' ', index=False)

    print("File saved as sncosmo_input.dat")