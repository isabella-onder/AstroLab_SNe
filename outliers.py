import csv
import numpy as np
import os
import shutil

SN_list = ['2025acsn', '2019np','2020ue', '2025aebt', '2025afwg', '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc', '2026dix','2025acsn_subbed', '2026kc_subbed' ]


def clean_csv(SN, sigma_threshold=3):
    """
    Removes rows from the photometry CSV where either:
      1. The zpt_mean (col 9) deviates by more than `sigma_threshold` sigma
         from the global mean of zpt_means across all rows.
      2. The (zpt_mean - zpt_median) difference (cols 9, 14) deviates by more
         than `sigma_threshold` sigma from the mean of that difference.

    The original CSV is overwritten in place after a backup is made.

    Parameters
    ----------
    SN : str
        Supernova name, used to locate ../Photometry_final/{SN}.csv
    sigma_threshold : float
        Number of standard deviations beyond which a row is flagged (default 3).
    """
    filepath = f'../Photometry_final/{SN}.csv'
    backup_path = f'../Photometry_final/{SN}_backup.csv'

    # --- Pass 1: read all rows and collect zpt statistics ---
    with open(filepath, mode='r') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    zpt_means   = []
    zpt_medians = []
    valid_indices = []  # rows where both columns are parseable

    for i, row in enumerate(rows):
        try:
            zpt_means.append(float(row[9]))
            zpt_medians.append(float(row[14]))
            valid_indices.append(i)
        except (ValueError, IndexError):
            # Rows that cannot be parsed are left untouched by the sigma cuts
            pass

    zpt_means_arr   = np.array(zpt_means)
    zpt_medians_arr = np.array(zpt_medians)
    differences     = zpt_means_arr - zpt_medians_arr  # signed, matching your diagnostic code

    # --- Compute statistics ---
    zpt_mean_of_means  = np.mean(zpt_means_arr)
    zpt_std_of_means   = np.std(zpt_means_arr)

    diff_mean = np.mean(differences)
    diff_std  = np.std(differences)

    print(f'\n=== Sigma-clipping statistics for SN{SN} ===')
    print(f'  zpt_mean  : mean = {zpt_mean_of_means:.4f},  std = {zpt_std_of_means:.4f}')
    print(f'  (mean-median): mean = {diff_mean:.4f},  std = {diff_std:.4f}')
    print(f'  Rejection bounds (plain zpt)  : '
          f'[{zpt_mean_of_means - sigma_threshold*zpt_std_of_means:.4f}, '
          f'{zpt_mean_of_means + sigma_threshold*zpt_std_of_means:.4f}]')
    print(f'  Rejection bounds (mean-median): '
          f'[{diff_mean - sigma_threshold*diff_std:.4f}, '
          f'{diff_mean + sigma_threshold*diff_std:.4f}]')

    # --- Identify rows to remove ---
    flagged_plain_zpt   = set()
    flagged_diff        = set()

    for j, i in enumerate(valid_indices):
        z = zpt_means_arr[j]
        d = differences[j]

        if abs(z - zpt_mean_of_means) > sigma_threshold * zpt_std_of_means:
            flagged_plain_zpt.add(i)

        if abs(d - diff_mean) > sigma_threshold * diff_std:
            flagged_diff.add(i)

    flagged_all = flagged_plain_zpt | flagged_diff

    print(f'\n  Rows flagged by plain zpt criterion    : {len(flagged_plain_zpt)}')
    print(f'  Rows flagged by mean-median criterion  : {len(flagged_diff)}')
    print(f'  Total unique rows to be removed        : {len(flagged_all)}')

    if flagged_all:
        print('\n  Flagged row details:')
        for i in sorted(flagged_all):
            row = rows[i]
            reasons = []
            if i in flagged_plain_zpt:
                reasons.append(f'plain zpt = {row[9]}')
            if i in flagged_diff:
                reasons.append(f'mean-median = {float(row[9]) - float(row[14]):.4f}')
            print(f'    Row {i+2} (MJD {row[2]}, band {row[4]}): {", ".join(reasons)}')

    # --- Write backup then overwrite original ---
    shutil.copy2(filepath, backup_path)
    print(f'\n  Backup written to: {backup_path}')

    cleaned_rows = [row for i, row in enumerate(rows) if i not in flagged_all]

    with open(filepath, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(cleaned_rows)

    print(f'  Cleaned CSV written to: {filepath}')
    print(f'  Rows removed by plain zpt only        : {len(flagged_plain_zpt - flagged_diff)}')
    print(f'  Rows removed by mean-median only       : {len(flagged_diff - flagged_plain_zpt)}')
    print(f'  Rows before: {len(rows)},  after: {len(cleaned_rows)},  removed: {len(rows) - len(cleaned_rows)}\n')


# --- Run ---
#for SN in SN_list:
clean_csv('2025advo', sigma_threshold=2)