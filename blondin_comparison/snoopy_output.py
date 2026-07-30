import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Read the .dat file, skipping comment lines beginning with '#'
# The file uses commas and spaces as mixed delimiters, so we preprocess

#model_name = 'ebv_model2_dm15'
models = ['ebv_model2_dm15', 'ebv_model2_st', 'max_model_st']
fig, ax = plt.subplots(figsize=(8, 5))
for model in models:
    filter = 'V'
    filename = f'{model}/SN2019np_lc_{filter}s_model.dat'

    # Load the .dat file produced by dump_lc()
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.replace(',', '').split()
            data.append([float(x) for x in parts])

    data    = np.array(data)
    t       = data[:, 0]
    mag     = data[:, 1]
    mag_err = data[:, 2]

    # Cubic spline through the model points (they are noiseless, so no smoothing needed)
    cs = CubicSpline(t, mag)
    t_dense   = np.linspace(t.min(), t.max(), 500)
    mag_dense = cs(t_dense)

   
    ax.plot(t_dense, mag_dense, color='steelblue', lw=1.8, label='Model (interpolated)')
    ax.errorbar(t, mag, yerr=mag_err,
                fmt='o', color='k', ecolor='grey',
                elinewidth=1.2, capsize=3, markersize=4,
                label='Model sample points', zorder=5)
    ax.invert_yaxis()
    ax.set_xlabel('MJD', fontsize=12)
    ax.set_ylabel('Magnitude', fontsize=12)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()