import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

data = np.loadtxt('star_mags.dat', comments='#').reshape(-1, 7)

gaia_mag_v = data[:, 2]
counts = data[:,3]
mag_v      = data[:, 5]
err_mag_v  = data[:, 6]

zeropoints = gaia_mag_v - mag_v
zeropoints = gaia_mag_v - 2.5*np.log10(1/counts)
print(f'Mean zpt:   {np.mean(zeropoints):.3f}')
print(f'Median zpt: {np.median(zeropoints):.3f}')
print(f'Std zpt:    {np.std(zeropoints):.3f}')

fig, ax = plt.subplots(figsize=(5, 4))

ax.hist(zeropoints, bins=12, color='silver', edgecolor='white')
ax.axvline(np.mean(zeropoints),   color='firebrick', linestyle=':', linewidth=3, label=f'Mean = {np.mean(zeropoints):.3f}')
ax.axvline(np.median(zeropoints), color='dodgerblue', linestyle='--',  linewidth=3, label=f'Median = {np.median(zeropoints):.3f}')
ax.axvline(26.5, color='white', linestyle='--',  linewidth=0.1, label=f'Difference = {(np.mean(zeropoints) - np.median(zeropoints)):.3f}')
ax.axvline(26.5, color='white', linestyle='--',  linewidth=0.1, label=f'Zpt error = 0.025')

ax.set_xlabel('Zeropoint (mag)', fontsize=16)
ax.set_ylabel('Star number', fontsize=16)

ax.minorticks_on()
ax.tick_params(which='both', direction='in', labelsize=16)
ax.tick_params(which='major', length=5, width=1.5)
ax.tick_params(which='minor', length=2, width=1.0)
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.legend(frameon=False, fontsize=14)
fig.tight_layout()
plt.savefig('zeropoint_histogram.png', dpi=150, bbox_inches='tight')
plt.show()