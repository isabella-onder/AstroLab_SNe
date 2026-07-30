import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

data = np.loadtxt('star_mags.dat', comments='#')
print(data.shape)
gaia_mag_v = data[:, 2]
mag_v      = data[:, 5]
err_mag_v  = data[:, 6]

fig, ax = plt.subplots(figsize=(5, 5))

ax.errorbar(gaia_mag_v, mag_v, yerr=err_mag_v, markersize = 14,
            fmt='.', color='dimgrey', capsize=5, label='GAIA stars')

# 1:1 line for reference
#lim = [min(gaia_mag_v.min(), mag_v.min()) - 0.5,
       #max(gaia_mag_v.max(), mag_v.max()) + 0.5]
lim = [(11,12), (16,16)]
xlim = [10.8,16]
ylim = [11.1,16]
ax.plot(xlim, ylim, color='firebrick', linewidth=3, linestyle='--', label = 'best fit')
#ax.set_xlim(lim)
#ax.set_ylim(lim)

ax.set_xlabel('Gaia $V$ magnitude', fontsize=16)
ax.set_ylabel('Instrumental $V$ magnitude', fontsize=16)



ax.minorticks_on()
ax.tick_params(which='both', direction='in', labelsize=16)
ax.tick_params(which='major', length=5, width=2.0)
ax.tick_params(which='minor', length=2, width=1.5)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.legend(frameon=False, fontsize=16)
fig.tight_layout()
plt.savefig('gaia_vs_instrumental.png', dpi=150, bbox_inches='tight')
plt.show()