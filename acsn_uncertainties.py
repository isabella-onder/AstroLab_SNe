import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D

# ── Data ──────────────────────────────────────────────────────────────────────
data_1a, current = [], None
with open('snpy_acsn/2025acsn.txt') as f:
    for line in f:
        parts = line.split()
        if 'filter' in parts:
            current = [[], [], []]
            data_1a.append(current)
        elif current and len(parts) == 3:
            for i, val in enumerate(parts):
                current[i].append(float(val))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.5, 9))

filters = ['I', 'R', 'V', 'B']
colours = ["#D7A929", 'firebrick', 'mediumseagreen', 'steelblue']
offsets = [-1.3, 0.5, 1.6, 2.5]
shapes  = ['o', 'v', '^', 'D']

Tmax_ue = 0

for band, colour, off, shape in zip(data_1a, colours, offsets, shapes):
    MJD    = [day - Tmax_ue for day in band[0]]
    mags   = [m + off for m in band[1]]
    errors = band[2]
    ax.scatter(MJD, mags, color=colour, s=40, marker=shape)
    ax.errorbar(MJD, mags, yerr=errors, fmt=',', color=colour)

#ax.set_xlim(xmin=-20, xmax=55)
ax.set_ylabel('Uncalibrated apparent magnitude', fontsize=17)
ax.set_xlabel('Time (days)', fontsize=16)
ax.invert_yaxis()

ax.minorticks_on()
ax.tick_params(which='both', direction='in', labelsize=17)
ax.tick_params(which='major', length=5, width=2.0)
ax.tick_params(which='minor', length=2, width=1.3)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    ax.errorbar([], [], yerr=[0.7], fmt=shape, color=colour,
                markersize=8, capsize=2, label=f)
    for f, colour, shape in zip(filters, colours, shapes)
]
ax.legend(handles=legend_elements, frameon=False, fontsize=14)

fig.tight_layout()
plt.savefig("Figures/acsn_uncertainties.png", bbox_inches='tight', dpi=150)
plt.show()
