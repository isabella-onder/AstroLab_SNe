import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
import curves_colour_output as cc


# data is an output reading from the txt file fed into SNooPy
# data = [[band1], [band2], ...] where band1 = [[MJD], [mags], [errors]]


# Extracting the data for SN types: Ia, Ic, and II
data_1a, current = [], None
with open('snpy_txt_clean/2020ue.txt') as f:
    for line in f:
        parts = line.split()
        if 'filter' in parts:
            current = [[], [], []]
            data_1a.append(current)
        elif current and len(parts) == 3:
            for i, val in enumerate(parts):
                current[i].append(float(val))


data_1c, current = [], None
with open('snpy_txt_clean/2025advo.txt') as f:
    for line in f:
        parts = line.split()
        if 'filter' in parts:
            current = [[], [], []]
            data_1c.append(current)
        elif current and len(parts) == 3:
            for i, val in enumerate(parts):
                current[i].append(float(val))

data_2, current = [], None
with open('snpy_txt_clean/2026dix.txt') as f:
    for line in f:
        parts = line.split()
        if 'filter' in parts:
            current = [[], [], []]
            data_2.append(current)
        elif current and len(parts) == 3:
            for i, val in enumerate(parts):
                current[i].append(float(val))

fig = plt.figure(figsize=(14, 6))


# 7-column grid:
#   col 0 — wide light curve (Ia)
#   col 1 — spacer
#   col 2 — narrow colour curves (Ia)
#   col 3 — wide light curve (Ic)
#   col 4 — spacer
#   col 5 — narrow colour curves (Ic)
#   col 6 — wide light curve (II)
#
# Spacer columns (1 and 4) have width 0.7, which is half the colour-curve
# column width of 1.4, so the gap between the light curve and colour curve
# panels is effectively doubled relative to the original layout.
gs = gridspec.GridSpec(
    4, 7,
    figure=fig,
    width_ratios=[1.8,  1.4, 0.7, 1.8,  1.4,0.7, 1.8],
    hspace=0.05,
    wspace=0.1,
)

# --- Column 0: single tall panel (Ia light curve) ---
ax0 = fig.add_subplot(gs[:, 0])

# --- Column 2: four stacked panels (Ia colour curves) ---
ax1a = fig.add_subplot(gs[0, 1])
ax1b = fig.add_subplot(gs[1, 1])
ax1c = fig.add_subplot(gs[2, 1])
ax1d = fig.add_subplot(gs[3, 1])

# --- Column 3: single tall panel (Ic light curve) ---
ax2 = fig.add_subplot(gs[:, 3])

# --- Column 5: four stacked panels (Ic colour curves) ---
ax3a = fig.add_subplot(gs[0, 4])
ax3b = fig.add_subplot(gs[1, 4])
ax3c = fig.add_subplot(gs[2, 4])
ax3d = fig.add_subplot(gs[3, 4])

# --- Column 6: single tall panel (II light curve) ---
ax4 = fig.add_subplot(gs[:, 6])


# Collect axes for convenience
axs = {
    "col0":  ax0,
    "col1": [ax1a, ax1b, ax1c, ax1d],
    "col2":  ax2,
    "col3": [ax3a, ax3b, ax3c, ax3d],
    "col4":  ax4,
}

col1 = [ax1a, ax1b, ax1c, ax1d]
col3 = [ax3a, ax3b, ax3c, ax3d]


#################################################################################################
# Light curve column plots for each SN type
#################################################################################################
filters  = ['I', 'R', 'V', 'B']
colours  = ["#D7A929", 'firebrick', 'mediumseagreen', 'steelblue']
offsets  = [-0.2, 0.5, 0.8, 1.0]
shapes   = ['o', 'v', '^', 'D']

#Tmax_acsn = 61001.2
Tmax_advo = 61014.8
Tmax_dix  = 61089  # no reliable published maximum for 2026dix; epoch relative to first observation
Tmax_acsn = 0
Tmax_ue = 58872.7



for band, colour, off, shape in zip(data_1a, colours, offsets, shapes):
    MJD    = [day - Tmax_ue for day in band[0]]
    mags   = [m + off for m in band[1]]
    errors = band[2]
    ax0.set_xlim(xmin=-20, xmax=55)
    ax0.scatter(MJD, mags, color=colour, s=15, marker=shape)
    ax0.errorbar(MJD, mags, yerr=errors, fmt=',', color=colour)
    ax0.text(0.05, 0.01, '2025acsn', transform=ax0.transAxes,
             fontsize=9, va="bottom", ha="left", color="black", fontweight='bold')

for band, colour, off, shape in zip(data_1c, colours, offsets, shapes):
    MJD    = [day - Tmax_advo for day in band[0]]
    mags   = [m + off for m in band[1]]
    errors = band[2]
    ax2.scatter(MJD, mags, color=colour, s=15, marker=shape)
    ax2.errorbar(MJD, mags, yerr=errors, fmt=',', color=colour)
    ax2.text(0.05, 0.01, '2025advo', transform=ax2.transAxes,
             fontsize=9, va="bottom", ha="left", color="black", fontweight='bold')

for band, colour, off, shape in zip(data_2, colours, offsets, shapes):
    MJD    = [day - Tmax_dix for day in band[0]]
    mags   = [m + off for m in band[1]]
    errors = band[2]
    ax4.scatter(MJD, mags, color=colour, s=15, marker=shape)
    ax4.errorbar(MJD, mags, yerr=errors, fmt=',', color=colour)
    ax4.text(0.05, 0.01, '2026dix', transform=ax4.transAxes,
             fontsize=9, va="bottom", ha="left", color="black", fontweight='bold')


# Axis labels and formatting
ax0.set_ylabel('Uncalibrated apparent magnitude')
ax0.set_xlabel('Time (days)')
ax2.set_xlabel('Time (days)')
ax4.set_xlabel('Time (days)')

ax0.set_title("Type Ia", fontsize=11, fontweight="bold")
ax2.set_title("Type Ic", fontsize=11, fontweight="bold")
ax4.set_title("Type II", fontsize=11, fontweight="bold")

ax0.invert_yaxis()
ax2.invert_yaxis()
ax4.invert_yaxis()

for ax in col1:
    ax.invert_yaxis()
    ax.yaxis.tick_right()
for ax in col3:
    ax.invert_yaxis()
    ax.yaxis.tick_right()


#####################################################################################################
# Colour curve panels for SN Ia (2025acsn) and SN Ic (2025advo)
#####################################################################################################
curves_ia = cc.get_colour_curves("2020ue", filepath="snpy_txt_clean/2020ue.txt")
curves_ic = cc.get_colour_curves("2025advo", filepath="snpy_txt_clean/2025advo.txt")
colour_pairs = ["B-V", "V-R", "R-I", "V-I"]

for colour_diff, ax, pair in zip(curves_ia, col1, colour_pairs):
    MJD   = [day - Tmax_ue for day in colour_diff[1]]
    mag   = colour_diff[2]
    error = colour_diff[3]
    ax.set_xlim(xmin=-10, xmax=55)
    ax.scatter(MJD, mag, marker='.', s=15)
    ax.errorbar(MJD, mag, yerr=error, fmt=',', color='dimgrey')
    ax.text(0.05, 0.01, pair, transform=ax.transAxes,
            fontsize=9, va="bottom", ha="left", color="black", fontweight='bold')

for colour_diff, ax, pair in zip(curves_ic, col3, colour_pairs):
    MJD   = [day - Tmax_advo for day in colour_diff[1]]
    mag   = colour_diff[2]
    error = colour_diff[3]
    ax.scatter(MJD, mag, marker='.', s=15)
    ax.errorbar(MJD, mag, yerr=error, fmt=',', color='dimgrey')
    ax.text(0.05, 0.01, pair, transform=ax.transAxes,
            fontsize=9, va="bottom", ha="left", color="black", fontweight='bold')


# plt.savefig("/mnt/user-data/outputs/figure_layout.pdf", bbox_inches="tight", dpi=150)
plt.show()
