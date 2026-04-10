import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
import curves_colour_output as cc




#data is an output reading from the txt file fed into SNooPy
#data = [[band1], [band2], ...] where band1 = [[MJD], [mags], [errors]]



#extracting the data for t types: Ia, Ic and II
data_1a, current = [], None
with open('snpy_txt_clean/2025acsn.txt') as f:
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




# 5-column grid: cols 1,3,5 are wide; cols 2,4 are narrow (with 4 subpanels each)
gs = gridspec.GridSpec(
    4, 5,
    figure=fig,
    width_ratios=[1.8, 1.4, 1.8, 1.4, 1.8],
    hspace=0.05,
    wspace=0.1,
)

# --- Column 0: single tall panel ---
ax0 = fig.add_subplot(gs[:, 0])

# --- Column 1: four stacked panels ---
ax1a = fig.add_subplot(gs[0, 1])
ax1b = fig.add_subplot(gs[1, 1])
ax1c = fig.add_subplot(gs[2, 1])
ax1d = fig.add_subplot(gs[3, 1])

# --- Column 2: single tall panel ---
ax2 = fig.add_subplot(gs[:, 2])

# --- Column 3: four stacked panels ---
ax3a = fig.add_subplot(gs[0, 3])
ax3b = fig.add_subplot(gs[1, 3])
ax3c = fig.add_subplot(gs[2, 3])
ax3d = fig.add_subplot(gs[3, 3])

# --- Column 4: single tall panel ---
ax4 = fig.add_subplot(gs[:, 4])


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
#section of code for plotting the large lightcurve column plots for each type
#################################################################################################
filters = ['I', 'R', 'V', 'B']
colours = ["#D7A929", 'firebrick', 'mediumseagreen', 'steelblue']
offsets = [0,0.4, 0.7, 0.9]
shapes = ['o', 'v', '^', 'D']
Tmax_acsn = 61001.2
Tmax_advo = 61014.8
Tmax_dix = 61089 #since cannot find a reliable source for dix max, just count from when we started to osberve
for band, colour, off,shape in zip(data_1a, colours, offsets, shapes):
    MJD = band[0]
    MJD = [day - Tmax_acsn for day in MJD]
    mags = band[1]
    errors = band[2]
    mags = [m+off for m in mags]
    ax0.set_xlim(xmin = -20, xmax = 100) #such that visually comparatively makes more sense
    ax0.scatter(MJD, mags, color = colour, s = 15, marker = shape)
    ax0.errorbar(MJD, mags,yerr=errors, fmt = ',', color = colour)
    ax0.text(0.05, 0.01,'2025acsn',transform=ax0.transAxes,fontsize=9,va="bottom", ha="left",color="black",fontweight = 'bold')


for band,colour,off,shape in zip(data_1c,colours, offsets,shapes):
    MJD = band[0]
    MJD = [day - Tmax_advo for day in MJD]
    mags = band[1]
    errors = band[2]
    mags = [m+off for m in mags]
    ax2.scatter(MJD, mags, color = colour, s = 15, marker = shape)
    ax2.errorbar(MJD, mags,yerr=errors, fmt = ',', color = colour)
    ax2.text(0.05, 0.01,'2025advo',transform=ax2.transAxes,fontsize=9,va="bottom", ha="left",color="black",fontweight = 'bold')


for band, colour, off,shape in zip(data_2, colours, offsets,shapes):
    MJD = band[0]
    MJD = [day - Tmax_dix for day in MJD]
    mags = band[1]
    errors = band[2]
    mags = [m+off for m in mags]
    ax4.scatter(MJD, mags, color = colour, s = 15, marker = shape)
    ax4.errorbar(MJD, mags,yerr=errors, fmt = ',', color = colour)
    ax4.text(0.05, 0.01,'2026dix',transform=ax4.transAxes,fontsize=9,va="bottom", ha="left",color="black",fontweight = 'bold')


#axes information

ax0.set_ylabel('Uncalibrated apparent magnitude')
ax0.set_xlabel('Time (days)')
ax2.set_xlabel('Time (days)')
ax4.set_xlabel('Time (days)')

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
#section of code for plotting the smaller, colour, plots for type Ia and Ib
#####################################################################################################
curves_ia = cc.get_colour_curves("2025acsn", filepath="snpy_txt_clean/2025acsn.txt")
curves_ic = cc.get_colour_curves("2025advo", filepath="snpy_txt_clean/2025advo.txt")
colour_pairs = ["B-V","V-R","R-I","V-I",]

for colour_diff, ax, pair in zip(curves_ia, col1, colour_pairs):
    diff_name = colour_diff[0]  #i.e. which colours we are subtracting
    MJD = colour_diff[1]
    MJD = [day - Tmax_acsn for day in MJD]
    mag = colour_diff[2]
    error = colour_diff[3]
    ax.set_xlim(xmin = -10, xmax = 70)
    ax.scatter(MJD, mag, marker = '.', s = 15)
    ax.errorbar(MJD, mag, yerr = error, fmt = ',', color = 'dimgrey')
    ax.text(0.05, 0.01,pair,transform=ax.transAxes,fontsize=9,va="bottom", ha="left",color="black",fontweight = 'bold')

for colour_diff, ax, pair in zip(curves_ic, col3, colour_pairs):
    diff_name = colour_diff[0]  #i.e. which colours we are subtracting
    MJD = colour_diff[1]
    MJD = [day - Tmax_advo for day in MJD]
    mag = colour_diff[2]
    error = colour_diff[3]
    ax.scatter(MJD, mag, marker = '.', s = 15)
    ax.errorbar(MJD, mag, yerr = error, fmt = ',', color = 'dimgrey')
    ax.text(0.05, 0.01,pair,transform=ax.transAxes,fontsize=9,va="bottom", ha="left",color="black",fontweight = 'bold')



#plt.savefig("/mnt/user-data/outputs/figure_layout.pdf",bbox_inches="tight", dpi=150)
plt.show()


'''
# ── Example: label each panel so the layout is visible ────────────────────────
for label, ax in [
    ("col 0", ax0), ("col 2", ax2), ("col 4", ax4),
]:
    #ax.text(0.5, 0.5, label, ha="center", va="center",
            transform=ax.transAxes, fontsize=12, color="white")
    #ax.set_facecolor("#1a6080")

for col_label, panel_axes in [("col 1", [ax1a, ax1b, ax1c, ax1d]),
                                ("col 3", [ax3a, ax3b, ax3c, ax3d])]:
    for i, ax in enumerate(panel_axes):
        #ax.text(0.5, 0.5, f"{col_label}\npanel {i+1}",
                ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="white")
        #ax.set_facecolor("#1a6080")
'''