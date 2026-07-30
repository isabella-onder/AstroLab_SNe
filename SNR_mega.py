import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

# ── Data ──────────────────────────────────────────────────────────────────────

common_axis = list(range(1, 21))

# Individual-axis SNe
ahkt_axis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
ahkt_SNR  = [x + 128 for x in [10.18685335, 16.39607544, 19.77092389, 19.86136633,
                                 19.35809652, 17.16384144, 10.88721495, 15.98270434,
                                 13.1844191,  14.70178074, 11.0173711,  11.0893994]]

ahsa_axis = [15, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
ahsa_SNR  = [x + 120 for x in [9.193185775, 9.517553861, 9.517553861, 9.978014966,
                                 10.55732944,  11.22598883, 12.11408395, 13.22310808,
                                 14.56384092,  16.35590867, 17.38526819, 18.91641164,
                                 19.7129928,   19.56117104, 17.89521626,  9.258671435]]

advo_SNR  = [x + 70 for x in [28.16150779, 48.74924861, 59.93177411, 64.4383758,
                                63.17908452,  59.57037138, 55.62126553, 51.59163229,
                                47.85236073,  44.74188035, 42.25624591, 39.57335381,
                                37.48226155,  35.90907032, 34.74755878, 33.27712981,
                                32.36183393,  31.45988678, 29.76863708, 28.20304604]]

acd_axis = [12, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
acd_SNR  = [83.42448157, 83.42448157, 88.74930679, 93.46571994, 99.24565902,
             105.4648466, 112.1282064, 119.0781104, 124.6266629, 125.6644832,
             118.0818102,  95.36917435]

ue_SNR   = [x + 15 for x in [46.385587,   76.66583083, 91.726019,   96.46793546,
                               96.1720307,  93.18290647, 88.9908769,  84.53577792,
                               79.96629993, 75.31909153, 71.19900228, 67.10756609,
                               63.20150009, 59.25109464, 56.23910065, 53.78934911,
                               51.09867863, 48.57581566, 45.9730353,  44.33983654]]

np_SNR   = [x + 5  for x in [42.90347357, 71.70360969, 93.74122462, 95.7073775,
                               95.09239865, 91.88370891, 87.52727603, 83.47883088,
                               79.14293365, 75.76579231, 71.69342944, 67.84999399,
                               64.43283064, 60.91877211, 57.76818097, 55.48397733,
                               52.76343267, 50.86499234, 49.42547139, 47.75206575]]

# "Bad" SNe
aebt_axis = [2, 3, 4, 5, 6, 7, 8, 9, 10]
aebt_SNR  = [19.19118236, 27.40537131, 33.69884771, 38.41838288, 43.19521957,
              48.024735,   54.40420774, 59.74311661, 65.38301222]

acsn_axis = list(range(2, 21))
acsn_SNR  = [x + 2 for x in [18.19794011, 22.59290124, 24.67082937, 26.69274534,
                               27.98810229, 29.51605876, 30.93479953, 32.93458699,
                               35.44497987, 38.28361653, 41.55336157, 45.45676553,
                               49.26553045, 53.23181231, 57.14085778, 60.84051491,
                               64.24999778, 66.68212744, 68.65464927]]

ahxd_axis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
ahxd_SNR  = [x - 20 for x in [17.29592563, 29.56156577, 37.1845584,  42.83555394,
                                46.32683957, 46.32683957, 43.49543754, 45.89377657,
                                42.06569204, 43.24932667, 48.37131539, 45.90383572]]

afwg_SNR  = [x - 5  for x in [6.455950847, 10.44709652, 11.76130788, 13.8690571,
                                15.97108563, 16.80559211, 17.67,       17.0057889,
                                18.44470188, 19.24972796, 19.74459339, 20.55310361,
                                21.36924636, 21.38724023, 21.35066613, 21.26525463,
                                21.73130031, 22.02224072, 22.42736882, 22.81122972]]

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 7))

# Pixel↔arcsec conversion (plate scale 1.78 pix/arcsec)
def arcsec_to_pix(x): return x * 1.78
def pix_to_arcsec(x): return x / 1.78

secax = ax.secondary_xaxis('top', functions=(arcsec_to_pix, pix_to_arcsec))
secax.set_xlabel("Aperture radius (pix)", size=19)
secax.tick_params(which='both', direction='in', labelsize=16)
secax.tick_params(which='major', length=6, width=1.7)
secax.tick_params(which='minor', length=3, width=1.2)
secax.minorticks_on()

# ── Good SNe (filled circle markers) ─────────────────────────────────────────
good_data = [
    ('2025ahkt', ahkt_axis,   ahkt_SNR,  'darkorchid'),
    ('2025ahsa', ahsa_axis,   ahsa_SNR,  'blue'),
    ('2025advo', common_axis, advo_SNR,  'dodgerblue'),
    ('2026acd',  acd_axis,    acd_SNR,   'chartreuse'),
    ('2020ue',   common_axis, ue_SNR,    'forestgreen'),
    ('2019np',   common_axis, np_SNR,    'silver'),
]

for name, axes, snr, colour in good_data:
    ax.plot([a / 1.78 for a in axes], snr,
            label=name, color=colour,
            marker='o', markersize=5,
            markerfacecolor='none', markeredgecolor=colour)

# ── Blank spacer in legend ────────────────────────────────────────────────────
blank = Line2D([], [], color='none', label='')

# ── Bad SNe (triangle markers) ────────────────────────────────────────────────
bad_data = [
    ('2025aebt', aebt_axis,   aebt_SNR,  'red'),
    ('2025acsn', acsn_axis,   acsn_SNR,  'darkorange'),
    ('2025ahxd', ahxd_axis,   ahxd_SNR,  'burlywood'),
    ('2025afwg', common_axis, afwg_SNR,  'saddlebrown'),
]


for name, axes, snr, colour in bad_data:
    ax.plot([a / 1.78 for a in axes], snr,
            label=name, color=colour,
            marker='^', markersize=5, markeredgecolor=colour)

# ── Legend (explicit handle order with blank spacer) ──────────────────────────
good_handles = ax.lines[:len(good_data)]
bad_handles  = ax.lines[len(good_data):]
ax.legend(handles=[*good_handles, blank, *bad_handles],
          frameon=False, loc=(1.04, 0))

# ── Aperture selection line ───────────────────────────────────────────────────
ax.axvline(x=4 / 1.78, linewidth=3, color='black', alpha=0.8)

# ── Axes formatting ───────────────────────────────────────────────────────────
ax.set_xlabel("Aperture radius ('')", size=19)
ax.set_ylabel('SNR', size=19)
ax.set_ylim(ymin=0)

ax.minorticks_on()
ax.tick_params(which='both',  direction='in', labelsize=15)
ax.tick_params(which='major', length=6,  width=1.5)
ax.tick_params(which='minor', length=3,  width=1.0)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

plt.rc('font',  size=18)
plt.rc('axes',  titlesize=14, labelsize=20)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

fig.tight_layout()
plt.savefig("Figures/SNR.png", dpi=150, bbox_inches='tight')
plt.show()
