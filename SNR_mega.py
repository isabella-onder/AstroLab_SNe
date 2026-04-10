import matplotlib.pyplot as plt
import seaborn as sns


acsn_axis = data = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
acsn_SNR = data2 = [18.19794011, 22.59290124, 24.67082937, 26.69274534, 27.98810229, 29.51605876, 30.93479953, 32.93458699, 35.44497987, 38.28361653, 41.55336157, 45.45676553, 49.26553045, 53.23181231, 57.14085778, 60.84051491, 64.24999778, 66.68212744, 68.65464927]

ahxd_axis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
ahxd_SNR = [17.29592563, 29.56156577, 37.1845584, 42.83555394, 46.32683957, 46.32683957, 43.49543754, 45.89377657, 42.06569204, 43.24932667, 48.37131539, 45.90383572]

aebt_axis = [2, 3, 4, 5, 6, 7, 8, 9, 10]
aebt_SNR = [19.19118236, 27.40537131, 33.69884771, 38.41838288, 43.19521957, 48.024735, 54.40420774, 59.74311661, 65.38301222]

ahkt_axis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
ahkt_SNR = [10.18685335, 16.39607544, 19.77092389, 19.86136633, 19.35809652, 17.16384144, 10.88721495, 15.98270434, 13.1844191, 14.70178074, 11.0173711, 11.0893994]
ahkt_SNR = [ahkt + 10 for ahkt in ahkt_SNR]

common_axis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
afwg_SNR = [6.455950847, 10.44709652, 11.76130788, 13.8690571, 15.97108563, 16.80559211, 31.10483793, 17.0057889, 18.44470188, 19.24972796, 19.74459339, 20.55310361, 21.36924636, 21.38724023, 21.35066613, 21.26525463, 21.73130031, 22.02224072, 22.42736882, 22.81122972]
advo_SNR = [28.16150779, 48.74924861, 59.93177411, 64.4383758, 63.17908452, 59.57037138, 55.62126553, 51.59163229, 47.85236073, 44.74188035, 42.25624591, 39.57335381, 37.48226155, 35.90907032, 34.74755878, 33.27712981, 32.36183393, 31.45988678, 29.76863708, 28.20304604]
advo_SNR = [advo - 5 for advo in advo_SNR]
np_SNR = [42.90347357, 71.70360969, 93.74122462, 95.7073775, 95.09239865, 91.88370891, 87.52727603, 83.47883088, 79.14293365, 75.76579231, 71.69342944, 67.84999399, 64.43283064, 60.91877211, 57.76818097, 55.48397733, 52.76343267, 50.86499234, 49.42547139, 47.75206575]
np_SNR = [np - 10 for np in np_SNR]
ue_SNR = [46.385587, 76.66583083, 91.726019, 96.46793546, 96.1720307, 93.18290647, 88.9908769, 84.53577792, 79.96629993, 75.31909153, 71.19900228, 67.10756609, 63.20150009, 59.25109464, 56.23910065, 53.78934911, 51.09867863, 48.57581566, 45.9730353, 44.33983654]

#brian unfortunately did not actually put the semi major axes on his - assuming he did one at a time
ahsa_SNR = [9.193185775, 9.517553861, 9.517553861, 9.978014966, 10.55732944, 11.22598883, 12.11408395, 13.22310808, 14.56384092, 16.35590867, 17.38526819, 18.91641164, 19.7129928, 19.56117104, 17.89521626, 9.258671435]
ahsa_axis= [15, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

gm_SNR = [14.98875752, 15.5756862, 14.98077551, 14.98077551, 14.43288356, 13.92958683, 12.96895606, 12.67011875, 11.03651978]
gm_axis = [3, 4, 5, 5, 6, 7, 8, 9, 10]

acd_SNR = [83.42448157, 83.42448157, 88.74930679, 93.46571994, 99.24565902, 105.4648466, 112.1282064, 119.0781104, 124.6266629, 125.6644832, 118.0818102, 95.36917435]
acd_axis = [12, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]


kc_SNR = [14.77745188, 14.77745188, 14.75757733, 14.77472335, 15.24401077, 16.02598146, 16.82022607, 17.67661842, 18.95105014, 20.50934181, 21.99571399, 22.48574526, 20.59094492]
kc_axis= [13, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]

mega_axes = [acsn_axis, ahxd_axis, aebt_axis, ahkt_axis, ahsa_axis, gm_axis, acd_axis, kc_axis ]
mega_SNR = [acsn_SNR, ahxd_SNR, aebt_SNR,ahkt_SNR, ahsa_SNR, gm_SNR, acd_SNR, kc_SNR]
names = ['2025acsn', '2025ahxd', '2025aebt','2025ahkt', '2025ahsa', '2026gm', '2026acd', '2026kc' ]

mega_SNR_2 = [afwg_SNR, advo_SNR, np_SNR, ue_SNR]
names_2 = ['2025afwg', '2025advo', '2019np', '2020ue']


bad_ones = ['2025ahxd', '2025acsn', '2025afwg', '2025aebt']
bad_axes = [ahxd_axis, acsn_axis, common_axis, aebt_axis]
bad_SNR = [ahxd_SNR, acsn_SNR, afwg_SNR, aebt_SNR]

colours = ["#4B7FC0", "#33B65D", "#E2C258", "#D47153"]
colours = sns.color_palette("tab20", n_colors=4)


fig, ax = plt.subplots(ncols=1, nrows = 2, sharex= True, figsize = (5, 8))

#plotting those with individual axes
for axis, snr, name, c in zip(bad_axes, bad_SNR, bad_ones, colours):
    axis = [a/1.78 for a in axis]
    ax[0].scatter(axis, snr, label = name, color = c)
ax[1].set_xlabel("Aperture radius ('')")
ax[0].legend(frameon=False)
ax[0].set_ylabel("SNR")
#plt.legend(frameon=False)
#plt.show()

colours = ["#7C5FCA", "#4B7FC0", "#33B65D", "#E0BA3D", "#DA7D61", "#E45858", "#6B6969"]
colours = ["#6B6969", "#E45858", "#DA7D61", "#E0BA3D" , "#33B65D", "#4B7FC0", "#7C5FCA"]
colours = sns.color_palette("tab20", n_colors=10)
colours[4] = "#E45858"
colours[5] = "#7C5FCA"
colours[6] = "#C161CA"
well_names = ['2026acd', '2020ue', '2019np','2025advo', '2025ahsa', '2026gm', '2025ahkt']
well_behaved = [acd_SNR, ue_SNR, np_SNR, advo_SNR, ahsa_SNR, gm_SNR, ahkt_SNR]
well_behaved_axes = [acd_axis, common_axis, common_axis, common_axis, ahsa_axis, gm_axis, ahkt_axis]
for axes, good, name,c in zip(well_behaved_axes, well_behaved, well_names, colours):
    axes = [a/1.78 for a in axes]
    ax[1].scatter(axes, good, label = name, color = c)
plt.axvline(x=4/1.78, linewidth = 3, color = '#13A413', alpha = 0.8) #this is the mean max across all good SNe
ax[1].legend(frameon = False)
plt.ylabel('SNR')
#plt.xlim(xmin = 0 )
fig.tight_layout()

plt.rc('font',   size=14)
plt.rc('axes',   titlesize=14, labelsize=16)
plt.rc('xtick',  labelsize=11)
plt.rc('ytick',  labelsize=11)


from matplotlib.ticker import AutoMinorLocator
for ax in ax:
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=6, width=1.5)
    ax.tick_params(which='minor', length=3, width=1.0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.savefig("Figures/SNR.png", dpi = 150, bbox_inches = 'tight')
plt.show()
