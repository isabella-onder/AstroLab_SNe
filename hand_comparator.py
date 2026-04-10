import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../Photometry_final/2025advo.csv')

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

bands = ['B', 'V']
auto_colours = {'B': 'steelblue', 'V': 'mediumseagreen'}
labels = {'B': 'B band', 'V': 'V band'}

for ax, band in zip(axes, bands):
    subset = df[df['Filter'] == band]

    # Automatic measurements
    auto = subset.dropna(subset=['magnitude', 'mag_err_upper'])
    ax.errorbar(
        auto['MJD'], auto['magnitude'],
        yerr=auto['mag_err_upper'],
        fmt='o', color=auto_colours[band],
        markersize=4, capsize=2, linewidth=0.8,
        label='Automatic'
    )

    # Hand measurements (only rows where hand mags exist)
    hand = subset.dropna(subset=['magnitude_hand', 'mag_err_hand'])
    ax.errorbar(
        hand['MJD'], hand['magnitude_hand'],
        yerr=hand['mag_err_hand'],
        fmt='s', color='grey',
        markersize=4, capsize=2, linewidth=0.8,
        label='By hand'
    )

    ax.set_title(labels[band])
    ax.set_xlabel('MJD')
    ax.invert_yaxis()  # Magnitudes increase downward
    ax.legend()
axes[0].invert_yaxis()
axes[0].set_ylabel('Magnitude')

plt.suptitle('Photometric Light Curves', y=1.02)
plt.tight_layout()
plt.savefig('lightcurves_BV.pdf', bbox_inches='tight')
plt.show()