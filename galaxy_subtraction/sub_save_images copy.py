from astropy.io import fits
from astropy.wcs import WCS 
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
import numpy as np


before_filename = "amosaic2019np_19_01_13_R.fits"
before_info = "Before: 2019 Typical R band"
after_filename = "amosaic2019np_20_10_26_R.fits"
after_info = "After: 20_10_26 R band"

with fits.open(before_filename) as before_hdul, fits.open(after_filename) as after_hdul:

    before_data = before_hdul[0].data
    after_data = after_hdul[0].data

    # --- Save raw images individually ---
    for data, label, fname in [
        (before_data, before_info + ' RAW', 'raw_before.png'),
        (after_data,  after_info  + ' RAW', 'raw_after.png'),
    ]:
        fig, ax = plt.subplots()
        norm = ImageNormalize(data, interval=ZScaleInterval())
        im = ax.imshow(data, cmap='gray', origin='lower', norm=norm)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel(label)
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close(fig)

    before_wcs = WCS(before_hdul[0].header)
    after_wcs = WCS(after_hdul[0].header)

    after_aligned, footprint = reproject_interp((after_data, after_wcs), before_wcs, shape_out=before_data.shape)

    print('a random point on both', before_data[398,260], after_aligned[398,260])

    mask = footprint > 0
    before_data = np.where(mask, before_data, np.nan)
    after_aligned = np.where(mask, after_aligned, np.nan)

    flat_before = before_data.flatten()
    flat_before = flat_before[~np.isnan(flat_before)]
    median_before = np.median(flat_before)

    flat_after = after_aligned.flatten()
    flat_after = flat_after[~np.isnan(flat_after)]
    median_after = np.median(flat_after)

    after_aligned = after_aligned - median_after
    before_data = before_data - median_before

    print('these are the medians, i.e. the background', median_before, median_after)

    brightest = int(0.005 * len(flat_after))
    bright_before = np.sort(flat_before)[-brightest:]
    bright_after = np.sort(flat_after[-brightest:])

    before_bright_mean = np.median(bright_before)
    after_bright_mean = np.median(bright_after)
    print('these are the bright before and after means', before_bright_mean, after_bright_mean)

    scaling_factor = before_bright_mean / after_bright_mean
    print('this is the scaling factor', scaling_factor)

    after_aligned = after_aligned * scaling_factor

    difference = before_data - after_aligned

    # --- Save aligned/difference images individually ---
    for data, label, fname in [
        (before_data,   before_info,  'aligned_before.png'),
        (after_aligned, after_info,   'aligned_after.png'),
        (difference,    'Difference', 'difference.png'),
    ]:
        fig, ax = plt.subplots()
        norm = ImageNormalize(data, interval=ZScaleInterval())
        im = ax.imshow(data, cmap='gray', origin='lower', norm=norm)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel(label)
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close(fig)

    hdu_diff = fits.PrimaryHDU(difference, header=before_hdul[0].header)
    hdu_diff.writeto('output.fits', overwrite=True)
    print("Difference image saved")