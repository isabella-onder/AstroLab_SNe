from astropy.io import fits
from astropy.wcs import WCS 
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
import numpy as np

from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.psf.matching import create_matching_kernel, TukeyWindow
from photutils.psf import extract_stars
from astropy.nddata import NDData
from astropy.table import Table
from scipy.signal import fftconvolve
import numpy as np


def fwhm_to_psf_size(fwhm, factor=10):
    """Derive PSF cutout size from FWHM. Must be odd."""
    size = int(np.ceil(fwhm * factor))
    if size % 2 == 0:
        size += 1
    return size


def build_empirical_psf(data, size=51, fwhm_guess=4.0, threshold_sigma=10):
    """
    Detect bright isolated stars and median-stack them into an empirical PSF.
    size: cutout size in pixels (must be odd)
    fwhm_guess: used only for star detection, not the PSF shape itself
    """
    data_clean = np.where(np.isnan(data), 0, data)
    _, median, std = sigma_clipped_stats(data_clean, sigma=3.0)

    daofind = DAOStarFinder(fwhm=fwhm_guess, threshold=threshold_sigma * std,
                            exclude_border=True)
    sources = daofind(data_clean - median)

    if sources is None or len(sources) == 0:
        raise ValueError("No stars detected — try lowering threshold_sigma")

    print(f"  Detected {len(sources)} stars")

    margin = size // 2 + 5
    h, w = data.shape
    mask = (
        (sources['xcentroid'] > margin) &
        (sources['xcentroid'] < w - margin) &
        (sources['ycentroid'] > margin) &
        (sources['ycentroid'] < h - margin)
    )
    sources = sources[mask]

    if len(sources) == 0:
        raise ValueError("All detected stars were too close to edges")

    sources.sort('peak')
    sources.reverse()
    sources = sources[:40]

    stars_tbl = Table()
    stars_tbl['x'] = sources['xcentroid']
    stars_tbl['y'] = sources['ycentroid']

    nddata = NDData(data=data_clean)
    stars = extract_stars(nddata, stars_tbl, size=size)

    cutouts = []
    for star in stars:
        arr = star.data.copy().astype(float)
        total = np.sum(arr)
        if total > 0:
            cutouts.append(arr / total)

    if len(cutouts) == 0:
        raise ValueError("No valid star cutouts extracted")

    psf = np.median(np.array(cutouts), axis=0)
    psf /= psf.sum()
    print(f"  PSF built from {len(cutouts)} stars")
    return psf


def match_psfs(before_data, after_aligned, fwhm_before, fwhm_after):
    """
    fwhm_before, fwhm_after: known FWHMs in pixels for each image.
    Convolves the sharper image to match the broader PSF.
    """
    psf_size = max(fwhm_to_psf_size(fwhm_before), fwhm_to_psf_size(fwhm_after))
    print(f"  Using PSF cutout size: {psf_size} px")

    print("Building PSF for BEFORE image...")
    psf_before = build_empirical_psf(before_data, size=psf_size, fwhm_guess=fwhm_before)

    print("Building PSF for AFTER image...")
    psf_after = build_empirical_psf(after_aligned, size=psf_size, fwhm_guess=fwhm_after)

    def psf_width(psf):
        y, x = np.indices(psf.shape)
        cy, cx = np.array(psf.shape) // 2
        r2 = (x - cx)**2 + (y - cy)**2
        return np.sqrt(np.sum(r2 * psf))

    width_before = psf_width(psf_before)
    width_after  = psf_width(psf_after)
    print(f"  PSF width proxy — before: {width_before:.3f} px, after: {width_after:.3f} px")

    window   = TukeyWindow(alpha=0.8)
    nan_mask = np.isnan(before_data) | np.isnan(after_aligned)

    def convolve_safe(image, kernel):
        img = image.copy()
        img[np.isnan(img)] = np.nanmedian(img)
        result = fftconvolve(img, kernel, mode='same')
        result[nan_mask] = np.nan
        return result

    if width_before <= width_after:
        print("  Before image is sharper → convolving BEFORE to match AFTER PSF")
        kernel = create_matching_kernel(psf_before, psf_after, window=window)
        return convolve_safe(before_data, kernel), after_aligned
    else:
        print("  After image is sharper → convolving AFTER to match BEFORE PSF")
        kernel = create_matching_kernel(psf_after, psf_before, window=window)
        return before_data, convolve_safe(after_aligned, kernel)
    
before_filename = "differences/amosaic2019np_19_01_13_R.fits"
before_info = "Before: 2019 Typical R band"
after_filename = "differences/amosaic2019np_20_10_26_R.fits"
after_info = "After: 20_10_26 R band "
fwhm_before_2019 = 5.575
fwhm_after_2019 =  5.03

#after_aligned = after_aligned * 1.5 #scaling by hand after_aligned


#before_filename = "differences/amosaic_2025aebt_25_11_25_R.fits"
#before_info = "Before: 25_11_25 R band"
#after_filename = "differences/amosaic_2025aebt_after.fits" 
#after_info = "After: 26_03_06 R band"
#    after_aligned = after_aligned * 0.7 #scaling by hand after_aligned


#before is 2019-06-12
#after is 2020-10-27
with fits.open(before_filename) as before_hdul, fits.open(after_filename) as after_hdul:

    before_data = before_hdul[0].data #extracting the data from the fits file
    after_data = after_hdul[0].data

    



    before_wcs = WCS(before_hdul[0].header) #applying astropy's world coordinate system
    after_wcs = WCS(after_hdul[0].header)

    after_aligned, footprint = reproject_interp((after_data, after_wcs), before_wcs, shape_out=before_data.shape) #aligning "after" to "before"

    print('a random point on both', before_data[398,260], after_aligned[398,260])


    #clipping them to ensure that they both show the same part of the sky: just putting nan since we are only using this to get values anyways
    mask = footprint > 0
    before_clipped = np.where(mask, before_data, np.nan)
    after_clipped = np.where(mask, after_aligned, np.nan)
    #renaming them back
    before_data = before_clipped
    after_aligned = after_clipped

    #removing background and getting the brightest ones
    flat_before = before_data.flatten()
    flat_before = flat_before[~np.isnan(flat_before)] 
    median_before= np.median(flat_before)

    flat_after = after_aligned.flatten()
    flat_after = flat_after[~np.isnan(flat_after)] 
    median_after= np.median(flat_after)

    after_aligned = after_aligned - median_after
    before_data = before_data - median_before



    print('these are the medians, i.e. the background', median_before, median_after)

   

    #subtract the background on both: rather subtract the median - more rigorously
    #median_after = np.median(after_aligned)  #want to take the median from all rather than the projected data
    #print('median before and after',median_before, median_after)
    #median_before = np.median(before_data)
    #after_aligned = after_aligned - median_after
    #before_data = before_data - median_after

    #getting the brightest ones
    brightest =int(0.005 * len(flat_after))#taking 0.5% of the total number of pixels
    bright_before= np.sort(flat_before)[-brightest:]
    bright_after = np.sort(flat_after[-brightest:])
    

    #taking their median and then scaling accordingly
    before_bright_mean = np.median(bright_before)
    after_bright_mean = np.median(bright_after)
    print('these are the bright before and after means', before_bright_mean, after_bright_mean)
    scaling_factor = before_bright_mean/after_bright_mean #how much I need to multiply after by to get before scale
    print('this is the scaling factor', scaling_factor)

   

    #scaling the after to match the before
    after_aligned = after_aligned * scaling_factor



    #print('a point that looks like a star', before_data[260,398], after_aligned[260,398])
    #print('a point that looks like bulge', before_data[267,570], after_aligned[267,570])
    #after_aligned = after_aligned * 1.5 #scaling by hand after_aligned
    #print('a random point on both', before_data[398,260], after_aligned[398,260])
    #print('a point that looks like a star', before_data[260,398], after_aligned[260,398])
    #print('a point that looks like bulge', before_data[267,570], after_aligned[267,570])

    #difference = before_data - after_aligned #subtracting the two arrays

      # --- PSF MATCHING ---

    print("################Matching PSFs...####################")
    before_psf_matched, after_psf_matched = match_psfs(
        before_data, after_aligned, fwhm_before_2019, fwhm_after_2019
    )


    fig, axs = plt.subplots(nrows =1, ncols = 2)
    #plot before
    norm = ImageNormalize(before_psf_matched, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info+' PSF matched')

    
    #plot after
    norm = ImageNormalize(after_psf_matched, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[1].colorbar()
    fig.colorbar(im, ax=axs[1])
    axs[1].set_xlabel(after_info+' PSF matched')

    plt.show()

    difference = before_psf_matched - after_psf_matched

    

    fig, axs = plt.subplots(nrows =1, ncols = 3)
    #plot before
    norm = ImageNormalize(before_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_data, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info)
    

    
    #plot after
    norm = ImageNormalize(after_aligned, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_aligned, cmap='gray', origin='lower', norm=norm)
    #axs[1].colorbar()
    fig.colorbar(im, ax=axs[1])
    axs[1].set_xlabel(after_info)

    #plotting the difference
    norm = ImageNormalize(difference, interval=ZScaleInterval())
    im = axs[2].imshow(difference, cmap='gray', origin='lower', norm=norm)
    #axs[2].colorbar()
    fig.colorbar(im, ax=axs[2])
    axs[2].set_xlabel('Difference')
    plt.show()

    hdu_diff = fits.PrimaryHDU(difference, header=before_hdul[0].header) #making a fits file out of it
    hdu_diff.writeto('output.fits', overwrite=True)
    print("Difference image saved")