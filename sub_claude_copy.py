from astropy.io import fits
from astropy.wcs import WCS 
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.psf.matching import create_matching_kernel, TukeyWindow
from photutils.psf import extract_stars
from astropy.nddata import NDData
from astropy.table import Table
from scipy.signal import fftconvolve
import numpy as np
from astropy.convolution import Gaussian2DKernel

def fwhm_to_gaussian_psf(fwhm, size=None):
    """Build a clean analytical Gaussian PSF from a FWHM in pixels."""
    sigma = fwhm / 2.355
    if size is None:
        size = int(np.ceil(fwhm * 10))
        if size % 2 == 0:
            size += 1
    psf = Gaussian2DKernel(sigma, x_size=size, y_size=size)
    arr = np.array(psf)
    return arr / arr.sum()


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


def match_psfs(before_data, after_aligned, psf_before, psf_after):
    """
    psf_before, psf_after: 2D numpy arrays of the PSF for each image.
                           Must be the same shape, normalised to sum to 1.
    """
    # Ensure both PSFs are the same size (pad the smaller one if needed)
    def pad_to_match(psf, target_shape):
        py = target_shape[0] - psf.shape[0]
        px = target_shape[1] - psf.shape[1]
        if py < 0 or px < 0:
            raise ValueError("Target shape is smaller than PSF — check your PSF arrays")
        return np.pad(psf, ((py//2, py - py//2), (px//2, px - px//2)))

    target_shape = (max(psf_before.shape[0], psf_after.shape[0]),
                    max(psf_before.shape[1], psf_after.shape[1]))
    # Shapes must be odd for create_matching_kernel
    target_shape = tuple(s + 1 if s % 2 == 0 else s for s in target_shape)

    psf_before = pad_to_match(psf_before, target_shape)
    psf_after  = pad_to_match(psf_after,  target_shape)

    # Renormalise after padding
    psf_before /= psf_before.sum()
    psf_after  /= psf_after.sum()

    def psf_width(psf):
        y, x = np.indices(psf.shape)
        cy, cx = np.array(psf.shape) // 2
        r2 = (x - cx)**2 + (y - cy)**2
        return np.sqrt(np.sum(r2 * psf))

    width_before = psf_width(psf_before)
    width_after  = psf_width(psf_after)
    print(f"  PSF width proxy — before: {width_before:.3f} px, after: {width_after:.3f} px")

    window  = TukeyWindow(alpha=0.4)
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
        
        print(f"  psf_after sum: {psf_after.sum():.6f}")
        print(f"  psf_before sum: {psf_before.sum():.6f}")
        print(f"  kernel sum: {kernel.sum():.6f}")
        print(f"  kernel min/max: {kernel.min():.4f} / {kernel.max():.4f}")
        print(f"  after_aligned nanmedian before convolve: {np.nanmedian(after_aligned):.4f}")

        return convolve_safe(before_data, kernel), after_aligned
    else:
        print("  After image is sharper → convolving AFTER to match BEFORE PSF")
        
        kernel = create_matching_kernel(psf_after, psf_before, window=window)

        print(f"  psf_after sum: {psf_after.sum():.6f}")
        print(f"  psf_before sum: {psf_before.sum():.6f}")
        print(f"  kernel sum: {kernel.sum():.6f}")
        print(f"  kernel min/max: {kernel.min():.4f} / {kernel.max():.4f}")
        print(f"  after_aligned nanmedian before convolve: {np.nanmedian(after_aligned):.4f}")
        return before_data, convolve_safe(after_aligned, kernel)
    
    # Diagnostics
        print(f"  psf_after sum: {psf_after.sum():.6f}")
        print(f"  psf_before sum: {psf_before.sum():.6f}")
        print(f"  kernel sum: {kernel.sum():.6f}")
        print(f"  kernel min/max: {kernel.min():.4f} / {kernel.max():.4f}")
        print(f"  after_aligned nanmedian before convolve: {np.nanmedian(after_aligned):.4f}")
        
    
#################################### 2019np TESTER ######################################################
before_filename = "differences/amosaic2019np_19_01_13_R.fits"
before_info = "Before: 2019 Typical R band"
after_filename = "differences/amosaic2019np_20_10_26_R.fits"
after_info = "After: 20_10_26 R band "
fwhm_before = 5.575
fwhm_after =  5.03

#after_aligned = after_aligned * 1.5 #scaling by hand after_aligned

################################### 2025aebt TESTER ######################################################
#before_filename = "differences/amosaic_2025aebt_25_11_25_R.fits"
#before_info = "Before: 25_11_25 R band"
#after_filename = "differences/amosaic_2025aebt_after.fits" 
#after_info = "After: 26_03_06 R band"
#fwhm_before = 3.83
#fwhm_after = 6.07
#after_aligned = after_aligned * 0.7 #scaling by hand after_aligned





with fits.open(before_filename) as before_hdul, fits.open(after_filename) as after_hdul:

    before_data = before_hdul[0].data #extracting the data from the fits file
    after_data = after_hdul[0].data
    before_wcs = WCS(before_hdul[0].header) #applying astropy's world coordinate system
    after_wcs = WCS(after_hdul[0].header)
    after_aligned, footprint = reproject_interp((after_data, after_wcs), before_wcs, shape_out=before_data.shape) #aligning "after" to "before"

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
    after_aligned_scaled = after_aligned * scaling_factor





    print("################Matching PSFs...####################")
    psf_before = fwhm_to_gaussian_psf(fwhm_before)
    psf_after  = fwhm_to_gaussian_psf(fwhm_after)

    before_psf_matched, after_psf_matched = match_psfs(
    before_data, after_aligned, psf_before, psf_after
    )

    fig, axs = plt.subplots(nrows =1, ncols = 2)
    #plot before
    norm = ImageNormalize(after_aligned, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info+'fresh PSF matched')
    #plot after
    norm = ImageNormalize(after_psf_matched, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[1].colorbar()
    fig.colorbar(im, ax=axs[1])
    axs[1].set_xlabel(after_info+'fresh PSF matched')

    plt.show()

    #appropriately scaling the matched psf ones

    flat_before_psf = before_psf_matched.flatten()
    flat_before_psf = flat_before_psf[~np.isnan(flat_before_psf)] 
    flat_after_psf = after_psf_matched.flatten()
    flat_after_psf = flat_after_psf[~np.isnan(flat_after_psf)] 
    
    #CHECKING THE MEAN OF AFTER_PSF
   
    #getting the brightest ones
    brightest_psf =int(0.005 * len(flat_before_psf))#taking 0.5% of the total number of pixels
    bright_before_psf= np.sort(flat_before_psf)[-brightest_psf:]
    bright_after_psf = np.sort(flat_after_psf[-brightest_psf:])
    
    #taking their median and then scaling accordingly
    before_bright_mean_psf = np.median(bright_before_psf)
    after_bright_mean_psf = np.median(bright_after_psf)
    print('these are the bright before and after means', before_bright_mean_psf, after_bright_mean_psf)
    scaling_factor_psf = np.abs(before_bright_mean_psf/after_bright_mean_psf) #how much I need to multiply after by to get before scale
    print('this is the PSF matched scaling factor', scaling_factor_psf)

    #scaling the after to match the before
    after_psf_matched = after_psf_matched * scaling_factor_psf
    difference_psf = before_psf_matched - after_psf_matched


    

#plotting the psf matched ones and scaled
    fig, axs = plt.subplots(nrows =1, ncols = 3)
    #plot before
    norm = ImageNormalize(before_psf_matched, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info+' PSF matched and scaled')

    
    #plot after
    norm = ImageNormalize(after_psf_matched, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_psf_matched, cmap='gray', origin='lower', norm=norm)
    #axs[1].colorbar()
    fig.colorbar(im, ax=axs[1])
    axs[1].set_xlabel(after_info+' PSF matched and scaled')

    #plotting the difference
    norm = ImageNormalize(difference_psf, interval=ZScaleInterval())
    im = axs[2].imshow(difference_psf, cmap='gray', origin='lower', norm=norm)
    #axs[2].colorbar()
    fig.colorbar(im, ax=axs[2])
    axs[2].set_xlabel('Difference')
    plt.show()

    plt.show()

    difference = before_data - after_aligned_scaled

    

    fig, axs = plt.subplots(nrows =1, ncols = 3)
    #plot before
    norm = ImageNormalize(before_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_data, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info)
    

    
    #plot after
    norm = ImageNormalize(after_aligned_scaled, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_aligned_scaled, cmap='gray', origin='lower', norm=norm)
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

    # Right before PSF matching, after all background/scaling is done
flat_after_final = after_aligned[~np.isnan(after_aligned)].flatten()
brightest = int(0.005 * len(flat_after_final))
after_bright_mean_check = np.median(np.sort(flat_after_final)[-brightest:])
print('after bright mean before PSF matching:', after_bright_mean_check)

# Then after PSF matching
flat_after_psf = after_psf_matched[~np.isnan(after_psf_matched)].flatten()
after_bright_mean_psf = np.median(np.sort(flat_after_psf)[-brightest:])
print('after bright mean after PSF matching:', after_bright_mean_psf)

# These should now be identical since after is returned unchanged