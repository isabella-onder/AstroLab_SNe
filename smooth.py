from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.nddata import Cutout2D
from astropy.convolution import convolve
from photutils.psf.matching import create_matching_kernel, TukeyWindow

# -----------------------
# 1. Load images
# -----------------------
good_data = fits.getdata("differences/Late_2019.fits")
poor_data = fits.getdata("differences/Typical_2019.fits")

# -----------------------
# 2. Star detection function
# -----------------------
def find_stars(data, fwhm=3.0, threshold_sigma=5.0):
    mean, median, std = sigma_clipped_stats(data)
    finder = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
    sources = finder(data - median)
    return sources

# -----------------------
# 3. Extract star cutouts
# -----------------------
def extract_cutouts(data, sources, size=25, max_stars=30):
    cutouts = []
    for i, star in enumerate(sources):
        if i >= max_stars:
            break
        x, y = star['xcentroid'], star['ycentroid']
        try:
            cutout = Cutout2D(data, (x, y), size)
            img = cutout.data

            # Normalize flux
            img = img / np.sum(img)

            cutouts.append(img)
        except:
            continue
    return np.array(cutouts)

# -----------------------
# 4. Build PSF
# -----------------------
def build_psf(cutouts):
    return np.mean(cutouts, axis=0)

# -----------------------
# 5. Generate PSFs
# -----------------------
sources_good = find_stars(good_data)
sources_poor = find_stars(poor_data)

cutouts_good = extract_cutouts(good_data, sources_good)
cutouts_poor = extract_cutouts(poor_data, sources_poor)

psf_good = build_psf(cutouts_good)
psf_poor = build_psf(cutouts_poor)

# -----------------------
# 6. Create matching kernel
# -----------------------
window = TukeyWindow(alpha=0.5)
kernel = create_matching_kernel(psf_good, psf_poor, window=window)

# -----------------------
# 7. Convolve image
# -----------------------
matched = convolve(good_data, kernel)

# -----------------------
# 8. Save result
# -----------------------
fits.writeto("matched_starstack.fits", matched, overwrite=True)

print("Done: star-stacked PSF matching complete.")