from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.nddata import NDData
from photutils.psf import extract_stars, EPSFBuilder
from astropy.table import Table
from astropy.convolution import convolve
from photutils.psf.matching import create_matching_kernel, TukeyWindow

# -----------------------
# 1. Load images
# -----------------------
good_data = fits.getdata("good_seeing.fits")
poor_data = fits.getdata("poor_seeing.fits")

# -----------------------
# 2. Detect stars
# -----------------------
def detect_sources(data, fwhm=3.0, threshold_sigma=5.0):
    mean, median, std = sigma_clipped_stats(data)
    finder = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
    sources = finder(data - median)
    return sources

# -----------------------
# 3. Convert to table format for EPSFBuilder
# -----------------------
def make_star_table(sources):
    return Table({
        'x': sources['xcentroid'],
        'y': sources['ycentroid']
    })

# -----------------------
# 4. Build ePSF
# -----------------------
def build_epsf(data, sources, size=25):
    nddata = NDData(data)
    stars_tbl = make_star_table(sources)

    stars = extract_stars(nddata, stars_tbl, size=size)

    epsf_builder = EPSFBuilder(oversampling=4, maxiters=10, progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)

    return epsf.data

# -----------------------
# 5. Generate PSFs
# -----------------------
sources_good = detect_sources(good_data)
sources_poor = detect_sources(poor_data)

psf_good = build_epsf(good_data, sources_good)
psf_poor = build_epsf(poor_data, sources_poor)

# -----------------------
# 6. Create matching kernel
# -----------------------
window = TukeyWindow(alpha=0.5)
kernel = create_matching_kernel(psf_good, psf_poor, window=window)

# -----------------------
# 7. Convolve better image
# -----------------------
matched = convolve(good_data, kernel)

# -----------------------
# 8. Save result
# -----------------------
fits.writeto("matched_epsf.fits", matched, overwrite=True)

print("Done: ePSF matching complete.")
