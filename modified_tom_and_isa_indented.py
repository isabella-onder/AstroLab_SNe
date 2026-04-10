import pandas as pd
import sys
from astropy.io import fits
import os
import numpy as np
import subprocess
from astropy.coordinates import Angle
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.wcs import WCS

# Load data
mosaic_list = pd.read_csv('mosaic_list.txt', sep=" ", header=None, names=['file'])
radec_var = np.loadtxt('var_sky_position.csv', delimiter=' ')
#file_path = pd.read_csv('filepath_supernovae.csv', header=None)
#file_path = file_path[0].iloc[0]


def checkxy(x, y, sz):
    flag = 0
    if x <= 0 or x > sz[1] or y <= 0 or y >= sz[0] or not np.isfinite(x) or not np.isfinite(y):
        flag = 1
        x = 0
        y = 0
    return x, y, flag


def getphot_revis(img, hdr, radec_var, rad_asec, r_in_asec, r_out_asec):
    spix = 1 / 1.78

    wcs = WCS(hdr)
    coord = SkyCoord(ra=radec_var[0] * u.deg, dec=radec_var[1] * u.deg)
    x_var, y_var = wcs.world_to_pixel(coord)

    rad = rad_asec / spix
    r_in = r_in_asec / spix
    r_out = r_out_asec / spix

    gain = 1.
    dark = 0
    rd = 15.

    sz = img.shape
    x_var, y_var, xy_var_bad = checkxy(x_var, y_var, sz)

    abort = 0
    if xy_var_bad == 1:
        abort = 1
        return 0, abort

    positions = [(x_var, y_var)]
    aperture = CircularAperture(positions, r=rad)
    annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
    annulus_masks = annulus_aperture.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(img)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)

    phot = aperture_photometry(img, aperture)
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture.area
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    phot['aper_sum_bkgsub_err'] = np.sqrt(
        phot['aperture_sum'] + phot['aper_bkg'] + dark + rd**2
    )

    for col in phot.colnames:
        phot[col].info.format = '%.8g'

    return phot, abort


def getphot(img, hdr, sources, rad_asec, r_in_asec, r_out_asec):
    spix = 1 / 1.78

    rad = rad_asec / spix
    r_in = r_in_asec / spix
    r_out = r_out_asec / spix

    gain = 1.
    dark = 0.
    rd = 15.
    linear_lim = 50000.

    wcs = WCS(hdr)
    coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    ra = coords.ra.deg
    dec = coords.dec.deg

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=rad)
    aperture = CircularAperture(positions, r=rad)
    annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

    annulus_masks = annulus_aperture.to_mask(method='center')
    phot_masks = apertures.to_mask(method='center')

    bkg_median = []
    max_val = []
    max_flag = []

    for mask in phot_masks:
        phot_data = mask.multiply(img)
        phot_data_1d = phot_data[mask.data > 0]
        max_tmp = np.nanmax(phot_data_1d)
        max_val.append(max_tmp)
        max_flag.append(1 if max_tmp >= linear_lim else 0)

    max_val = np.array(max_val)
    max_flag = np.array(max_flag)

    for mask in annulus_masks:
        annulus_data = mask.multiply(img)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)

    phot = aperture_photometry(img, aperture)
    phot['max_val'] = max_val
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture.area
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    phot['aper_sum_bkgsub_err'] = np.sqrt(
        phot['aperture_sum'] + phot['aper_bkg'] + dark + rd**2
    )
    phot['aper_SNR'] = phot['aper_sum_bkgsub'] / phot['aper_sum_bkgsub_err']
    phot['max_flag'] = max_flag
    phot['RA'] = ra
    phot['Dec'] = dec

    for col in phot.colnames:
        phot[col].info.format = '%.8g'

    return phot


# Lists
magnitudes = []
mag_err_uppers = []
mag_err_lowers = []
mjds = []
reg_date = []
airmasses = []
no_stars = []
telescopes = []
filters = []
ras = []
decs = []
counts_list = []
counts_errs = []
fwhm_list = []

df = pd.read_excel('2019np.xlsx')
zpt_means = df['zpt_mean'].tolist()
zpt_errs = df['zpt_err'].tolist()

# Main loop
for i in range(len(mosaic_list)):
    mosaic = mosaic_list['file'].iloc[i]
    mosaic_filename = f"{mosaic}"

    with fits.open(mosaic_filename) as hdul:
        hdr = hdul[0].header
        img = hdul[0].data

        date = hdr['DATE-OBS']
        reg_date.append(date)

        t = Time(str(date), format='isot', scale='utc')
        mjds.append(t.mjd)

        airmasses.append(hdr['AIRMASS'])
        telescopes.append(hdr['TELESCOP'])
        filters.append(hdr['FILTER'])
        ras.append(hdr['RA'])
        decs.append(hdr['DEC'])
        fwhm_list.append(hdr['PSF_FWHM'])

        phot, abort = getphot_revis(img, hdr, radec_var, 5, 15, 25)

        mean = zpt_means[i]
        zpt_err = zpt_errs[i]

        if phot == 0:
            magnitudes.append(np.nan)
            mag_err_uppers.append(np.nan)
            mag_err_lowers.append(np.nan)
            no_stars.append(np.nan)
            counts_list.append(np.nan)
            counts_errs.append(np.nan)
        else:
            counts = phot['aper_sum_bkgsub'][0]
            counts_err = phot['aper_sum_bkgsub_err'][0]

            mag = mean - (2.5 * np.log10(counts))

            upper_err = abs(((mean + zpt_err) - (2.5 * np.log10(counts - counts_err))) - mag)
            lower_err = abs(((mean - zpt_err) - (2.5 * np.log10(counts + counts_err))) - mag)

            magnitudes.append(mag)
            mag_err_uppers.append(upper_err)
            mag_err_lowers.append(lower_err)
            counts_list.append(counts)
            counts_errs.append(counts_err)
            no_stars.append(np.nan)

# Plot
plt.figure()
plt.errorbar(mjds, magnitudes, yerr=[mag_err_lowers, mag_err_uppers], fmt='.')
plt.xlabel('MJD date')
plt.ylabel('absolute magnitude')
plt.savefig("mosaic_mag_vs_time.png")
plt.show()