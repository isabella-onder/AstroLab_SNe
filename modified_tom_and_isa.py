#/home/astrolab/Desktop/local_database_files/query_output_files/query_results.csv  => output of the terminal query

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
from astropy.coordinates import SkyCoord
import astropy.units as u

#want to import info from the collection of images instead

#somehow put the collection of images/files into a list

mosaic_list = pd.read_csv('mosaic_list.txt', sep=" ", header=None, names = ['file'])
#mosaic_list = mosaic_list[mosaic_list['file'] == '26_01_30/pt5m/I/amosaic.fits']
#load in the RA and Dec Positions
radec_var = np.loadtxt('var_sky_position.csv', delimiter=' ')

file_path = pd.read_csv('filepath_supernovae.csv', header=None)
file_path = file_path[0].iloc[0]


#not sure if this fct is still needed
    

#revised checkxy (no printlog needed due to issues)
#revised checkxy (no printlog needed due to issues)
def checkxy(x,y,sz):
    flag = 0
    if x <= 0 or x > sz[1] or y <= 0 or y >= sz[0] or np.isfinite(x)==False or np.isfinite(y)==False:
        flag = 1
        x = 0
        y = 0
        #rintlog('WARNING: RA/Dec are not in image range')
    return x,y,flag

#slightly revised get_phot function (gives the specific values needed)
def getphot_revis(img,hdr,radec_var,rad_asec,r_in_asec,r_out_asec):
    # get pixel scale of image
    #spix = getscale(hdr) #imported #we are using only pt5m
    #printlog('pixel scale is: '+str(spix))
    # (hdr)
    #print(radec_var)
    spix = 1 / 1.78  #we are using only pt5m

    # convert RA/Dec of variable star to x/y coordinates on this image
    #x_var, y_var = ad2xy(hdr, radec_var[0], radec_var[1])  # imported
    #printlog('x/y of var, cal1, cal2: '+str(x_var)+' '+str(y_var)+ '')
    wcs = WCS(hdr)
    coord = SkyCoord(ra=radec_var[0]*u.deg, dec=radec_var[1]*u.deg)
    x_var, y_var = wcs.world_to_pixel(coord)    

    # convert aperture sizes from arcseconds to pixels
    rad = rad_asec / spix        # pixels
    r_in = r_in_asec / spix      # pixels
    r_out = r_out_asec / spix    # pixels

    # CCD parameters (used for noise calculation)
    gain = 1.
    dark = 0
    rd = 15.

    # check if x/y coordinates are on the frame.
    sz = img.shape
    x_var, y_var, xy_var_bad = checkxy(x_var, y_var, sz)  # imported
    print(x_var, y_var, xy_var_bad)
    # if any of the x/y positions of the variable are not in the imagegetphot_revis(img,hdr,radecvar,rad_asec,r_in_asec,r_out_asec,file), no photometry is calculated.
    abort = 0
    if xy_var_bad == 1:
        #word = file+' ' + str(x_var) +' '+ str(y_var)+ ''
        #g.write(word+'\n')
        abort = 1
        return 0, abort

    if abort == 0:
        # create tuple of x/y positions
        positions = [(x_var, y_var)]
        # create circular aperture (for photometry calculation)
        aperture = CircularAperture(positions, r=rad)  # imported in photutils
        # create annulus aperture (for sky calculation)
        annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        # create mask
        annulus_masks = annulus_aperture.to_mask(method='center')

        # calculate the data and (sigma clipped) sky
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(img)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        # write results to output object
        phot = aperture_photometry(img, aperture)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * aperture.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        # noise calculation
        phot['aper_sum_bkgsub_err'] = np.sqrt((phot['aperture_sum']) + (phot['aper_bkg']) + dark + rd**2)
        for col in phot.colnames:
            phot[col].info.format = '%.8g'  # for consistent table output
        #printlog(phot)
    return(phot,abort)

#regular getphot function (no printlog fct)
def getphot(img,hdr,sources,rad_asec,r_in_asec,r_out_asec):

    # get pixel scale of image
    #spix = getscale(hdr)
    spix = 1 /1.78
    #printlog('pixel scale is: '+str(spix))

    # convert aperture sizes from arcseconds to pixels
    rad = rad_asec / spix        # pixels
    r_in = r_in_asec / spix      # pixels
    r_out = r_out_asec / spix    # pixels
    
    # CCD parameters (used for noise calculation)
    gain = 1.
    dark = 0.
    rd = 15.
    linear_lim = 50000. # linearity limit.  ignore stars if any pixel within aperture exceeds linear_lim

    #ra,dec = xy2ad(hdr,sources['xcentroid'], sources['ycentroid'])
    wcs = WCS(hdr)
    coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    ra = coords.ra.deg
    dec = coords.dec.deg
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=rad)
    ## create circular aperture (for photometry calculation)
    aperture = CircularAperture(positions, r=rad)
    ## create annulus aperture (for sky calculation)
    annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
    ## create mask
    annulus_masks = annulus_aperture.to_mask(method='center')
    phot_masks = apertures.to_mask(method='center')
    ## calcullate the data and (sigma clipped) sky
    bkg_median = []

    # loop over apertures and find maximum value of pixels within aperture (used to reject no linearity)
    max_val = []
    max_flag = []    
    for mask in phot_masks:
        phot_data = mask.multiply(img)
        phot_data_1d = phot_data[mask.data > 0]
        max_tmp = np.nanmax(phot_data_1d)
        max_val.append(max_tmp)
        if max_tmp >= linear_lim:
            max_flag.append(1)
        if max_tmp < linear_lim:
            max_flag.append(0)        
    max_val = np.array(max_val)
    max_flag = np.array(max_flag)

    # loop over sky annuli and get sig-clipped sky value
    for mask in annulus_masks:
        annulus_data = mask.multiply(img)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)

    # write results to output object
    phot = aperture_photometry(img, aperture)
    phot['max_val'] = max_val
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture.area
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    # noise calculation
    phot['aper_sum_bkgsub_err'] = np.sqrt((phot['aperture_sum']) + (phot['aper_bkg']) + dark + rd**2)
    phot['aper_SNR'] =  phot['aper_sum_bkgsub'] / phot['aper_sum_bkgsub_err']
    #print(len(max_val),len(max_flag))
    phot['max_flag'] =  max_flag

    phot['RA'] = ra
    phot['Dec'] = dec

    for col in phot.colnames:
        phot[col].info.format = '%.8g'  # for consistent table output

    return(phot)


#refcat, no printlog
def get_refcat(hdr,refcat_rad,maglim):
    naxis1 = float(hdr['NAXIS1'])
    naxis2 = float(hdr['NAXIS2'])
    xc = naxis1/2.
    yc = naxis2/2.

    rac,decc = xy2ad(hdr,xc,yc)
    #printlog('RA/Dec center of frame: '+str(rac)+' '+str(decc))
    command = '/mnt/64bin/auto_astrom/refcat '+str(rac)+' '+str(decc)+' -rad '+str(refcat_rad)+' -all > refcat_output'
    #print(command)
    #sys.exit()
    refcat_output = subprocess.getoutput(command)
    #refcat_output = subprocess.getoutput('cat mm > refcat_output')
    data = np.genfromtxt('refcat_output')
    ra = np.array(data[:,0])
    dec = np.array(data[:,1])
    #pm_ra = np.array(data[:,4])
    #pm_dec = np.array(data[:,5])
    #ra,dec = update_for_propermotions(hdr,ra,dec,pm_ra,pm_dec)
    
    gaia_gmag = np.array(data[:,21])
    gaia_gmag_err = np.array(data[:,22])

    gaia_rmag = np.array(data[:,25])
    gaia_rmag_err = np.array(data[:,26])

    gaia_imag = np.array(data[:,29])
    gaia_imag_err = np.array(data[:,30])

    ps_b    =  0.213 + 0.587 * (gaia_gmag - gaia_rmag) + gaia_gmag
    ps_v    =  0.006 + 0.474 * (gaia_gmag - gaia_rmag) + gaia_rmag
    ps_r_c  = -0.138 - 0.131 * (gaia_gmag - gaia_rmag) + gaia_rmag
    ps_i_c  = -0.367 - 0.149 * (gaia_gmag - gaia_rmag) + gaia_imag

    bad = np.where(gaia_gmag >=maglim)
    ra[bad] = 0
    dec[bad] = -90
    
    refcat = Table()
    refcat['RA'] = ra
    refcat['Dec'] = dec
    refcat['bmag'] = ps_b
    refcat['bmag_err'] = gaia_gmag_err
    refcat['vmag'] = ps_v
    refcat['vmag_err'] = gaia_gmag_err
    refcat['gmag'] = gaia_gmag
    refcat['gmag_err'] = gaia_gmag_err
    refcat['rmag'] = ps_r_c
    refcat['rmag_err'] = gaia_rmag_err
    refcat['imag'] = ps_i_c
    refcat['imag_err'] = gaia_imag_err
    for col in refcat.colnames:
        refcat[col].info.format = '%.8g'  # for consistent table output

    #print(refcat)
    return(refcat)

#get_automag_pars revised with no printlog
def get_automag_pars():

	#file = '/mnt/64bin/auto_astrom/automag_driver'

	#if not os.path.exists('automag_driver'):
		#printlog('No automag_driver in current directory.  Using default at /mnt/64bin/auto_astrom/automag_driver')
	#if not os.path.exists(file):
	#	#printlog('No default autmag_driver file.  Will not continue.  Aborting...')
	#	sys.exit()
     
    #initialising
    file = ''

	if os.path.exists('automag_driver'):
		file = 'automag_driver'
    else:
        print('There is no automag driver file that was found')

	#printlog('automag_driver: using file '+file)
	data = open(file,'r')
	n = len(open(file).readlines(  ))
	for i in range(0,n):
		line = data.readline()
		word = line.split()
		if word[0] == 'radius':
			rad = float(word[2])
		if word[0] == 'radius_1':
			r_in = float(word[2])
		if word[0] == 'radius_2':
			r_out = float(word[2])
	#printlog('aperture parameters; rad, r_in, r_out: '+str(rad)+', '+str(r_in)+', '+str(r_out))
	return rad,r_in,r_out

#lists for use in the final table
zpt_means = []
zpt_errs = []
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
zpt_medians = []
counts_list = []
counts_errs = []
fwhm_list = []

import pandas as pd
df = pd.read_excel('2019np.xlsx')
zpt_means = df['zpt_mean'].tolist()
zpt_errors = df['zpt_err'].tolist()

#calculations for each row of the table
for i in range(0, len(mosaic_list)):
	
	mosaic = mosaic_list['file'].iloc[i]
	
	mosaic_filename = f"/mnt/share/supernovae/{file_path}/{mosaic}"
	#(f"magnitude of {name_of_obj}")
    # Open that specific FITS file
	with fits.open(mosaic_filename) as hdul:
		hdr = hdul[0].header  
		img = hdul[0].data
        
        #retriving the date and finding MJD
        
		date = hdr['DATE-OBS']
		reg_date.append(date)
		#import practice file
		#output_list.txt
		#'2025-09-26T00:33:02.840'
		#date_rep = date.replace("T", " ")
		#print(date_rep)
		t= Time(str(date), format='isot', scale='utc')
		mjd_time = t.mjd
		
		mjds.append(mjd_time)
        
        #retriving airmass
		airmass = hdr['AIRMASS']
		airmasses.append(airmass)
		
		#retriving telescope
		
		telescope = hdr['TELESCOP']
		telescopes.append(telescope)
		
		filterx = hdr['FILTER']
		filters.append(filterx)
		
		ra = hdr['RA']
		ras.append(ra)
		
		dec = hdr['DEC']
		decs.append(dec) #(f"magnitude of sn2026GM")
		
		fwhm = hdr['PSF_FWHM']
		fwhm_list.append(fwhm)
		
		#zeropoint calculation -- I CANNOT USE THIS LINE BUT I DO NOT NEED TO
		#mean,median,mean_weighted,std,zpt_err,zpt_gd,mag_phot_gd,counts_gd,mag_phot,counts, aper_SNR, nstars =zeropoint_revis(img,hdr, radec_var)
		#zpt_means.append(mean)
		#zpt_medians.append(median)
		#zpt_errs.append(zpt_err)
        
        phot, abort = getphot_revis(img,hdr,radec_var,5,15,25)
     
        mean = zpt_means[i]
        zpt_err = zpt_errs[i]
        
        if phot == 0:
            #accounting for results where photometry can't be calculated
            magnitudes.append(np.nan)
            mag_err_uppers.append(np.nan)
            mag_err_lowers.append(np.nan)
            no_stars.append(np.nan)
            counts_list.append(np.nan)
            counts_errs.append(np.nan)

       
     
        else: 
            #calculating magnitudes
            counts = phot['aper_sum_bkgsub'][0]
            counts_err = phot['aper_sum_bkgsub_err'][0]
            
            no_stars.append(nstars)
            print(counts)
                            
            mag = mean - (2.5 * np.log10(counts))
            
            upper_err = np.abs(((mean + zpt_err) - (2.5 * np.log10(counts - counts_err))) - mag)
            lower_err = np.abs(((mean - zpt_err) - (2.5 * np.log10(counts + counts_err))) - mag)
            
            
            
            magnitudes.append(mag)
            mag_err_uppers.append(upper_err)
            mag_err_lowers.append(lower_err)
            counts_list.append(counts)
            counts_errs.append(counts_err)
			
plate_scales = []
for telescope in telescopes:
    if telescope == 'pt5m':
        plate_scale = 1/1.78
        plate_scales.append(plate_scale)
    elif telescope =='14-inch West Dome':
        plate_scale = 0.91
        plate_scales.append(plate_scale)
    elif telescope == '16-inch Far-East Dome': #seems like online it is 20 but it in 2019np it calls it like this
        plate_scale =1.07
        plate_scales.append(plate_scale)
    elif telescope == '14-inch DRACO2 Dome':
        plate_scale = 0.983
        plate_scales.append(plate_scale)
    elif telescope == '16-inch East Dome':
        plate_scale = 0.93
        plate_scales.append(plate_scale)
    else:
        print('this telescope is not included!')
fwhm_arcsecs = [plate_scale * fwhm for plate_scale, fwhm in zip(plate_scales,fwhm_list)]

mosaic_list['zpt_mean'] = zpt_means
mosaic_list['zpt_err'] = zpt_errs
mosaic_list['magnitude'] = magnitudes
mosaic_list['mag_err_upper'] = mag_err_uppers
mosaic_list['mag_err_lower'] = mag_err_lowers
mosaic_list['MJD'] = mjds
mosaic_list['airmass'] = airmasses
mosaic_list['#_stars'] = no_stars
mosaic_list['date'] = reg_date
mosaic_list['Telescope'] = telescopes
mosaic_list['plate_scale'] = plate_scales
mosaic_list['Filter'] = filters
mosaic_list['RA/deg'] = ras
mosaic_list['Dec/deg'] = decs
mosaic_list['zpt_median'] = zpt_medians
mosaic_list['counts_phot'] = counts_list
mosaic_list['counts_phot_err'] = counts_errs
mosaic_list['psf_fwhm_pix'] = fwhm_list
mosaic_list['psf_fwhm_arcsecs'] = fwhm_arcsecs


#delete any observations where the zpoint could not be found
df_1= mosaic_list.drop(mosaic_list.index[(mosaic_list["zpt_mean"] == 'fail')],axis=0)
         
#final table to be printed
df_final = df_1[['file','Telescope','plate_scale', 'MJD', 'date', 'Filter','RA/deg', 'Dec/deg', '#_stars', 'airmass', 'zpt_mean', 'zpt_err', 'magnitude', 'mag_err_upper', 'mag_err_lower', 'zpt_median', 'counts_phot', 'counts_phot_err', 'psf_fwhm_pix', 'psf_fwhm_arcsecs']]

#name_find = df_1.iloc[0]picture
#name_of_obj = name_find['Object']

mean_mags = [np.sum(np.array(df_final['magnitude']))/len(df_final['MJD'])] * len(df_final['MJD'])

#plotting and saving simple fig of mag vs time
plt.figure()
plt.errorbar(df_final['MJD'], df_final['magnitude'], yerr = [df_final['mag_err_lower'], df_final['mag_err_upper']], fmt = '.')
plt.xlabel('MJD date')
plt.ylabel('absolute magnitude')
#plt.title(f"magnitude of {name_of_obj}")
plt.savefig(f"mosaic_mag_vs_time.png")
plt.show()



df_final.to_csv(f"mosaic_photometry_darkmagic.csv", index=False)


#contents of a hdr
#SIMPLE  =                    T / conforms to FITS standard              
#        BITPIX  =                  -32 / array data type                 
#                       NAXIS   =                    2 / number of array dimensions
#                                            NAXIS1  =                 1508                                                  
#                                            NAXIS2  =                 1208
#                                           EXTEND  =                    T                                                  
#                                           ORIGIN  = 'Astrolab-Durham'    / Originating Institution                        
#                                           TELESCOP= '16-inch East Dome'  / Name of telescope                             
#                                            INSTRUME= 'GM2000  '           / Type of telescope                              
#                                            LONG_OBS=               1.5717 / Observatory longitude in degrees               
#                                            LAT_OBS =              54.7667 / Observatory latitude in degrees                
#                                            DETECTOR= 'G3-16200'           / Type of CCD                                   
#                                             GAIN    =                  1.0 / Default gain of detector                      
#                                              RDNOISE =                 18.0 / Default readoutof detector                   
#                                               CCD_TEMP=            -10.00611 / Temperature of CCD                             
#                                               CCD_TSET=                -10.0 / Requested temperature of CCD                   
#                                               HEATSINK=              15.7462 / Heat sink temperature of CCD                  
 #                                               EXPTIME =                180.0 / Total exposure time (seconds)                 
  #                                               FILTER  = 'V       '                                                           
   #                                               DATE-OBS= '2025-10-24T03:21:31.0484' / UTC date of the observation             
    #                                               OBJ_NAME= 'NGC1275 '           / Name of object                                
     #                                               OBSERVER= 'Astrolab'           / Name of observer                             
      ##                                                RA      =             49.98178 / Telescope RA (degrees)                      
        #                                                 DEC     =             41.47984 / Telescope Dec (degrees)                   
         ##                                                     AZIMUTH =             184.5247 / Telescope Azimuth                             
           #                                                    ALTITUDE=               76.807 / Telescope Altitude                             
            #                                                   HA      =               0.0909 / Telescope HA                                   
             ##                                                  AIRMASS =               1.0271 / Telescope Airmass                              
               #                                                TCFPOS  =                 6012 / TCF Focuser position                           
                #                                               TCFTEMP =                  5.9 / TCF temperature                                
                 ##                                              PSF_FWHM=                 3.21 / PSF FWHM in pixels                            
                   #                                             PSF_EL  =                1.073 / PSF elongation                                
                    #                                             PSF_PA  =                   46 / PSF PA in degrees                             
                     ##                               picture             RADECSYS= 'FK5     '                                                            
                       #                                          CTYPE1  = 'RA---TAN'                                                            
                        ##                                         CTYPE2  = 'DEC--TAN'                                                            
                          #                                       EQUINOX =               2000.0                                                  
                           #                                      LONPOLE =                  180                                                 
                            #                                      CRVAL1  =           49.9817833                                                 
                             ##                                      CRVAL2  =           41.4798361                                                 
                               #                                    
                                ##                                    CRPIX1  =          753.9861589                                                
                                  #                                    CRPIX2  =          603.9852305                                              
                                   #                                       CUNIT1  = 'deg     '                                                    
                                    #                                              CUNIT2  = 'deg     '                                                  
                                     #                                                       CD1_1   =       -0.00025483639                              
                                      #                                                                          CD2_1   =        3.3388918E-05          
                                       ##                                                                 CD1_2   =         3.342824E-05                                                  
                                         #                                                               CD2_2   =         0.0002548269                                                  
                                          #                                                              AST_CAL = 'Refcat2 '                                                            
                                           #                                                             N_ASTARS=                  275                                                 
                                            #                                                             RMS_X   =                0.109                                                 
                                             #                                                             RMS_Y   =                0.116                                                  
                                              #                                                            RMS_R   =                0.159                                                  
                                               #                                                           PSCALE  =               0.9252 / arcsec                                           
##
