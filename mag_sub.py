import numpy as np

#relevant stars:
# 2020ue: g = 18.3908 ± 0.0011, b = 18.6739 ± 0.0182, r = 17.9513 ± 0.0249
# 2025acsn : g = 19.5104 ± 0.0024, b = 19.9801 ± 0.0564	mag, r= 18.5000 ± 0.0455

acsn_g = 19.5104
acsn_r = 18.5000

ue_g = 18.3908
ue_r = 17.9513

def gaia_to_JC(gaia_gmag, gaia_rmag, band):
#################################################################################
#inputs:
#gaia_mags: magnitude of star as given in the gaia database: g or r
#band: band in which that magnitude we want to output - Johnsons V,B or Cousins R (no I since gaia does not produce I)

#outputs:
#mag_star: magnitude of star in our Johnson Counsin band, using the transformations in the astrolab code
################################################################################
    if band == 'V':
        ps_v = 0.006 + 0.474 * (gaia_gmag - gaia_rmag) + gaia_rmag 
        return ps_v
    elif band == 'B':
        ps_b = 0.213 + 0.587 * (gaia_gmag - gaia_rmag) + gaia_gmag 
        return ps_b
    elif band =='R':
        ps_r_c = -0.138 - 0.131 * (gaia_gmag - gaia_rmag) + gaia_rmag 
        return ps_r_c
    else:
        print('the requested band does not have a transfomation', band)



def mag_sub(photo_mag, star_mag_gaia_r, star_mag_gaia_v,zeropoint_mag, band):
#################################################################################
#inputs:
#photo_mag : magnitude as given by the photometry directly (with default radius) 
#star_mag_gaia_r/v : magnitude of the start within the aperture, as given by GAIA, in r or v band
#zeropoint_mag : magnitude of the zeropoint on that given day ( = 1 count)
#band: V, B or R which we want to obtain in Johsons VB or Cousins IR, from Gaia

#output
#sn_mag: the magnitude of the supernova without the star
#################################################################################

    #convert the star mag from gaia to 
    star_mag = gaia_to_JC(star_mag_gaia_v, star_mag_gaia_r, band)


    #calculate the count contribution from the star
    star_cts = 10**((zeropoint_mag - star_mag)/2.5)

    #convert the total photometric magnitude into counts
    total_cts = 10**((zeropoint_mag - photo_mag)/2.5)

    #subtract star from total for SN contribution only (background already accounted for)
    sn_cts = total_cts - star_cts

    #calculate magnitude from SN contribution only - will be lower
    sn_mag = zeropoint_mag - 2.5*np.log10(sn_cts)
    #print('This is the difference in magnitude (after - before)', sn_mag - photo_magz)
    return sn_mag






#need to do this for 2020ue and 2025acsn i.e. when star in aperture but no subtraction is being done