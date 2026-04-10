from astropy.io import fits
from astropy.wcs import WCS 
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
import numpy as np


before_filename = "differences/amosaic2019np_19_01_13_R.fits"
before_info = "Before: 2019 Typical R band"
after_filename = "differences/amosaic2019np_20_10_26_R.fits"
after_info = "After: 20_10_26 R band "
after_aligned = after_aligned * 1.5 #scaling by hand after_aligned


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

    fig, axs = plt.subplots(nrows =1, ncols = 2)
    #plot before
    norm = ImageNormalize(before_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_data, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info+' RAW')

    
    #plot after
    norm = ImageNormalize(after_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[1].imshow(after_data, cmap='gray', origin='lower', norm=norm)
    #axs[1].colorbar()
    fig.colorbar(im, ax=axs[1])
    axs[1].set_xlabel(after_info+' RAW')

    plt.show()



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
    difference = before_data - after_aligned #subtracting the two arrays

    

    fig, axs = plt.subplots(nrows =1, ncols = 3)
    #fig2, ax2 = plt.subplots()
    #plot before
    norm = ImageNormalize(before_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = axs[0].imshow(before_data, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig.colorbar(im, ax=axs[0])
    axs[0].set_xlabel(before_info)

    '''
    norm = ImageNormalize(before_data, interval=ZScaleInterval()) #such that can actually see features, otherwise grey
    im = ax2.imshow(before_data, cmap='gray', origin='lower', norm=norm)
    #axs[0].colorbar()
    fig2.colorbar(im, ax=axs[0])
    ax2.set_xlabel(before_info)
    fig2.savefig("Differences_before", dpi = 150, bbox_inches = 'tight')
    '''




    
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