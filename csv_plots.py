#creating plots from the svc files

import csv
import numpy as np
import matplotlib.pyplot as plt
import mag_sub as ms

from matplotlib.ticker import AutoMinorLocator


SN_list = ['2025acsn', '2019np','2020ue', '2025aebt', '2025afwg', '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc', '2026dix' ]


#the column structure in the csv file is as follows
#file	Telescope	MJD	  date	Filter	RA/deg	Dec/deg	#_stars	airmass	zpt_mean  zpt_err	magnitude	mag_err_upper	mag_err_lower zpt_median counts_phot counts_phot_error fwhm(pix) fwhm(arcsecs) platescale
#0        1          2     3      4        5       6       7      8       9          10        11             12            13            14          15            16          17             18           19

#the current supernovae that we have analysed are
#2025acsn, 2025afwg, 2025ahqr, 2025ahxd, 2026acd, 2026gm, 2026kc

SNe = '2025acsn' #here, put the name of the supernova to analyse
zpt_plot_tick = False
star_nb = False
plain_zpt = False
psf = False

report_plot = True

mega_subplot = False

def plot(SN, plot):
    try:
        with open(f'../ZTF_tables/{SN}.csv') as file:
            table = csv.reader(file)
            next(table)

            g_ztf = []
            r_ztf = []
            i_ztf = []

            row_nb = 1
            for row in table:
                row_nb += 1
                #print('this is the problematic row 3', row[3], 'see how many spaces')
                #only keep ztf observations with errors
                if row[1] == 'g':
                    try:
                        g_mag_info = [float(row[0]), float(row[2]), float(row[3])]
                        g_ztf.append(g_mag_info)
                    except ValueError: 
                        print('the stack on row', row_nb, 'is most certainly missing something')
                        continue
                elif row[1] == 'r':
                    try:
                        r_mag_info = [float(row[0]), float(row[2]), float(row[3])]
                        r_ztf.append(r_mag_info)
                    except ValueError:
                        print('the stack on row', row_nb, 'is most certainly missing something')
                        continue
                elif row[1] == 'i':
                    try:
                        i_mag_info = [float(row[0]), float(row[2]), float(row[3])]
                        i_ztf.append(i_mag_info)
                    except ValueError:
                        print('the stack on row', row_nb, 'is most certainly missing something')
                        continue
                else:
                    print('the band is not in g or r')

            #info in each observation which then goes into overargcing: MJD (time), mag, mag_error
            overarching_ztf = [g_ztf, r_ztf, i_ztf]
            colours_ztf = ['darkseagreen', 'rosybrown', 'grey']
            labels_ztf = ['g', 'r', 'i']
    except FileNotFoundError:
        print('there is no ztf file for this supernova')
        overarching_ztf, colours_ztf, labels_ztf = [], [], []
        
    




    with open(f'../Photometry_final_uncleaned_backup/{SN}.csv', mode = 'r') as file:
    #with open(f'../Photometry_dark_magic/mosaic_photometry_table_{SN}.csv', mode = 'r') as file:
        table = csv.reader(file)
        next(table) #discarding the first row to avoid the headers


        
        v = []
        b = []
        r = []
        i = []

        zpt_medians = []
        zpt_means = []
        zpt_errors=[]
        mag_errors_up = []
        mag_errors_down = []
        stars = []
        tot_error_v = []
        tot_error_b = []
        tot_error_r = []
        tot_error_i = []
        days = []
        psfs = []

        row_nb = 1
        for stacked in table:  #each line is a different stacked image
            row_nb += 1
            #print('this is what is meant to be the band', stacked[4])


            #if incomplete data, do not include
            if stacked[11] == '' or stacked[12] == '' or stacked[13] == '':
                print('the stack on line', row_nb, 'has no magnitude(or mag error), it will be disregarded: \n', stacked)
                continue
            #sort according to band
            if stacked[4] == 'V':
                v_mag_info = [float(stacked[2]), float(stacked[11]), float(stacked[12]),float(stacked[13])] 
                error_v = [np.sqrt(float(stacked[12])**2+float(stacked[13])**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
                if SNe == '2020ue':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.ue_r, ms.ue_g, float(stacked[14]),'V')
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])]
                if SNe == '2025acsn':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.acsn_r, ms.acsn_g, float(stacked[14]),'V') 
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])] 
                tot_error_v.append(error_v)
                v.append(v_mag_info)
            elif stacked[4] == 'B':
                b_mag_info = [float(stacked[2]), float(stacked[11]), float(stacked[12]),float(stacked[13])] 
                error_b = [np.sqrt(float(stacked[12])**2+float(stacked[13])**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
                if SNe == '2020ue':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.ue_r, ms.ue_g, float(stacked[14]),'B')
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])]
                if SNe == '2025acsn':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.acsn_r, ms.acsn_g, float(stacked[14]),'B') 
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])]  
                tot_error_b.append(error_b)
                b.append(b_mag_info)
            elif stacked[4] == 'R' or stacked[4] == 'r': #this is to account for both our .5 and east16 measurements
                r_mag_info = [float(stacked[2]), float(stacked[11]), float(stacked[12]),float(stacked[13])] 
                error_r = [np.sqrt(float(stacked[12])**2+float(stacked[13])**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
                if SNe == '2020ue':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.ue_r, ms.ue_g, float(stacked[14]),'R')
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])]
                if SNe == '2025acsn':
                    subtracted_mag = ms.mag_sub(float(stacked[11]), ms.acsn_r, ms.acsn_g, float(stacked[14]),'R') 
                    v_mag_info = [float(stacked[2]), subtracted_mag, float(stacked[12]),float(stacked[13])]  
                tot_error_r.append(error_r)
                r.append(r_mag_info)
            elif stacked[4] == 'I':
                i_mag_info = [float(stacked[2]), float(stacked[11]), float(stacked[12]),float(stacked[13])]  
                error_i = [np.sqrt(float(stacked[12])**2+float(stacked[13])**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
                tot_error_i.append(error_i)
                i.append(i_mag_info)
            else:
                print('the band is not V, B, R or I')
            
            zpt_medians.append(float(stacked[14]))
            zpt_means.append(float(stacked[9]))
            zpt_errors.append(float(stacked[10]))
            mag_errors_down.append(float(stacked[13]))
            mag_errors_up.append(float(stacked[12]))
            stars.append(float(stacked[7]))
            days.append(float(stacked[2]))
            psfs.append(float(stacked[18]))

            

            zpt_mean_of_means = np.mean(zpt_means)
            zpt_std_of_means = np.std(zpt_means)

            fwhm_mean = np.mean(psfs)
            fwhm_std = np.std(psfs)

        if mega_subplot:
            fig, axs = plt.subplots(2,2)
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print(tot_error, 'this is tot error')
        
            axs[0][0].scatter(stars, tot_error, c= psfs, cmap = 'RdYlGn_r')
            axs[0][0].set_xlabel('nb stars')

            mean_median = [abs(mean - median) for mean, median in zip(zpt_means, zpt_medians)]
            mean_median_raw = [mean - median for mean, median in zip(zpt_means, zpt_medians)]
            mean_difference = np.mean(mean_median_raw)
            difference_std = np.std(mean_median_raw)
            print('this is the mean and std of the zpt', mean_difference, difference_std)
            print('and this is 3 sigma difference',mean_difference-3*difference_std,mean_difference+3*difference_std )
            axs[0][1].scatter(mean_median, tot_error,c= psfs, cmap = 'RdYlGn_r')
            axs[0][1].set_xlabel('Zpt mean - median')

            print('this is the mean and std of the plain zpt', zpt_mean_of_means, zpt_std_of_means)
            print('and this is 3 sigma difference',zpt_mean_of_means-3*zpt_std_of_means,zpt_mean_of_means+3*zpt_std_of_means )

            axs[1][0].scatter(zpt_medians, tot_error, c= psfs, cmap = 'RdYlGn_r')
            axs[1][0].set_xlabel('plain zpt')

            print('this is the mean and std of the psf', fwhm_mean, fwhm_std)
            print('and this is 3 sigma difference',fwhm_mean-3*fwhm_std,fwhm_mean+3*fwhm_std )
            axs[1][1].scatter(psfs, tot_error)
            axs[1][1].set_xlabel('psf')
            plt.show()
            #plt.colorbar()

            return
        
        '''
        if report_plot:
            fig, axs = plt.subplots(ncols=1, nrows=3, sharex=True, figsize=(3, 4.5))
            tot_error = [np.sqrt(m_1**2 + m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]

            mean_median = [abs(mean - median) for mean, median in zip(zpt_means, zpt_medians)]
            mean_median_raw = [mean - median for mean, median in zip(zpt_means, zpt_medians)]
            mean_difference = np.mean(mean_median_raw)
            difference_std = np.std(mean_median_raw)

            vmin, vmax = min(psfs), max(psfs)  # shared colour scale across all panels

            S = 27
            from matplotlib.colors import LinearSegmentedColormap
            rg = LinearSegmentedColormap.from_list('rg', ['mediumseagreen', 'firebrick'])
            #RdYlGn_r'
            rg = 'YlOrRd'
            edge = 'silver'
            sc0 = axs[0].scatter(tot_error, zpt_medians,   c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s = S, edgecolor = edge)
            sc1 = axs[1].scatter(tot_error, mean_median,   c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s = S, edgecolor = edge)
            sc2 = axs[2].scatter(tot_error, psfs,          c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s = S, edgecolor = edge)

            axs[0].set_ylabel('Zeropoint')
            axs[1].set_ylabel('Zpt mean - median')
            axs[2].set_ylabel("FWHM ('')")
            axs[2].set_xlabel('Error in magnitude')

            

            for ax in axs:
                ax.minorticks_on()
                ax.tick_params(which='both', direction='in', labelsize=9)
                ax.tick_params(which='major', length=5, width=1.5)
                ax.tick_params(which='minor', length=2, width=1.0)
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_major_locator(plt.MaxNLocator(5))
                ax.yaxis.set_label_position('right')
                ax.tick_params(which='both', left=True, right=True, labelleft=False, labelright=True)   
            

            # Single colourbar for all three panels, placed to the right
            plt.colorbar(sc0, ax=axs[0], label='FWHM',  location = 'top' )

            fig.tight_layout()
            plt.savefig('Figures/acsn_outliers.png', bbox_inches = 'tight', dpi = 150)
            plt.show()
            '''
        if report_plot:
            fig, axs = plt.subplots(ncols=1, nrows=3, sharex=True, figsize=(4, 6))
            tot_error = [np.sqrt(m_1**2 + m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]

            mean_median = [abs(mean - median) for mean, median in zip(zpt_means, zpt_medians)]
            mean_median_raw = [mean - median for mean, median in zip(zpt_means, zpt_medians)]
            mean_difference = np.mean(mean_median_raw)
            difference_std = np.std(mean_median_raw)

            vmin, vmax = min(psfs), max(psfs)

            S = 38
            from matplotlib.colors import LinearSegmentedColormap
            rg = 'YlOrRd'
            edge = 'silver'
            sc0 = axs[0].scatter(tot_error, zpt_medians, c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s=S, edgecolor=edge)
            sc1 = axs[1].scatter(tot_error, mean_median, c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s=S, edgecolor=edge)
            sc2 = axs[2].scatter(tot_error, psfs,        c=psfs, cmap=rg, vmin=vmin, vmax=vmax, s=S, edgecolor=edge)

            # --- Zeropoint panel: mean and 3σ shading ---


            zpt_mean_val = 25.2951
            zpt_std_val  = 0.4765
            axs[0].axhline(zpt_mean_val, color='silver', linewidth=2, linestyle='--', zorder=3)
            #axs[0].axhspan(zpt_mean_val + 3*zpt_std_val, axs[0].get_ylim()[1],
            #            color='grey', alpha=0.45, zorder=0)
            axs[0].axhspan(axs[0].get_ylim()[0], zpt_mean_val - 3*zpt_std_val,
                        color='silver', alpha=0.45, zorder=0)

            # --- Mean–median panel: mean and 3σ shading ---
            mm_mean_val = 0.00142
            mm_std_val  = 0.04765
            axs[1].axhline(mm_mean_val, color='grey', linewidth=2, linestyle='--', zorder=3)
            axs[1].axhspan(mm_mean_val + 3*mm_std_val, axs[1].get_ylim()[1],
                        color='silver', alpha=0.45, zorder=0)
            #axs[1].axhspan(axs[1].get_ylim()[0], mm_mean_val - 3*mm_std_val,
             #           color='steelblue', alpha=0.15, zorder=0)

            labelsi = 14
            padsi = 7
            axs[0].set_ylabel('Zeropoint', size = labelsi,labelpad = padsi)
            axs[1].set_ylabel('Zpt mean - median', size = labelsi,labelpad = padsi)
            axs[2].set_ylabel("FWHM ('')", size = labelsi, labelpad = padsi +4)
            axs[2].set_xlabel('Error in magnitude', size = labelsi, labelpad = padsi)

            for ax in axs:
                ax.minorticks_on()
                ax.tick_params(which='both', direction='in', labelsize=labelsi)
                ax.tick_params(which='major', length=5, width=1.5)
                ax.tick_params(which='minor', length=2, width=1.0)
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_major_locator(plt.MaxNLocator(5))
                ax.yaxis.set_label_position('right')
                ax.tick_params(which='both', left=True, right=True, labelleft=False, labelright=True)

            plt.colorbar(sc0, ax=axs[0], label='FWHM', location='top')

            fig.tight_layout()
            plt.savefig('Figures/acsn_outliers.png', bbox_inches='tight', dpi=150)
            plt.show()


        if star_nb:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print(tot_error, 'this is tot error')
        
            plt.scatter(stars, tot_error, c= psfs, cmap = 'RdYlGn_r')
            plt.colorbar()
            plt.xlabel('nb of stars')
            plt.ylabel('Error on magnitude')
            plt.show()


            return
        
        if zpt_plot_tick:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print(tot_error, 'this is tot error')
            mean_median = [abs(mean - median) for mean, median in zip(zpt_means, zpt_medians)]
            mean_median_raw = [mean - median for mean, median in zip(zpt_means, zpt_medians)]

            mean_difference = np.mean(mean_median_raw)
            difference_std = np.std(mean_median_raw)

            print(np.argmax(mean_median))

            print('this is the mean and std of the zpt', mean_difference, difference_std)
            print('and this is 3 sigma difference',mean_difference-3*difference_std,mean_difference+3*difference_std )

            #plt.errorbar(mean_median, tot_error,xerr=zpt_errors,fmt='none',color='grey',alpha=0.5,zorder=1)
            plt.scatter(mean_median, tot_error,c= psfs, cmap = 'RdYlGn_r')
            plt.colorbar()
            plt.xlim(xmin = 0)
            plt.xlabel('|Zeropoint mean - median|')
            plt.ylabel('Error on magnitude')
            plt.show()
            return 

        if plain_zpt:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            #plt.errorbar(zpt_medians, tot_error,xerr=zpt_errors,fmt='none',color='grey',alpha=0.5,zorder=1)
            print('this is the mean and std of the zpt', zpt_mean_of_means, zpt_std_of_means)
            print('and this is 3 sigma difference',zpt_mean_of_means-3*zpt_std_of_means,zpt_mean_of_means+3*zpt_std_of_means )

            plt.scatter(zpt_medians, tot_error, c= psfs, cmap = 'RdYlGn_r')
            plt.colorbar()
            #plt.xlim(xmin = 0)
            plt.xlabel('Zeropoint median')
            plt.ylabel('Error on magnitude')
            plt.show()
            return
        
        if psf:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print('this is the mean and std of the psf', fwhm_mean, fwhm_std)
            print('and this is 3 sigma difference',fwhm_mean-3*fwhm_std,fwhm_mean+3*fwhm_std )

            plt.scatter(psfs, tot_error)
            plt.xlabel('psf')
            plt.ylabel('Error on magnitude')
            plt.show()
            return
        
        #structure for each element in overarching (i.e. each observation): MJD, magnitude, mag_err_upper, mag_err_lower

        overarching = [i,r,v,b]
        colours = ['firebrick', 'indianred', 'mediumseagreen', 'cornflowerblue' ]
        labels = ['I', 'R', 'V', 'B']



                

        

    index = 0
    fig, axs = plt.subplots(2,2)
    axs_list = [axs[0,0], axs[0,1], axs[1,0],axs[1,1]]
    for band,axs in zip(overarching,axs_list):
        #print(band, index)
        if band == []:
            print('there were no observations in band', labels[index])
            continue
        time, mags, mag_error_up, mag_error_down = map(list, zip(*band)) #turns rows into columns, and map(list) is simply to ensure it is a list rather than tuple (otherwise just = zip(*band))
        mags = [m * (-1) for m in mags] # such as to get the y axis the correct way up
        axs.errorbar(time, mags,(mag_error_up, mag_error_down), fmt = '.',color = colours[index], label = labels[index])


        if labels[index] == 'V' and overarching_ztf != []:
            try:
                band_ztf = overarching_ztf[0]
                time_ztf, mag_ztf, mag_error_ztf = map(list, zip(*band_ztf))
                mag_ztf = [-m for m in mag_ztf]
                axs.errorbar(time_ztf, mag_ztf, mag_error_ztf, fmt = '.', color = colours_ztf[0], label = labels_ztf[0])
            except ValueError: #when there is no observation in that band from ztf but there are still some ztf values that were in fact taken
                print('oops, no ztf g band')
        if labels[index] == 'R' and overarching_ztf != []:
            try:
                band_ztf = overarching_ztf[1]
                time_ztf, mag_ztf, mag_error_ztf = map(list, zip(*band_ztf))
                mag_ztf = [-m for m in mag_ztf]
                axs.errorbar(time_ztf, mag_ztf, mag_error_ztf, fmt = '.', color = colours_ztf[1], label = labels_ztf[1])
            except ValueError: #when there is no observation in that band from ztf but there are still some ztf values that were in fact taken
                print('oops, no ztf r band')
        if labels[index] == 'I' and overarching_ztf != []:
            try: 
                band_ztf = overarching_ztf[2]
                time_ztf, mag_ztf, mag_error_ztf = map(list, zip(*band_ztf))
                mag_ztf = [-m for m in mag_ztf]
                axs.errorbar(time_ztf, mag_ztf, mag_error_ztf, fmt = '.', color = colours_ztf[2], label = labels_ztf[2])
            except ValueError: #when there is no observation in that band from ztf but there are still some ztf values that were in fact taken
                print('oops, no ztf i band')
        
        axs.set_xlabel('Modified Julian Date')
        axs.set_ylabel('Apparent magntitude')
        #axs.invert_yaxis()
        index +=1
    fig.suptitle(f'SN{SN}')
    plt.legend()
        
    if plot:
        plt.show()
    #print (overarching)
    return overarching, overarching_ztf


plot(SNe, True)

    



        
#could make a rejection threshold on the basis of zpt (e.g. if it deviates by x amount from the mean zeropoint across all observations)
#rather do a rejection based on mean vs median of zeropoint
        