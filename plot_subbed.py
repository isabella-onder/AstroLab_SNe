#creating plots from the svc files

import csv
import numpy as np
import matplotlib.pyplot as plt
import mag_sub as ms

#the column structure in the csv file is as follows
#file	Telescope	MJD	  date	Filter	RA/deg	Dec/deg	#_stars	airmass	zpt_mean  zpt_err	magnitude	mag_err_upper	mag_err_lower zpt_median counts_phot counts_phot_error psf
#0        1          2     3      4        5       6       7      8       9          10        11             12            13            14          15            16          17

#the current supernovae that we have analysed are
#2025acsn, 2025afwg, 2025ahqr, 2025ahxd, 2026acd, 2026gm, 2026kc

SNes = ['2026kc', '2026kc_subbed'] #here, put the name of the supernova to analyse
zpt_plot_tick = False
star_nb = False
plain_zpt = False
psf = False
SN = SNes[0]
SNe = SNes[0]

def plot(SN):
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
        
    




    with open(f'../Photometry_final/{SN}.csv', mode = 'r') as file:
    #with open(f'../Photometry_dark_magic/mosaic_photometry_table_{SN}.csv', mode = 'r') as file:
        table = csv.reader(file)
        next(table) #discarding the first row to avoid the headers


        
        v = []
        b = []
        r = []
        i = []

        zpt_medians = []
        zpt_means = []
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
            mag_errors_down.append(float(stacked[13]))
            mag_errors_up.append(float(stacked[12]))
            stars.append(float(stacked[7]))
            days.append(float(stacked[2]))
            psfs.append(float(stacked[17]))

        if star_nb:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print(tot_error, 'this is tot error')
        
            plt.scatter(stars, tot_error)
            plt.xlabel('nb of stars')
            plt.ylabel('Error on magnitude')
            plt.show()
            return
        
        if zpt_plot_tick:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            print(tot_error, 'this is tot error')
            mean_median = [abs(mean - median) for mean, median in zip(zpt_means, zpt_medians)]
            print(mean_median, 'this is mean median')
            print(np.argmax(mean_median))
            plt.scatter(mean_median, tot_error)
            plt.xlabel('|Zeropoint mean - median|')
            plt.ylabel('Error on magnitude')
            plt.show()
            return 

        if plain_zpt:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            plt.scatter(zpt_medians, tot_error)
            plt.xlabel('Zeropoint median')
            plt.ylabel('Error on magnitude')
            plt.show()
            return
        
        if psf:
            tot_error = [np.sqrt(m_1**2+m_2**2) for m_1, m_2 in zip(mag_errors_down, mag_errors_up)]
            plt.scatter(psfs, tot_error)
            plt.xlabel('psf')
            plt.ylabel('Error on magnitude')
            plt.show()
            return
        
        #structure for each element in overarching (i.e. each observation): MJD, magnitude, mag_err_upper, mag_err_lower

        overarching = [i,r,v,b]
        
    return overarching, overarching_ztf

colours = ['firebrick', 'indianred', 'mediumseagreen', 'cornflowerblue' ]
labels = ['I', 'R', 'V', 'B']
colours_ztf = ['darkseagreen', 'rosybrown', 'grey']
labels_ztf = ['g', 'r', 'i']

overarching_1, overarching_ztf = plot(SNes[0])
overarching_2, overarching_ztf = plot(SNes[1])

overarchings = [overarching_1, overarching_2]

over_colours = [['grey', 'grey', 'grey', 'grey'], colours]   

index = 0
fig, axs = plt.subplots(2,2)
axs_list = [axs[0,0], axs[0,1], axs[1,0],axs[1,1]]
for overarching,colours in zip(overarchings,over_colours):
    index = 0
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








#could make a rejection threshold on the basis of zpt (e.g. if it deviates by x amount from the mean zeropoint across all observations)
#rather do a rejection based on mean vs median of zeropoint
