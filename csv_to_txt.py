#converting it into the correct format for SNoopy
import csv_plots as cp
import constants as c
import numpy as np

SN_list = ['2025acsn', '2019np','2020ue', '2025aebt', '2025afwg', '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc', '2026dix', '2026kc_subbed', '2025acsn_subbed']
SN_list = ['2025advo']

SN = 0 #put supernova name
def text(SN):
    overarching, overarching_ztf = cp.plot(SN, False)
    bands = ['I', 'R', 'V', 'B']
    bands_ztf = ['g','r', 'i']
    type, redshift, RA, DEC, _ = c.constants(SN)
    print(redshift)
    name = 'SN'+SN

    #need to modify the intaken overarching such that it has the data that I want
    #need to put in a place at the top, like a constants file with all the info for all the SN

    with open(f"snpy_txt_clean/{SN}.txt", "w") as f:
        f.write(f"{name} {redshift} {RA} {DEC}\n") #this is the header
        for band, band_name in zip(overarching, bands):
            f.write(f'filter {band_name}s \n')
            for stack in band:
                error = np.sqrt(stack[2]**2+stack[3]**2)
                new_stack = [stack[0], stack[1], error]
                line = " ".join(str(element) for element in new_stack)
                f.write(line + "\n")
    print(f'### The csv_to_txt.py file was run for {SN}. ###')

    '''
    #the same for ztf
    with open(f"snpy_txt_ztf/output_{SN}.txt", "w") as f:
        f.write(f"{name} {redshift} {RA} {DEC}\n") #this is the header
        for band, band_name in zip(overarching_ztf, bands_ztf):
            f.write(f'filter {band_name} \n')
            for stack in band:
                #time_ztf, mag_ztf, mag_error_ztf = map(list, zip(*band_ztf))
                new_stack = [stack[0], stack[1], stack[2]]
                line = " ".join(str(element) for element in new_stack)
                f.write(line + "\n")
    print(f'### The csv_to_txt.py file was run for ZTF data of {SN}. ###')
    '''

#text('2019np')
for SN in SN_list:
    text(SN)