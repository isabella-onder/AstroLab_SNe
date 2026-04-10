import pandas as pd
import numpy as np

SN = '2019np'


#for when I want to get the Hubble constant for the data done with ZTF
# i.e. both the raw ZTF data and from my own data but constrained

ZTF_list = ['2019np', '2025acsn', '2025ahkt',  '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc' ] #all those which have been observed by ZTF too
my_SNe = ['2019np','2020ue', '2025acsn',  '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc'] #all those which we observed
my_SNe_new = ['2019np','2020ue', '2025ahsa',  '2026acd', '2026atb', '2026gm', '2026kc_subbed', '2025acsn_subbed', '2025ahkt'] #all those which we observed


ZTF = False
EBV = False
dm15 = False
st = False
completely_automatic = True


def hubble(DM, DM_error,z, SN_specific): #here, z is z_CMB

        
    H0 = (3*10**5 * z)*(1+z)/(10**((DM-25)/5))    

    alpha = abs( (3*10**5 * z)*(1+z)/(10**((DM+DM_error-25)/5)) - (3*10**5 * z)*(1+z)/(10**((DM-25)/5)))
    print(f'this was the inferred H0 for {SN_specific} in km/s/Mpc', H0, ' p/m ', alpha)
    return [H0,alpha]

def weighted_mean(values, uncertainties):
    weights = 1 / np.array(uncertainties)**2
    weighted_avg = np.sum(weights * np.array(values)) / np.sum(weights)
    uncertainty  = 1 / np.sqrt(np.sum(weights))
    return weighted_avg, uncertainty




#for ZTF data, row is as follows:
#zcmb    DM     DM error    Tmax    Tmax error   
#0       1         2         3          4
# title: dm15data    DM     DM error    dm15    dm15 error  
#5                   6         7         8        9       
# title:st data      DM      DM error    st     st error
#10                  11          12       13         14

#for my data, rows as follows
#zcmb   dm     dmerror      dm15       dm15error  tmax    tmax_error
#0        1      2           3         4            5       6
#title for st   dm  dmerror   st   st error  tmax   tmax_error
#7              8     9        10   11          12      13          

if ZTF: 
    df = pd.read_excel('SNe_analysis.xlsx', sheet_name = 'ZTF')
    ZTF_hubbles = []
    ZTF_errors = []
    for ztf_SN in ZTF_list:
        #df = pd.read_excel('SNe_analysis.xlsx', sheet_name = 'ZTF')
        ztf_outputs = df[ztf_SN].tolist()
        pure_ZTF_hubble = hubble(ztf_outputs[1],ztf_outputs[2],ztf_outputs[0], ztf_SN)
        ZTF_hubbles.append(pure_ZTF_hubble[0])
        ZTF_errors.append(pure_ZTF_hubble[1])
    mean_hubble = np.mean(ZTF_hubbles)
    weighted_hubble, weighted_uncertainty = weighted_mean(ZTF_hubbles, ZTF_errors)
    print('############## THE MEAN AND WEIGHTED H0 for this method are', mean_hubble, 'and', weighted_hubble, 'p/m', weighted_uncertainty)

if EBV:
    df = pd.read_excel('SNe_analysis.xlsx', sheet_name = 'EBV_fixed')
    EBV_hubbles = []
    EBV_errors = []

    for my_SN in my_SNe:
        my_outputs = df[my_SN].tolist()
        if dm15:
            EBV_dm15_my_hubble = hubble(my_outputs[1], my_outputs[2], my_outputs[0], my_SN)
            EBV_hubbles.append(EBV_dm15_my_hubble[0])
            EBV_errors.append(EBV_dm15_my_hubble[1])
        if st:
            EBV_st_my_hubble = hubble(my_outputs[8], my_outputs[9], my_outputs[0], my_SN)
            EBV_hubbles.append(EBV_st_my_hubble[0])
            EBV_errors.append(EBV_st_my_hubble[1])
    mean_hubble = np.mean(EBV_hubbles)
    weighted_hubble, weighted_uncertainty = weighted_mean(EBV_hubbles, EBV_errors)
    print('############## THE MEAN AND WEIGHTED H0 for this method are', mean_hubble, 'and', weighted_hubble, 'p/m', weighted_uncertainty)
        
 
if completely_automatic:
    df = pd.read_excel('snpy_txt_clean/snoopy_results_EBV2_st.xlsx')
    EBV_hubbles = []
    EBV_errors = []

    for my_SN in my_SNe_new:
        my_SN = 'SN'+my_SN
        my_outputs = df[my_SN].tolist()
        EBV_dm15_my_hubble = hubble(my_outputs[1], my_outputs[2], my_outputs[0], my_SN)
        EBV_hubbles.append(EBV_dm15_my_hubble[0])
        EBV_errors.append(EBV_dm15_my_hubble[1])
        
    mean_hubble = np.mean(EBV_hubbles)
    mean_error_list = [error**2 for error in EBV_errors]
    mean_error = np.sqrt(np.sum(mean_error_list))
    weighted_hubble, weighted_uncertainty = weighted_mean(EBV_hubbles, EBV_errors)
    print('############## THE MEAN AND WEIGHTED H0 for this method are', mean_hubble, 'p/m', mean_error,'and', weighted_hubble, 'p/m', weighted_uncertainty)
        
    

         