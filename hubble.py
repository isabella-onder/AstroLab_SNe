from constants import constants
import constants as c
#given the redshift to the host galaxy (which can be inferred in other ways) & distance modulus from fitting
#infer hubble's constant. Small z assumption is made, errors are not yet accounted for

ZTF = False
#using EBV_2 and the non-sys error: need to see what it actually changes because seems to be dominated by sys
def hubble_constants(SN):
    if SN =='2019np':
        mu =  31.576
        mu_error = 0.03





        H0 = [55.82365620972563,  1.944811169891075] #H0 and error
        mu_EBV =  32.671
        mu_error_EBV =  0.037
        H0_EBV = [39.81239254518524,0.6726218011685887]
        EBV = 0.020058351568198397
        
        if ZTF:
            mu = 32.505  #bet the exact same result if I take all the data into account
            mu_error = 0.008
            H0 = [42.97523671233206, 0.1580]
            #Tmax = 58510.098  +/-  0.031  +/- 0.340 (sys)
            mu_EBV = 32.574
            mu_error_EBV = 0.007
            H0_EBV = [41.6311375,0.1339 ]
            #Tmax = 58510.183  +/-  0.044  +/- 0.340 (sys)

        return [mu, mu_error, H0, mu_EBV, mu_error_EBV, H0_EBV]
    elif SN =='2020ue':
        mu = 31.0890
        mu_error = 0.0222
        H0 = [54.66929110685479,  0.5560629925670497]
        mu_EBV = 31.068
        mu_error_EBV =  0.018
        H0_EBV = [55.20055477181215, 0.4556830508921763]
        EBV = 0.028446389496717725
        if ZTF:
            mu = 0
            mu_error = 0

            mu_EBV = 0
            mu_error_EBV = 0
        return [mu, mu_error,H0,mu_EBV, mu_error_EBV,H0_EBV]
    elif SN == '2025acsn':
        mu = 34.312
        mu_EBV = 34.004
        mu_error_EBV = 0.053
        mu_error =  0.054
        H0 = [66.05949492830327, 2.0130174347656435]
        H0_EBV = [76.12639233752994, 1.8355557634170339]
        EBV =  0.025893508388037927
        
        if ZTF:
            mu = 34.584
            mu_error = 0.043
            H0 = [58.28210767000783, 1.1427638265000013 ]
            # Tmax = 61001.197  +/-  0.309  +/- 0.340 (sys)
            mu_EBV = 34.790
            mu_error_EBV = 0.015
            H0_EBV = [53.0072, 0.3649]
            # Tmax = 60999.869  +/-  0.384  +/- 0.340 (sys)

        return [mu, mu_error,H0, mu_EBV, mu_error_EBV,H0_EBV]
    
    #20205aebt is not Ia

    #2025afwg get the following error:
    #RuntimeError: All weights for filter Is are zero. The fitter is in a part of parameter space where the model is not valid or there is no useful data.

    #2025ahkt get the following error:
    #RuntimeError: All weights for filter Bs are zero. The fitter is in a part of parameter space where the model is not valid or there is no useful data.
    #ahkt does not work on its own due to the error, but does with ztf - it worked once I deleted tht singular B band
    
    elif SN =='2025ahkt': 
        mu = 35.088
        mu_error = 0.455
        H0 = [0,0]
        mu_EBV = 35.174
        mu_error_EBV = 0.288
        H0_EBV = [74.52985729179098, 9.257350868304826]
        EBV = 0.025893508388037927

        #using both Tmax and E(B-V) to constrain - much better fit
        mu =35.216
        mu_error = 0.008
        H0 = [74.52985729179098,9.257350868304826]


        if ZTF:
            mu = 35.181
            mu_error = 0.029
            H0 = [74.28998824899324, 0.985546680206852]
            #Tmax = 61044.553  +/-  0.194  +/- 0.340 (sys)
            mu_EBV = 35.231
            mu_error_EBV = 0.009
            H0_EBV = [72.59894166128588,0.3002747363401568]
            # Tmax = 61044.439  +/-  0.186  +/- 0.340 (sys)
        return[mu,mu_error, H0, mu_EBV, mu_error_EBV, H0_EBV]

    elif SN =='2025ahqr': #it's an awful fit
        mu = 34.829
        mu_error = 0.475
        H0 = [7.573804545184397, 1.4880547488005789]
        mu_EBV =  34.375
        mu_error_EBV = 0.326
        H0_EBV = [9.334995076480887, 1.3013231336133533]
        EBV = 0.016776075857038657
        if ZTF: 
            mu = 0
            mu_error = 0
            mu_EBV = 0
            mu_error_EBV = 0
            
        return[mu,mu_error, H0, mu_EBV, mu_error_EBV, H0_EBV]
    
    
    
    elif SN =='2025ahsa':
        mu = 33.992
        mu_error = 0.049
        H0 = [72.65718304229135, 1.6211736436609954]
        mu_EBV =  33.990
        mu_error_EBV =  0.022
        H0_EBV = [72.7241336080621,0.733075623838829]
        EBV = 0.04412837345003647
        if ZTF:
            mu = 33.973
            mu_error =  0.058
            H0 = [73.29570998447538, 1.0056548420126603]
            #Tmax = 61044.553  +/-  0.194  +/- 0.340
            mu_EBV = 33.911
            mu_error = 0.030
            H0_EBV = [75.41861557965255,  0.7602366081538179]
            #Tmax_EBV =  61035.094  +/-  0.246  +/- 0.340

        return [mu, mu_error,H0, mu_EBV, mu_error_EBV, H0_EBV]
    
    elif SN =='2026acd':
        mu =  32.554
        mu_error = 0.062
        H0 = [70.67568048961604, 1.9894000493273296]
        mu_EBV =  32.636
        mu_error_EBV =  0.054
        H0_EBV = [68.05656083610121, 1.6715548554360566]
        EBV = 0.03646973012399708
        mu = 32.137
        mu_error = 0.335
        if ZTF: 
            mu = 0
            mu_error = 0
            mu_EBV = 0
            mu_error_EBV = 0
            return
        return [mu, mu_error,H0, mu_EBV, mu_error_EBV, H0_EBV]
    
    elif SN =='2026atb':
        mu = 34.100
        mu_error = 0.149
        H0 = [82.9259, 5.4993]
        mu_EBV = 0
        mu_error_EBV = 0
        H0_EBV = [0,0]
        EBV = 0.05397520058351568
        if ZTF: 
            mu = 0
            mu_error = 0
            mu_EBV = 0
            mu_error_EBV = 0
        return[mu, mu_error, H0, mu_EBV, mu_error_EBV, H0_EBV]

    elif SN =='2026gm':
        mu = 34.424
        mu_error =  0.268
        H0 = [73.69771393182003, 8.55677465070896]
        H0_EBV = [83.30172915305859,3.4187935205235647]
        mu_EBV = 34.158
        mu_error_EBV = 0.091
        EBV = 0.3931436907366886

        #using both EBVhost and Tmax to constrain, same values
        if ZTF:
            mu = 0
            mu_error = 0
            mu_EBV = 31.864
            mu_error_EBV = 0.099
            H0_EBV = [239.58190645580228, 10.63331011707412]
            #Tmax = 61068.717  +/-  1.477  +/- 0.340 (sys)
        return [mu, mu_error,H0, mu_EBV, mu_error_EBV, H0_EBV]
    
    elif SN =='2026kc':
        mu =  34.280
        mu_error = 0.128
        H0 = [97.46458121878791, 5.579115459295437]
        mu_EBV = 34.392
        mu_error_EBV = 0.137
        H0_EBV = [92.56500332589144, 5.659591672617793]
        EBV = 0.03391684901531729

        #and using constraining EBV and Tmax
        #mu =  34.405
        #mu_error =  0.09
        
        if ZTF:
            mu = 35.051
            mu_error = 0.026
            H0 = [68.33556996398131, 0.8133331011707412]
            #Tmax = 61064.781  +/-  0.102  +/- 0.340 (sys)
            mu_EBV = 35.280   
            mu_error_EBV = 0.015 
            H0 = [61.49599334668057,0.4233354405051415]
            # Tmax = 61064.476  +/-  0.174  +/- 0.340
        return [mu, mu_error,H0, mu_EBV, mu_error_EBV, H0_EBV]



#error_propagation = True
def hubble_prelim(SN): #mu is a list such that

    mu = hubble_constants(SN)
    if mu == None:
        print(f'there was an error: most likely this {SN} does not have a defined distance modulus')
        return
    type, z, _, _, _ = constants(SN)
    if type != 'Ia':
        print('NB: this supernova is not Ia!')
        
    H0 = (3*10**5 * z)*(1+z)/(10**((mu[0]-25)/5))    
    if len(mu) > 1:
        alpha = abs( (3*10**5 * z)*(1+z)/(10**((mu[0]+mu[1]-25)/5)) - (3*10**5 * z)*(1+z)/(10**((mu[0]-25)/5)))
        print(f'this was the inferred H0 for {SN} in km/s/Mpc', H0, ' p/m ', alpha)
        
    else:
        alpha = ' -- (no error was input)'
        print(f'this was the inferred H0 for {SN} in km/s/Mpc', H0)
    
    

    if len(mu)>3:
        H0_EBV = (3*10**5 * z)*(1+z)/(10**((mu[3]-25)/5)) #by taking the correction of EBV into account (fixed as a prior)
        alpha_EBV = abs( (3*10**5 * z)*(1+z)/(10**((mu[3]+mu[4]-25)/5)) - (3*10**5 * z)*(1+z)/(10**((mu[3]-25)/5)))
        print(f'this was the inferred H0 for {SN} in km/s/Mpc with E(B-V) correction', H0_EBV, ' p/m ', alpha_EBV)
    
    return [H0,alpha]

hubble_prelim('2019np')


#for SN in c.SN_list:
    #hubble_prelim(SN)
'''
def print_hubble():
    for SN in c.SN_list:
        hubbles = hubble_constants(SN)
        if hubbles == None:
            continue
        H0 = hubbles[2]
        H0_EBV = hubbles[5]
        print('######################################')
        print(f'{SN}')
        print(f'H0 = ', H0[0],'p/m', H0[1])
        print(f'with E(B-V): H0 = ', H0_EBV[0],'p/m', H0_EBV[1])
'''
#print_hubble()
SN_list_EBV = ['2025acsn', '2019np','2020ue', '2025aebt', '2025afwg', '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc', '2026dix', '2026kc_subbed', '2025acsn_subbed' ]

def extinction(SN):
    _,_,_,_,A_ned = constants(SN)
    R_V =  2.742 #from the Schlaffly paper for the Landolt V Band assuming usually R_V = 3.1 as they suggest
    E_BV = A_ned / R_V
    print(f'this is the E(B-V) for the host galaxy of {SN}', E_BV)
    return E_BV

large_EBV =  []
SN_list_EBV = ['2025acsn', '2019np', '2020ue', '2025ahkt',
            '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc',
            '2026kc_subbed', '2025acsn_subbed']
for SN in SN_list_EBV:
    EBV = extinction(SN) 
    large_EBV.append(EBV)
print(large_EBV) 
extinction('2025aebt')