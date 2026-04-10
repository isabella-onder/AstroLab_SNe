#code to regroup all the constants necessary for the snoopy text file header
#info comes from TNS
SN_list = ['2019np','2020ue', '2025acsn', '2025advo', '2025aebt', '2025afwg', '2025ahkt', '2025ahqr', '2025ahsa', '2026acd', '2026atb', '2026gm', '2026kc', '2026dix' ]


#for the time being taking A_lambda Landolt V as host_NED --> this is correct

#spouts out type, redshit, coordinates, extinction
#adding also the CMB redshift which is given with cz
def constants(SN):
    type, z, RA, DEC, host_NED, discovery_mag = '', 0, 0 , 0 , 0, 0 #initialising
    if SN == '2025acsn':
        type = 'Ia'
        z = 0.015791
        RA = 346.8283843
        DEC = 42.1988628
        discovery_mag = 19.3635
        host_NED = 0.435 #CGCG 532-008 NED01 not in TNS but ZTF
        return type, z, RA, DEC, host_NED
    if SN=='2025acsn_subbed':
        type = 'Ia'
        z = 0.015791
        RA = 346.8283843
        DEC = 42.1988628
        discovery_mag = 19.3635
        host_NED = 0.435 #CGCG 532-008 NED01 not in TNS but ZTF
        return type, z, RA, DEC, host_NED
    elif SN == '2025ahxd':
        type = 'Ia'
        z = 0.016455
        RA = 12.9005888
        DEC = 29.711097590909
        discovery_mag = 20.1472
        host_NED = 0
        return type, z, RA, DEC, host_NED
    elif SN == '2019np':
        type = 'Ia'
        z = 0.00452
        RA = 157.341583
        DEC = 29.510639
        discovery_mag = 17.8
        host_NED = 	0.055 #NGC 3254
        return type, z, RA, DEC, host_NED
    elif SN == '2025advo':
        type = 'Ic'
        z = 0.014083
        RA = 124.416181802
        DEC = 4.61476787169
        discovery_mag = 16.895
        host_NED = 0.057 #UGC 4316
        return type, z, RA, DEC, host_NED
    elif SN == '2025aebt':
        type = 'Ib'
        z = 0.017649
        RA = 107.3281207
        DEC = 20.6387489
        discovery_mag = 19.173
        host_NED = 	0.198 #ngc 2342 -- not on TNS but in ZTF
        return type, z, RA, DEC, host_NED
    elif SN == '2025afwg':
        type = 'Ic'
        z = 0.016621
        RA = 42.8178164 
        DEC = 50.8743039
        discovery_mag = 20.0512
        host_NED = 1.054 # or as checked by Isabella by eye: CGCG 554-012
        return type, z, RA, DEC, host_NED
    elif SN == '2025ahkt':
        type = 'Ia'
        z = 0.026228
        RA = 141.426836087 
        DEC = 22.3349335633
        discovery_mag = 20.6289
        host_NED = 0.071 #CGCG 121-087
        return type, z, RA, DEC, host_NED
    elif SN == '2025ahqr':
        type = 'Ia' #Iax[02cx-like]
        z = 0.002328
        RA = 167.858238192 
        DEC = 55.6749353975
        discovery_mag = 20.0269
        host_NED = 0.046 # M108
        return type, z, RA, DEC, host_NED
    elif SN == '2025ahsa':
        type = 'Ia'
        z = 0.015
        RA = 218.750875871 
        DEC = 0.986217887085
        discovery_mag = 14.9301
        host_NED = 	0.121 #WISEA J143459.95+005908.8
        return type, z, RA, DEC, host_NED
    elif SN == '2026acd':
        type = 'Ia'
        z = 0.00758
        RA = 183.0531491 
        DEC = 13.221751
        discovery_mag = 18.876
        host_NED = 	0.100 #
        return type, z, RA, DEC, host_NED
    elif SN == '2026atb':
        type = 'Ia'
        z = 0.017941
        RA = 218.613723685
        DEC = 21.9496798525
        discovery_mag = 19.7674
        host_NED = 0.148 #SDSS J143428.51+215656.2
        return type, z, RA, DEC, host_NED
    elif SN == '2026gm':
        type = 'Ia'
        z = 0.0185
        RA = 66.1558857143 
        DEC = 33.8621121429
        discovery_mag = 18.08
        host_NED = 	1.078 #UGC 3028
        return type, z, RA, DEC, host_NED
    elif SN == '2026kc' or SN=='2026kc_subbed':
        type = 'Ia'
        z = 0.0228
        RA = 155.60448352
        DEC = 4.00043379084
        discovery_mag = 18.78
        host_NED = 0.093 #UGC05607
        return type, z, RA, DEC, host_NED
    elif SN =='2020ue':
        type = 'Ia'
        z = 0.003
        RA = 190.694917 
        DEC = 2.659481
        discovery_mag = 15
        host_NED = 0.078 #NGC 4636
        return type, z, RA, DEC, host_NED
    elif SN =='2026dix':
        type = 'IIb'
        z = 0.003185
        RA = 177.6559502
        DEC = 55.3535889
        discovert_mag = 16.2
        host_NED = 0.035
        return type, z, RA, DEC, host_NED
    
#next: add host_NED as an output, calculate E(V-B) for all the hosts by using
#A_lambda V = R_V*E(B-V) with Rv given in paper for the band we are using and 
# then use as prior in calculating H0 with SNooPy