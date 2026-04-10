import sncosmo
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import constants as c

data = Table.read("sncosmo_input.dat", format="ascii")

#data = sncosmo.load_example_data()

#print(data)

# create a model
model = sncosmo.Model(source='salt2')
SN_list = ['2019np','2020ue','2025ahsa','2026acd','2026atb','2026gm','2026kc_subbed','2025acsn_subbed']

z_cmb = [0.0055,	0.0041,	0.0158	,0.0087,	0.0186,	0.0182,	0.0240,	0.0148]




for SN,z in zip(SN_list,z_cmb):
    #type, z, RA, DEC, host_NED = c.constants(SN)
    model.set(z)  # set the model's redshift.
    result, fitted_model = sncosmo.fit_lc(data, model,
                                        ['t0', 'x0', 'x1', 'c'])
    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
    print("Number of chi^2 function calls:", result.ncall)
    print("Number of degrees of freedom in fit:", result.ndof)
    print("chi^2 value at minimum:", result.chisq)
    print("model parameters:", result.param_names)
    print("best-fit values:", result.parameters)
    print("The result contains the following attributes:\n", result.keys())
    plt.show()

    x0 = result.parameters[2]
    x1 = result.parameters[3]
    c  = result.parameters[4]

    # SALT2 magnitude
    mB = -2.5 * np.log10(x0) + 10.635

    # nuisance parameters (typical literature values)
    alpha = 0.14
    beta  = 3.1
    M     = -19.3

    # distance modulus
    mu = mB + alpha*x1 - beta*c - M

    print("Distance modulus μ =", mu)



