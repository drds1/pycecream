import numpy as np
import pycecream
import astropy_stark.myfake as mf


'''
A test example on using the pycecream code. Section 1 deals with generating a set of synthetic data using
the myfake function from my astro_python library. This generates a list where each index contains an Nx3 numpy 
array for each light curve where each of the 3 vertical axis are the time, light curve and noise respectively
(query synthetic_data['echo lightcurves'][0] for an example of the format required when inputting
your own light curve data).

Skip to Section 2 if you allready have your light curve data as a numpy N x 3 array for each light curve 
and want to begin using cream immediately
'''


''' Section 1
Generate synthetic light curves (4 disk light curves and 2 top-hat lagged light curve with 20-day lag).

The top-hat light curves illustrate an example where you may have two of the same emission line 
observed from different telescopes, so want to share the same resaponse function but optimize the 
noise model seperately.

'''
synthetic_data = mf.myfake(
    [4000.0,5000.0,5000.0,7000.0,-1.0,-1.0],
    [50.0,50.0,10.0,50.0,50,10.],
    [1.0,1.0,2.0,1.0,1.0,3.0],
    thcent = 20.0
)
dat = synthetic_data['echo light curves']



''' Section 2
Now the synthetic data is generated. Input each into cream and set the number of iterations and mdot step size
'''
a = pycecream.pycecream()
a.output_directory = 'fit_synthetic_lightcurves'
a.add_lc(dat[0], name = 'continuum 4000')
a.add_lc(dat[1], name = 'continuum 5000')
a.add_lc(dat[2], name = 'continuum 5000 2')
a.add_lc(dat[3], name = 'continuum 7000')
a.add_lc(dat[4],name='test line 1',kind='line')
a.add_lc(dat[5],name='test line 1 (shared)',kind='line',share_previous_lag=True)
a.N_iterations=40
a.p_accretion_rate_step = 0.1

'''
list input settings and run
'''
print(a.lightcurve_input_params)
a.run()



'''
Get the mcmc chains
'''
output_chains = a.get_MCMC_chains(location = None)


'''
Get the output light curves.
This is a dictionary of 3 data frames

output_lightcurves['model'] standard time, model, error envelope for each file

output_lightcurves['merged model'] AS above but with the error bars, vertical and horrizontal scalings
applied relative to the reference model. Not sure but I think the reference model defaults to the first occurence
of a particular wavelength in the order that it was added in self.add_lc

output_lightcurves['merged data'] DICTIONARY (since the input data light curves can
be different sizes)
The same transformations but applied to the input light curve data.
useful if using cream only to merge the orriginal light curves from different telescopes
to a new scale for further study elsewhere
'''
output_lightcurves = a.get_light_curve_fits(location=None)
print('output light curves',output_lightcurves.keys())

# how to install python 3 environment (skip the netcdf4 line) matplotlib should be ok now
# https://salishsea-meopar-docs.readthedocs.io/en/latest/work_env/python3_conda_environment.html