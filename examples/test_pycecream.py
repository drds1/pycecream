#!/usr/bin/env python
# coding: utf-8


import astropy_stark.myfake as mf
import matplotlib.pylab as plt
import numpy as np

output_directory = 'fit_synthetic_lightcurves'

'''
mf.myfake arguments are

wavelengths: enter the wavelengths (-1 indicates an emission line light curve modelled with a top-hat response),

snr: set the signal-to-noise relative to light curve rms

cadence:set the mean cadence

top hat centroid: set the centroid for the top-hat (I think thats what this does but the line lag 
thing is still newish so Im used to just making continuum light curve)
'''


synthetic_data = mf.myfake(
    [4000.0,5000.0,5000.0,7000.0,-1.0,-1.0],
    [50.0,50.0,10.0,50.0,50,10.],
    [1.0,1.0,2.0,1.0,1.0,3.0],
    thcent = 20.0
)

'''This recovers the synthetic data'''
dat = synthetic_data['echo light curves']


# #  Section 2: Settup and run PyceCREAM
# 
# 

# In[2]:


import pycecream

#instantiate a pycecream object
a = pycecream.pycecream()

'''
If you use a fortran compiler other than gfortran please indicate here.
I just re-enter gfortran here for demonstration purposes even though 
this is unecassary as gfortran is the default argument.
'''
a.fortran_caller = 'gfortran'



'''Choose an output directory in which to save the results. 
This will be a new directory that you have not previously created (pycecream will make it automatically).

NOTE: Each new cream simulation must have a new name for "output_directory argument below 
otherwise an excpetion is raised. This is to prevent accidentally overwriting previous simulations. 
I might change this in a future version 
'''
a.project_folder = output_directory

#test the merging by adding offset to dat1
d1 = np.array(dat[1])
d1[:,1] = d1[:,1] - np.mean(d1[:,1]) + 232.
dat[1] = d1

'''
Add each of the light curves in the simulation. 
In this case we are using the "dat" output from the synthetic data above.
'''
a.add_lc(dat[0],
         kind='continuum',
         wavelength=4000.,
         name = 'continuum 4000')
         #background_offset_start=[10.0,0.0],
         #vertical_scaling_start=[2.0,0.5])
a.add_lc(dat[1],
         name = 'continuum 5000',
         kind='continuum',
         wavelength=5000.)
         #background_offset_start=[10.0,0.0],
         #vertical_scaling_start=[2.0,0.5])
a.add_lc(dat[2],
         name = 'continuum 5000 (b)',
         kind='continuum',
         wavelength = 5000.)
         #background_offset_start=[10.0,0.0],
         #vertical_scaling_start=[2.0,0.5])

a.add_lc(dat[3],
         name = 'continuum 7000',
         kind='continuum',
         wavelength=7000.)
         #background_offset_start=[10.0,0.0],
         #vertical_scaling_start=[2.0,0.5])

#If adding a line light curve, must indicate using the "kind" argument
a.add_lc(dat[4],name='test line 1',kind='line',
         background_offset_start=[10.0,0.0],
         extra_variance_prior = [0.1,1.0],
         multiplicative_errorbar_prior = [10.0,0.0000001],
         vertical_scaling_start=[2.0,0.5],
         vertical_scaling_prior=[0.0,0.1],
         background_offset_prior=[5.0,0.0001],
         tophat_width_prior=[0.0, -0.1],
         tophat_centroid_prior=[12.4, 0.00000001]
         )

#If we want the same line response function model, set "share_previous_lag"=True
a.add_lc(dat[5],name='test line 1 (shared)',kind='line',share_previous_lag=True,background_offset_start=[10.0,3.3],vertical_scaling_start=[2.0,0.5])



'''
specify the numnber of MCMC iterations. Normally at least several thousand are necessary but shorter numbers 
can be used just to check everything is working is done here.
'''
a.N_iterations=100

'''
specify the step sizes for the fit parameters. 
Here we are setting the accretion rate step size to vary by ~ 0.1 solar masses per year.
'''
a.p_accretion_rate_step = 0.1

'''
Check the input settings are ok prior to running
'''
print(a.lightcurve_input_params)

'''
RUN!
'''
a.run()
op = a.get_flux_flux_analysis(plotfile='fluxflux.pdf',xlim=[-4,4])
plt.show()


'''
get chains
'''
chains = a.get_MCMC_chains()
fourier_chains = a.get_MCMC_fourier_chains()
cols = list(chains.columns)
fcols = [c for c in cols if 'noise m ' in c]
fchains = chains[fcols]


'''
clean up output directory DONT DO FOR REAL SIMULATION 
AS THIS DELETES ALL RESULTS
'''
import os
os.system('rm -rf '+output_directory)
