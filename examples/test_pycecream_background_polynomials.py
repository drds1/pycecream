#!/usr/bin/env python
# coding: utf-8


import astropy_stark.myfake as mf
import matplotlib.pylab as plt
import numpy as np

output_directory = 'fit_synthetic_lightcurves_background_polynomial'

'''
same as test_pycecream.py but using a background polynomial light curve
'''


synthetic_data = mf.myfake(
    [4000.0,5000.0,5000.0,7000.0],
    [50.0,50.0,10.0,50.0],
    [1.0,1.0,2.0,1.0],
    thcent = 20.0
)

'''This recovers the synthetic data'''
dat = synthetic_data['echo light curves']


''' Now append a simple polynomial model'''
tlo = dat[0][:,0]
thi = dat[0][:,-1]
n = len(dat[0])
t = np.arange(n)/n*(thi - tlo) + tlo
x = 5.0 + 2.0*t + 0.1*t**2
x = (x - np.mean(x))/np.std(x)
sig = np.ones(n)*0.05
x = x + np.random.randn(n)*sig
background_data = np.zeros((n,3))
background_data[:,0] = t
background_data[:,1] = x
background_data[:,2] = sig
dat.append(background_data)

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


a.add_lc(dat[4],
         name = 'continuum (background polynomial)',
         kind='continuum',
         wavelength=8000.,
         background_polynomials = [0.1,0.1,0.1])
         #background_offset_start=[10.0,0.0],
         #vertical_scaling_start=[2.0,0.5])


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
