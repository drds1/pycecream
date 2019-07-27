#!/usr/bin/env python
# coding: utf-8

# # PyceCREAM
# 
# Here is a python implementaiton of my accretion disc and emission line lightcurve-fitting code (previously CREAM). This guide briefly covers generating synthetic data and calling a new pycecream object to ingest and fit the accretion disc model (or emission line model) to a set of input light curves. I also demonstrate how to access the output of the pycecream fit. The output includes the fitted light curves, any new light curve data points after merging, fitted response functions and parameter MCMC chain histories for the disc and/or tophat response parameters.
# 
# Most of these features are used in some form or another from a previous f90 version of this code (CREAM) in the following literature
# 
# * Grier et al in prep
# * Grier et al 2018    https://iopscience.iop.org/article/10.3847/1538-4357/aa98dc/pdf
# * Starkey et al 2017  https://ui.adsabs.harvard.edu/#abs/arXiv:1611.06051
# * Starkey et al 2016  https://ui.adsabs.harvard.edu/#abs/arXiv:1511.06162
# 
# Please send questions to ds207@st-andrews.ac.uk. Though I am currently taking a break from academia and may take some time to respond, I will try to do so as soon as possible.
# 
# 
# ## Requirements & Installation
# 
# Please ensure that you have a fortran compiler installed. I use Gfortran. If you have an alternate (e.g ifort), please indicate the standard command used to call the fortran compiler using the ```fortran_caller``` argument (default is ```fortran_caller = gfortran```).
# 
# 
# command These are fairly easy to install from macports or wget etc. Also a Python version is required (I am using 3.7 but even 2 should be fine). The it's just...
# 
# ```
# pip install pycecream
# ```
# 

# #  Section 1: Generate Synthetic Light Curves
# 
# In this example we generate 4 disk light curves and 2 emission-line light curves modelled as a top-hat with a 20-day lag. The code below generates a list where each index contains an Nx3 numpy array for each light curve. The 3 vertical axis for each light curve are the time, flux and noise respectively (query synthetic_data['echo lightcurves'][0] for an example of the format required when inputting your own light curve data).
# 
# The example below combines continuum and line light curves and illustrates a case in which you may have two of the same emission line (and so want to fit with the same response function model) but observed from different telescopes that require seperate noise models.

# In[1]:


import astropy_stark.myfake as mf
import matplotlib.pylab as plt

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
a.output_directory = 'fit_synthetic_lightcurves'



'''
Add each of the light curves in the simulation. 
In this case we are using the "dat" output from the synthetic data above.
'''
a.add_lc(dat[0],
         kind='continuum',
         wavelength=4000.,
         name = 'continuum 4000',
         background_offset_start=[10.0,0.0],
         vertical_scaling_start=[2.0,0.5])
a.add_lc(dat[1],
         name = 'continuum 5000',
         kind='continuum',
         wavelength=5000.,
         background_offset_start=[10.0,0.0],
         vertical_scaling_start=[2.0,0.5])
a.add_lc(dat[2],
         name = 'continuum 5000 (b)',
         kind='continuum',
         wavelength = 5000.,
         background_offset_start=[10.0,0.0],
         vertical_scaling_start=[2.0,0.5])

a.add_lc(dat[3],
         name = 'continuum 7000',
         kind='continuum',
         wavelength=7000.,
         background_offset_start=[10.0,0.0],
         vertical_scaling_start=[2.0,0.5])

#If adding a line light curve, must indicate using the "kind" argument
a.add_lc(dat[4],name='test line 1',kind='line',background_offset_start=[10.0,0.0],vertical_scaling_start=[2.0,0.5],
         vertical_scaling_prior=[0.0,0.1],
         background_offset_prior=[5.0,0.0001]
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

#
## # Examine the output
##
## There are 2 output dataframes.
##
## ## 1) output_lightcurves = a.get_light_curve_fits():
## This a dictionary of 3 data frames.
##
##     1.1) output_lightcurves['model']: standard time, model, error envelope for each file
##
##     1.2) output_lightcurves['merged model'] AS above but with the error bars, vertical and horrizontal scalings applied relative to the reference model. Not sure but I think the reference model defaults to the first occurence of a particular wavelength in the order that it was added in self.add_lc
##
##     1.3) output_lightcurves['merged data'] DICTIONARY (since the input data light curves can be different sizes) The same transformations but applied to the input light curve data. useful if using cream only to merge the orriginal light curves from different telescopes to a new scale for further study elsewhere
##
## ## 2) output_chains = a.get_MCMC_chains():
## These are the MCMC chains for each parameter.
##
#
## In[3]:
#
#
#'''
#Get the mcmc chains and output fits.
#Each of these arguments come with a "location" argument where you can point to a
#previous simulation and recover the outputs.
#If this is left blank we default to the current simulation
#'''
#output_chains = a.get_MCMC_chains(location = None)
#output_lightcurves = a.get_light_curve_fits(location = None)
#
#'''
#make figures of the fit, posterior, light curves etc. file prefix tells the code where you want to save the output.
#The figure plotting is somewhat primitive and is a relic of when I still used cream. You may prefer to use your own
#output figures with the output of the "get_MCMC_chains" and "get_light_curve_fits" functions above.
#'''
#a.plot_results(file_prefix='fit_figures')
#
#
#
#
#'''
#figures can also be made on an indivdual basis with axes objects returned from python plotting functions
#'''
##plot the fitted light curves.
#a.plot_lightcurves()
#plt.show()
#
#
##plot the driving light curve
#a.plot_driver()
#plt.show()
#
#
##plot the parameter trace plots
#a.plot_trace()
#plt.show()
#
#
##plot the covariance parameter plot for the disc parameters
#a.plot_posterior()
#plt.show()
#
#
#
## In[4]:
#
#
## how to install python 3 environment (skip the netcdf4 line) matplotlib should be ok now
## https://salishsea-meopar-docs.readthedocs.io/en/latest/work_env/python3_conda_environment.html
#
#