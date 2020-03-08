# PyceCREAM

Here is a python implementaiton of my accretion disc and emission line lightcurve-fitting code (previously CREAM). This guide briefly covers generating synthetic data and calling a new pycecream object to ingest and fit the accretion disc model (or emission line model) to a set of input light curves. I also demonstrate how to access the output of the pycecream fit. The output includes the fitted light curves, any new light curve data points after merging, fitted response functions and parameter MCMC chain histories for the disc and/or tophat response parameters.

Most of these features are used in some form or another from a previous f90 version of this code (CREAM) in the following literature

* Grier et al in prep
* Grier et al 2018    https://iopscience.iop.org/article/10.3847/1538-4357/aa98dc/pdf
* Starkey et al 2017  https://ui.adsabs.harvard.edu/#abs/arXiv:1611.06051
* Starkey et al 2016  https://ui.adsabs.harvard.edu/#abs/arXiv:1511.06162

Please send questions to ds207@st-andrews.ac.uk. Though I am currently taking a break from academia and may take some time to respond, I will try to do so as soon as possible.


## Requirements & Installation

Please ensure that you have a fortran compiler installed. I use Gfortran. If you have an alternate (e.g ifort), please indicate the standard command used to call the fortran compiler using the ```fortran_caller``` argument (default is ```fortran_caller = gfortran```).


command These are fairly easy to install from macports or wget etc. Also a Python version is required (I am using 3.7 but even 2 should be fine). The it's just...

```
pip install pycecream
```


#  Section 1: Generate Synthetic Light Curves

In this example we generate 4 disk light curves and 2 emission-line light curves modelled as a top-hat with a 20-day lag. The code below generates a list where each index contains an Nx3 numpy array for each light curve. The 3 vertical axis for each light curve are the time, flux and noise respectively (query synthetic_data['echo lightcurves'][0] for an example of the format required when inputting your own light curve data).

The example below combines continuum and line light curves and illustrates a case in which you may have two of the same emission line (and so want to fit with the same response function model) but observed from different telescopes that require seperate noise models.


```python
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
```

#  Section 2: Settup and run PyceCREAM




```python
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
a.project_folder = 'fit_synthetic_lightcurves'



'''
Add each of the light curves in the simulation. 
In this case we are using the "dat" output from the synthetic data above.
'''
a.add_lc(dat[0], name = 'continuum 4000')
a.add_lc(dat[1], name = 'continuum 5000')
a.add_lc(dat[2], name = 'continuum 5000 (b)')
a.add_lc(dat[3], name = 'continuum 7000')

#If adding a line light curve, must indicate using the "kind" argument
a.add_lc(dat[4],name='test line 1',kind='line')

#If we want the same line response function model, set "share_previous_lag"=True
a.add_lc(dat[5],name='test line 1 (shared)',kind='line',share_previous_lag=True)



'''
specify the numnber of MCMC iterations. Normally at least several thousand are necessary but shorter numbers 
can be used just to check everything is working is done here.
'''
a.N_iterations=40

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
```

    pycecream path... /Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages/pycecream
    copying file...
    /Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages/pycecream
                       name  type  wavelength            noise model  \
    0        continuum 4000  line        -1.0  [var, multiplicative]   
    0        continuum 5000  line        -1.0  [var, multiplicative]   
    0    continuum 5000 (b)  line        -1.0  [var, multiplicative]   
    0        continuum 7000  line        -1.0  [var, multiplicative]   
    0           test line 1  line        -1.0  [var, multiplicative]   
    0  test line 1 (shared)  line        -1.0  [var, multiplicative]   
    
      share previous lag temporary file name      mean  standard deviation  \
    0              False          line_0.dat  3.800485            0.796293   
    0              False          line_1.dat  3.590073            0.675132   
    0              False          line_2.dat  3.593360            0.696788   
    0              False          line_3.dat  3.277164            0.524617   
    0              False          line_4.dat -0.000530            1.001517   
    0               True          line_5.dat  0.003286            1.025419   
    
       tophat centroid  tophat centroid step  tophat width  tophat width step  
    0              0.0                   5.0           2.0                0.0  
    0              0.0                   5.1           2.0                0.0  
    0              0.0                   5.2           2.0                0.0  
    0              0.0                   5.3           2.0                0.0  
    0              0.0                   5.4           2.0                0.0  
    0              0.0                   5.4           2.0                0.0  


# Examine the output

There are 2 output dataframes.

## 1) output_lightcurves = a.get_light_curve_fits():
This a dictionary of 3 data frames.

    1.1) output_lightcurves['model']: standard time, model, error envelope for each file

    1.2) output_lightcurves['merged model'] AS above but with the error bars, vertical and horrizontal scalings applied relative to the reference model. Not sure but I think the reference model defaults to the first occurence of a particular wavelength in the order that it was added in self.add_lc

    1.3) output_lightcurves['merged data'] DICTIONARY (since the input data light curves can be different sizes) The same transformations but applied to the input light curve data. useful if using cream only to merge the orriginal light curves from different telescopes to a new scale for further study elsewhere

## 2) output_chains = a.get_MCMC_chains(): 
These are the MCMC chains for each parameter.



```python
'''
Get the mcmc chains and output fits. 
Each of these arguments come with a "location" argument where you can point to a 
previous simulation and recover the outputs. 
If this is left blank we default to the current simulation
'''
output_chains = a.get_MCMC_chains(location = None)
output_lightcurves = a.get_light_curve_fits(location = None)


'''
NEW: 11/12/2019 Now the fourier chains are available as a pandas 
dataframe.
Stats on the sine and cosine parameters are also available for each 
freuqency accessible in the `fourier_stats` dictionary element of this
`get_MCMC_fourier_chains` function.
'''
output_fourier_chains = a.get_MCMC_fourier_chains(location=None)
fourier_chains = output_fourier_chains['fourier_chains']
fourier_stats = output_fourier_chains['fourier_stats']


'''
make figures of the fit, posterior, light curves etc. file prefix tells the code where you want to save the output.
The figure plotting is somewhat primitive and is a relic of when I still used cream. You may prefer to use your own
output figures with the output of the "get_MCMC_chains" and "get_light_curve_fits" functions above.
'''
a.plot_results(file_prefix='fit_figures')




'''
figures can also be made on an indivdual basis with axes objects returned from python plotting functions
'''
#plot the fitted light curves.
a.plot_lightcurves()
plt.show()


#plot the driving light curve
a.plot_driver()
plt.show()


#plot the parameter trace plots
a.plot_trace()
plt.show()


#plot the covariance parameter plot for the disc parameters
a.plot_posterior()
plt.show()



```

    cream_lcplot plotting results from... fit_synthetic_lightcurves/simulation_files
    -15.825983 [3.80048513 3.80048513 3.80048513 3.80048513 3.80048513] 0
    -15.825983 [3.59007335 3.59007335 3.59007335 3.59007335 3.59007335] 1
    -15.825983 [3.59336019 3.59336019 3.59336019 3.59336019 3.59336019] 2
    -15.825983 [3.27716422 3.27716422 3.27716422 3.27716422 3.27716422] 3
    making posterior plot.... posterior_fit_figures__1.pdf
    unable to make covariance plot for disc posteriors. Please check at least some of these are set to varyin the fit.
    fit_synthetic_lightcurves/simulation_files/output_20190406_001/G_plot.pdf
    Nth  6  Ndisk 1
    cream_lcplot plotting results from... fit_synthetic_lightcurves/simulation_files/output_20190406_001
    -15.825983 [3.80048513 3.80048513 3.80048513 3.80048513 3.80048513] 0
    -15.825983 [3.59007335 3.59007335 3.59007335 3.59007335 3.59007335] 1
    -15.825983 [3.59336019 3.59336019 3.59336019 3.59336019 3.59336019] 2
    -15.825983 [3.27716422 3.27716422 3.27716422 3.27716422 3.27716422] 3



![png](test_pycecream_tophat_0lag_files/test_pycecream_tophat_0lag_6_1.png)



![png](test_pycecream_tophat_0lag_files/test_pycecream_tophat_0lag_6_2.png)


    cream_lcplot plotting results from... fit_synthetic_lightcurves/simulation_files/output_20190406_001



![png](test_pycecream_tophat_0lag_files/test_pycecream_tophat_0lag_6_4.png)


    cream_lcplot plotting results from... fit_synthetic_lightcurves/simulation_files/output_20190406_001



![png](test_pycecream_tophat_0lag_files/test_pycecream_tophat_0lag_6_6.png)


    cream_lcplot plotting results from... fit_synthetic_lightcurves/simulation_files/output_20190406_001



```python
# how to install python 3 environment (skip the netcdf4 line) matplotlib should be ok now
# https://salishsea-meopar-docs.readthedocs.io/en/latest/work_env/python3_conda_environment.html

```
