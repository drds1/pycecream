### Script to plot an irradiated disk spectrum for a static disk

import numpy as np
import cgs
from pylab import *
from scipy import *
from plank import *


####Define all the default paramters

r_in = 3        ## ISCO radius in schwarzschild radii
r_out = 1e5     ## desired outer radius of accretion disk 
delta_r = 10  ## desired radial resolution 
M = 1e8         ## Black hole mass (solar masses)
D = 1 			## redshift distance to black hole
inc = 0 ## desired inclination (radians)
M_dot = 1 		## Accretion rate (solar masses per year)
beta = 0
up_wav  = 1e-6
low_wav = 1e-9
delta_wav = 1e-9

## Convert to SI units
M = M*const.Msun
rs = const.r_sch(M)
r_in = r_in*rs
r_out = r_out*rs
delta_r = delta_r*rs
M_dot=M_dot * const.Msun
D=D*const.c/const.H0




##define the temperature radius profile for a steady state accretion disk in parsecs
##viscous dissipation due gives rise to a factor of 3/8pi
def T_R(M,M_dot,R,R_in):
	T=(3*const.G*M*M_dot*(1-(R_in/R)**0.5)/(8*np.pi*const.sb*R**3))**0.25
	return(T)

	
## flat disk spectrum irradiated
def T_D(A,lum,H0,x):
	T_D = (1-A)*lum/(4*np.pi*const.sb*(x**2+H**2)) * x/(x**2 + H**2)**0.5
	return(T_D)

	
## plank spectrum (frequency) but in terms of wavelength. given wavelgnth, and temperature

def btu(temp,wave):
	plank=2*const.plank*const.c/wave**3 * (np.exp(const.plank*const.c/(wave*const.k*temp)) - 1)**-1
	return(plank)
	


### solid angle flat disk

def sa_flat(r,D,inc,delta_r):
	s_a=2*np.pi*r*delta_r*np.cos(inc)/D**2	
	return(s_a)
	


### Luminosity of the central source (powered by acrretion. May need to introduce efficiency
##parameter later
Lum = const.G * M * M_dot / r_in

##populate the x and y dimensions of the array with radial and azimuthal coordinates 

rad_inf=1.*np.arange(np.ceil((r_out-r_in)/delta_r))/np.ceil((r_out-r_in)/delta_r)*(r_out-r_in) + r_in


	

### Obtain the temperature radius relation, wavelengths of the corresponding peak wavelengths
## and solid angle element of each radii
#s_a=2*np.pi*rad_inf*delta_r*np.cos(inc)/D**2
#T=T_R(M,M_dot,rad_inf,r_in)	
#wavelength = const.w_c/T
wavelength = 1.*np.arange(np.ceil((up_wav-low_wav)/delta_wav))/np.ceil((up_wav-low_wav)/delta_wav)*(up_wav-low_wav) + low_wav





