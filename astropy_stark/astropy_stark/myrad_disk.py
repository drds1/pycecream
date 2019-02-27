### Script to plot an irradiated disk spectrum for a static disk

import numpy as np
import const
from pylab import *
from scipy import *


planck = 6.626307004e-34
c      = 2.99792458e8
boltz  = 1.38064852e-23
emsun  = 1.98855e30
gnewt  = 6.6740831e-11
sb     = 5.6700367e-8

####Define all the default paramters

#r_in = 3        ## ISCO radius in schwarzschild radii
#r_out = 1e5     ## desired outer radius of accretion disk 
#delta_r = 0.1  ## desired radial resolution 
#M = 1e8         ## Black hole mass (solar masses)
#D = 1 			## redshift distance to black hole
#inclination = 0 ## desired inclination (radians)
#M_dot = 1 		## Accretion rate (solar masses per year)

## Convert to SI units
#M = M*emsun
#rs = 3000*M/emsun
#r_in = r_in*rs
#r_out = r_out*rs
#delta_r = delta_r*rs
#M_dot=M_dot * const.Msun





##define the temperature radius profile for a steady state accretion disk in parsecs
##viscous dissipation due gives rise to a factor of 3/8pi
def T_R(M,M_dot,R,R_in):
	T=(3*gnewt*M*M_dot*(1-(R_in/R)**0.5)/(8*np.pi*sb*R**3))**0.25
	return(T)




###############################################################################################
##make the following code a definition
def flat_disk(r_out=1e5,r_in=1,inc=0,D=75.,embh=1.e7,emdot=1.0):



	####Define all the default paramters
	if default=='d':
		r_in = 3        ## ISCO radius in schwarzschild radii
		r_out = 1e5     ## desired outer radius of accretion disk 
		delta_r = 0.1  ## desired radial resolution 
		M = 1e8         ## Black hole mass (solar masses)
		D = 1 			## redshift distance to black hole
		inclination = 0 ## desired inclination (radians)
		M_dot = 1 		## Accretion rate (solar masses per year)
		                ## Total luminosity of point x ray source (see frank, king , rainne)

	## Convert to SI units
	M = embh*emsun
	rs = 3000*embh
	r_in = r_in*rs
	r_out = r_out*rs
	delta_r = delta_r*rs
	M_dot=emdot * emsun

	
	### Luminosity of the central source (powered by acrretion. May need to introduce efficiency
	##parameter later
	Lum = gnewt * M * M_dot / r_in
	
	##populate the x and y dimensions of the array with radial and azimuthal coordinates 
	rad_inf=1.*np.arange(np.ceil((r_out-r_in)/delta_r))*(r_out-r_in) + r_in
	

	

	### Obtain the temperature radius relation, wavelengths of the corresponding peak wavelengths
	## and solid angle element of each radii
	s_a=2*np.pi*rad_inf*delta_r*np.cos(inc)/D**2
	T=T_R(M,M_dot,rad_inf,r_in)	
	wavelength = const.w_c/T


	##plank spectrum as a function of wavelength have to sum over temperature contributions from all of disk
	B_T=np.zeros(len(rad_inf))
	for i in range(len(rad_inf)):
		B_T[i]=(2*const.plank*const.c**2/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T)) - 1)**-1 * s_a).sum()


	
	
	
	####Now define the emission associated by irradiation from point source
	## A thick disc has a depth (H) as a function of radius (x) given by 
	alph=9./7
	
	h_inf=rad_inf**9./7
	diff_h=9./7*rad_inf**2/7
	
	
	##solid angle element
	##see note book for derivation
	sa2=2*np.pi*rad_inf/D**2*(np.sin(inc)-np.cos(inc)**diff_H)
	
	##Also different temperature profile (see note book or frank, king and rainne chapt 5.10)
	
	T_r=(Lum/(4*np.pi*rad_inf**2*const.sb)*h_inf/rad_inf*(1.*alpha-1.)*(1-beta))**0.25
	
	## flux. As with plank spectrum
	for i in range(len(rad_inf)):
		B_T_r[i]=(2*const.plank*const.c**2/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T_r)) - 1)**-1 * sa2).sum()
	
	
	
	
	
	
	
	## Just to plot everything

	output=0
	if output == 1:
		plot(math.log(wavelength,10),math.log(B_T,10))
		xlabel('Wavelength / m')
		ylabel('flux Wm^-1')
		xticks(fontsize=15)
		yticks(fontsize=15)
		title ('Z = '+str(text[0])+' .  r_in/r_out (pc)=' +str(text[1])+'/'+str(text[2])+'.  M='+str(text[4])+'.  M dot='+str(text[3])+'.  Inc='+str(text[5]))
		show()
	
	return(wavelength,B_T,B_T_r)
