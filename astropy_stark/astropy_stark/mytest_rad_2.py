### Script to plot an irradiated disk spectrum for a static disk

import numpy as np
import const
from pylab import *
from scipy import *
import scipy.integrate


####Define all the default paramters

r_in = 3        ## ISCO radius in schwarzschild radii
r_out = 1e4     ## desired outer radius of accretion disk 
delta_r = 1     ## desired radial resolution 
M = 1e8         ## Black hole mass (solar masses)
D = 1 			## redshift distance to black hole
inc = 0         ## desired inclination (radians)
M_dot = 1 		## Accretion rate (solar masses per year)
beta = 0

## Convert to SI units
M = M*const.Msun
rs = const.r_sch(M)
r_in = r_in*rs
r_out = r_out*rs
delta_r = delta_r*rs
M_dot=M_dot * const.Msun





##define the temperature radius profile for a steady state accretion disk in parsecs
##viscous dissipation due gives rise to a factor of 3/8pi
def T_R(M,M_dot,R,R_in):
	T=(3*const.G*M*M_dot*(1-(R_in/R)**0.5)/(8*np.pi*const.sb*R**3))**0.25
	return(T)



### Luminosity of the central source (powered by acrretion. May need to introduce efficiency
##parameter later
Lum = const.G * M * M_dot / r_in

##populate the x and y dimensions of the array with radial and azimuthal coordinates 
rad_inf=1.*np.arange(np.ceil((r_out-r_in)/delta_r))/np.ceil((r_out-r_in)/delta_r)*(r_out-r_in) + r_in


	

### Obtain the temperature radius relation, wavelengths of the corresponding peak wavelengths
## and solid angle element of each radii
s_a=2*np.pi*rad_inf*delta_r*np.cos(inc)/D**2
T=T_R(M,M_dot,rad_inf,r_in)	
wavelength = const.w_c/T


##plank spectrum as a function of wavelength have to sum over temperature contributions from all of disk
#B_T=np.zeros((len(rad_inf),len(wavelength)))
#B_T=np.zeros(10)
#B_T_2=np.zeros(10)

#sum_ar=np.arange(len(rad_inf))
#sum_ar=np.arange(10)
#for i in range(10):
#	B_T[i]=(2*const.plank*const.c**2/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T)) - 1)**-1 * s_a).sum()

#	if i%100 != 0:
#		print ' done '+str(i)+ ' of ' + str(len(rad_inf))


#test_ar=np.zeros((len(rad_inf),len(T),10))
#test_ar[:
#B_T_2[:]=
#(2*const.plank*const.c**2/wavelength[sum_ar[:]]**3

#B_T_2[sum_ar[:]]=(2*const.plank*const.c**2/wavelength[sum_ar[:]]**3 * (np.exp(const.plank*const.c/(wavelength[sum_ar[:]]*const.k*T)) - 1)**-1 * s_a).sum()

#x=np.arange(len(rad_inf))
#y=np.arange(len(T))

#for z in range(len(wavelength)):
def plank(R,wav):
	#wav=1e-7
	flux=2*const.plank*const.c/wav**3 * (np.exp(const.plank*const.c/(wav*const.k*(3*const.G*M*M_dot*(1-(r_in/R)**0.5)/(8*np.pi*const.sb*R**3))**0.25)) - 1)**-1 * 2*np.pi*R*np.cos(inc)/D**2
	return(flux)

wav=1e-7


flux= lambda R: 2*const.plank*const.c/wavelength[:]**3 * (np.exp(const.plank*const.c/(wavelength[:]*const.k*(3*const.G*M*M_dot*(1-(r_in/R)**0.5)/(8*np.pi*const.sb*R**3))**0.25)) - 1)**-1 * 2*np.pi*R*np.cos(inc)/D**2
#output=scipy.integrate.quad(flux,rad_inf[1],rad_inf[-1])
	
	#d=np.zeros(10)
#d[z]=
#show()
wav=1e-7

def plank(R,wav):
	test=(2*const.plank*const.c**2/wav**3 * (np.exp(const.plank*const.c/(wav*const.k*(3*const.G*M*M_dot*(1-(r_in/rad_inf[:])**0.5)/(8*np.pi*const.sb*rad_inf[:]**3))**0.25)) - 1)**-1 * 2*np.pi*rad_inf[:]*delta_r*np.cos(inc)/D**2).sum()
	return(flux)

result=plank(rad_inf[:],wavelength[:])



	
### Define BTR

B_T_r = np.zeros((len(rad_inf), len(wavelength)))	
	
####Now define the emission associated by irradiation from point source
## A thick disc has a depth (H) as a function of radius (x) given by 
alpha=9./7
	
h_inf=rad_inf**9./7
diff_h=9./7*rad_inf**2/7
	
	
##solid angle element
##see note book for derivation
sa2=2*np.pi*rad_inf/D**2*(np.sin(inc)-np.cos(inc)**diff_h)
	
##Also different temperature profile (see note book or frank, king and rainne chapt 5.10)


T_r=(Lum/(4*np.pi*rad_inf**2*const.sb)*h_inf/rad_inf*(1.*alpha-1.)*(1-beta))**0.25

## flux. As with plank spectrum
#for i in range(len(rad_inf)):
#	if i%100 != 0:
#		print ' done '+str(i)+ ' of ' + str(len(rad_inf))
#	B_T_r[i]=(2*const.plank*const.c**2/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T_r)) - 1)**-1 * sa2).sum()
	
	
	