### Script to plot an irradiated disk spectrum for a static disk

import numpy as np
import const
from pylab import *
from scipy import *



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
wavelength = 1.*np.arange(np.ceil((up_wav-low_wav)/delta_wav))/np.ceil((up_wav-low_wav)/delta_wav)*(up_wav-low_wav) + low_wav

##plank spectrum as a function of wavelength have to sum over temperature contributions from all of disk
B_T=np.zeros(len(wavelength))
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

x=np.arange(len(rad_inf))
#y=np.arange(len(T))

for z in range(len(wavelength)):
	B_T[z]=(2*const.plank*const.c/wavelength[z]**3 * (np.exp(const.plank*const.c/(wavelength[z]*const.k*T[:])) - 1)**-1 * s_a[:]).sum()

#d=np.zeros(10)
#d[z]=
#show()





	
### Define BTR

B_T_r = np.zeros(len(wavelength))	
	
####Now define the emission associated by irradiation from point source
## A thick disc has a depth (H) as a function of radius (x) given by 
alpha=0
h0=3*rs
h_inf=h0*(rad_inf-r_in)**(alpha)
diff_h=alpha*h0*(rad_inf-r_in)**(alpha-1)
	
	
##solid angle element
##see note book for derivation
sa2=2*np.pi*rad_inf*delta_r/D**2*(np.cos(inc)-np.sin(inc)*diff_h)
#s_a=2*np.pi*rad_inf*delta_r*np.cos(inc)/D**2
##Also different temperature profile (see note book or frank, king and rainne chapt 5.10)
	
c_theta=1/diff_h /(1+1/diff_h**2)**0.5
n_check=np.where(np.isnan(c_theta) ==True)[0].astype(int)   # if you have a flat disk, diff = 0 which causes the cos_theta term to return a nan. This removes
##the effect by recognising that a flat disk would in fact return a cosine of 1

c_theta[n_check]=1
c_phi = rad_inf / (rad_inf**2 + (h0 - h_inf)**2)**0.5
#T_r=(Lum/(4*np.pi*rad_inf**2+*const.sb)*h_inf/rad_inf*(1.*alpha-1.)*(1-beta))**0.25
T_r=(Lum/(4*np.pi*rad_inf**2*const.sb)*c_theta*c_phi*(diff_h+(h0-diff_h)/rad_inf)*(1-beta))**0.25
## flux. As with plank spectrum
#for i in range(len(wavelength)):
#	if i%100 != 0:
#		print ' done '+str(i)+ ' of ' + str(len(rad_inf))
#	B_T_r[i]=(2*const.plank*const.c/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T_r)) - 1)**-1 * sa2).sum()
	
	
	
T_test= ((1-beta)*Lum/(4*np.pi*const.sb*(rad_inf**2 + r_in**2))*(rad_inf/(rad_inf**2 + r_in**2)**0.5))**0.25

for i in range(len(wavelength)):
	if i%100 != 0:
		print ' done '+str(i)+ ' of ' + str(len(rad_inf))
	B_T_r[i]=(2*const.plank*const.c/wavelength[i]**3 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T_test)) - 1)**-1 * sa2).sum()
	