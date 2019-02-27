###T(R) prof


import numpy as np
import const
from pylab import *

##Define the array containing the radius and angular information
size=1000
temp_inf=np.zeros((size,2,2))
#delta_r=(r_out-r_in)/size

##define the parameters
##default parameters
#delta_r=0.1
#r_out=1
#r_in=0.1
#inc=0
#D=1
#M=1*10**8
#M_dot=1
#R_star=3


##make the following code a definition
def flat_disk(r_out,r_in,inc,D,M,M_dot,R_star):


	##set this to 1 if you want to specify more parameters


	##text for plots
	delta_r=(r_out-r_in)/size
	text=np.zeros(6)
	text[0]=D
	text[1]=r_in
	text[2]=r_out
	text[3]=M_dot
	text[4]=M
	text[5]=inc

	##black hole scwarzschild radius
	def r_sch(M):
		rs=2*const.G*M/const.c**2
		return(rs)


	M=M*const.Msun
	M_dot=M_dot*const.Msun
	D=D*const.c/const.H0
	delta_r=delta_r*const.ps
	r_out=r_out*const.ps
	r_in=r_in*const.ps
	R_0=3.*r_sch(M)#*R_star
	#r_in=3.*r_sch(M)

	##populate the x and y dimensions of the array with radial and azimuthal coordinates 
	temp_inf[:,0,0]=1.*np.arange(size)/size * (r_out-r_in)  +r_in
	temp_inf[:,1,0]=1.*np.arange(size)/size * 2*np.pi
	#delta_r=(r_out-r_in)/size

	## solid angle  ## radii are stored in temp_inf[:,0,0]
	s_a=2*np.pi*temp_inf[:,0,0]*delta_r*np.cos(inc)/D**2



	##define the temperature radius profile for a steady state accretion disk in parsecs
	##viscous dissipation due gives rise to a factor of 3/8pi
	def T_R(M,M_dot,R,R_star):
		#M=M*const.Msun
		#M_dot=M_dot*const.Msun
		#R=R*const.ps
		#R_0=r_sch(M)*R_star
		T=(3*const.G*M*M_dot*(1-(R_0/R)**0.5)/(8*np.pi*const.sb*R**3))**0.25
		return(T)
	    


	

	###Now need to multiply by the solid angle element of each annuli


	T=T_R(M,M_dot,temp_inf[:,0,0],R_star)	

	##
	#define an array which stores wavelengths from the wein temperature of the lowest radi to highest radi

	wavmin=const.w_c/T[0]
	wavmax=const.w_c/T[size-1]
	wavelength=1.*np.arange(size)/size * (wavmax-wavmin) +wavmin   ## store the wavelength ranges over to calculate disk spectrum 	
	
	##plank spectrum as a function of wavelength have to sum over temperature contributions from all of disk
	B_T=np.zeros(size)
	for i in range(size):
		B_T[i]=(2*const.plank*const.c**2/wavelength[i]**5 * (np.exp(const.plank*const.c/(wavelength[i]*const.k*T)) - 1)**-1 * s_a).sum()


	
	
	## Just to plot everything

	output=0
	if output == 1:
		plot(wavelength,B_T)
		xlabel('Wavelength / m')
		ylabel('flux Wm^-1')
		xticks(fontsize=15)
		yticks(fontsize=15)
		title ('Z = '+str(text[0])+' .  r_in/r_out (pc)=' +str(text[1])+'/'+str(text[2])+'.  M='+str(text[4])+'.  M dot='+str(text[3])+'.  Inc='+str(text[5]))
		show()
	
	return(wavelength,B_T,text)
