import numpy as np

#define function to convert sdss magnitudes to panstarrs magnitudes (from Morganson et al 2014)
#input gmag_sub_imag g magnitude minus i magnitude
#MUst have g and i band magnitudes in sdss

def diff_magps(gmag_sub_imag):
 #gmagsdss = magsdss[0]
 #imagsdss = magsdss[1]
 a0 = np.array([0.00128,-0.00518,0.00585,0.00144])
 a1 = np.array([-0.10699,-0.03561,-0.01287,0.07379])
 a2 = np.array([0.00392,0.02359,0.00707,-0.0336])
 a3 = np.array([0.00152,-0.00447,-0.00178,0.00765])
 gi = 1.*gmag_sub_imag#gmagsdss - imagsdss
 mpsout = a0 + a1*gi + a2*gi**2 + a3*gi**3
 return(mpsout)
 
 
#alternative approach by finkbeiner + 2015 now need difference in panstarrs g and i as input 
#this one has ugrizy filters
def diff_magps_fink(gmag_sub_imag):
 #gmagps = magps[0]
 #imagps = magps[1]
 a0 = np.array([0.04438,-0.01808,-0.01836,0.01170,-0.01062,0.08924])
 a1 = np.array([-2.26095,-0.13595,-0.03577,-0.00400,0.07529,-0.20878])
 a2 = np.array([-0.13387,0.001941,0.02612,0.00066,-0.03592,0.10360])
 a3 = np.array([0.27099,-0.00183,-0.00558,-0.00058,0.00890,-0.02441])
 gi = 1.*gmag_sub_imag#gmagsdss - imagsdss
 mpsout = a0 + a1*gi + a2*gi**2 + a3*gi**3
 return(mpsout)
 
 
 
 
 
 
def sdss_mag_2_flux(mag,sigmag,filter):
 
 if (filter == 'u'):
  b =1.4e-10
 elif (filter == 'g'):
  b = 0.9e-10
 elif (filter == 'r'):
  b = 1.2e-10
 elif (filter == 'i'):
  b = 1.8e-10
 elif (filter == 'z'):
  b = 7.4e-10
 
 b2 = 2*b
 a = (-np.log(10)/2.5*mag - np.log(b))#*2*b
 siga = np.abs(-np.log(10.)/2.5*sigmag)
  
 
 f_f0 = np.sinh(a)*b2
 sig_f_f0 = np.abs(np.cosh(a)*siga)*b2
 
 flux = 3631.*10**6 * f_f0
 sigflux = 3631*10**6 * sig_f_f0
 return(flux,sigflux)