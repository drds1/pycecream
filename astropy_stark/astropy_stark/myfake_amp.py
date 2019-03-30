#import mean flux in mjy
import numpy as np

fluxmean = 2.0

wav = 5000.





def abmag(fmjy):
 return(-2.5 * np.log10(fmjy/1000.) + 8.9)


def absmag(abm,dMpc):
 return( abm - 5*np.log10(dMpc) - 25. )

#uses empirical results obtained by Macleod, Ivezic et al 2010 
def drwtau_mc(embh,wav, absmag):
 A = 2.4
 B = 0.17
 C = 0.03
 D = 0.21
 taulog = A + B*np.log10(wav/4000.) + C*(absmag + 23.) + D*np.log10(embh/1.e9)
 return(10**taulog)




def drwsfinf_mc(embh, wav, absM, z = 0.):
 #! uses empirical results obtained by Macleod, Ivezic et al 2010 (DIFFERENT NUMBERS FROM drwtau_mc.f90)
 A = -0.56
 B = -0.479
 C = 0.111
 D = 0.11
 E = 0.07
 sfinflog = A + B*np.log10(wav/4000.) + C*(absM + 23.) + D*np.log10(embh/1.e9) + E*np.log10(1.+z)
 return(10**sfinflog)





def mfamp(embh,wav,fmjy,tlen,dMpc, z= 0.0):
 abm    = abmag(fmjy)
 absM   = absmag(abm,dMpc)
 taumc  = drwtau_mc(embh,wav,absM)
 sfinf  = drwsfinf_mc(embh, wav, absM, z = z)
 if taumc == 0:
  sf_inf = sfinf
 else:
  sf_inf = sfinf * np.sqrt( (1. - np.exp(-tlen/taumc)) )

 sdmag = np.abs(sf_inf)/np.sqrt(2)
 #this is the variance in ab magnitudes. Need to change to fluxes
 a = 1./1000/3631
 sdmjy = np.abs( -0.4 * 10**(-abm*0.4) / a * np.log(10.)  * sdmag )
 
  
 #print('wav, mean_AB, rms_AB ', wav, abm, sdmag,'  <-- report myfake_amp --> mean_mJy, rms_mJy ',fmjy, sdmjy)
 
 return(sdmjy)
 
 
 
 
 
 
 
 
 