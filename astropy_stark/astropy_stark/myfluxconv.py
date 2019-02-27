import numpy as np




# go from difference of magnitudes to ratio of fluxes
# default option converts to mjy
def dmag2ratf(dmag,multiply=3631000.0):
    ratf=multiply*10**(-0.4*dmag)
    return(ratf)




## convert mags to fluxes (mJy)
def mabtonu(mab,sigmab = []):
    a = 3631000.0
    fnu=a*10**(-0.4*mab)
    if (sigmab == []):
     return(fnu)
    else:
     sigfnu = fnu*0.4*np.log(10.)*sigmab
     return(fnu,sigfnu)

#what if flux in mjy and error in mags  
def sigmab2sigfnu(fnumjy,sigmab):
 return(fnumjy*0.4*np.log(10)*sigmab)


# convert flam to fnu in jansky, flamin erg cm^-2 s^-1 A-1, wav in angsrom
def flam2fnu(flam,wavang):
 return(3.34e4 * wavang*wavang * flam)
 

#convert fnu to flam units as above jansky - erg cm^-2 s^-1 A-1
def fnu2flam(fnu,wavang):
 return(2.99401e-5 /(wavang*wavang) * fnu)
 

#convert erg s-1 cm -2 A-1 into watts/m2/hz
def flam2si(flam,wavang):
 a = 2.9792458e14
 return(flam/a/(wavang*wavang))