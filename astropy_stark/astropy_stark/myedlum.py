## script to work out the eddington luminosity for an AGN with an input Mdot

import numpy as np

pi=np.pi
G=6.67384e-11
msun=1.9891e30
mp=1.67262e-27
c=2.9979e8
thomson_e=6.6524e-29
year=24.0*3600*365
watt2erg=1.0e7

## input umdot: accretion rate (M0 /yr)
##       um   : BH mass (M0)
## output ledd (Eddington luminosity in ergs s_1)
##        eddrat(Eddington ratio l/ledd)


def edd(um,umdot,eta):
#    ledd=4*pi*G*um*msun*mp/thomson_e *watt2erg * c
    ledd2=1.26e31*um*watt2erg
    l=eta*umdot/year*msun*c**2 *watt2erg
    
    eddrat= l/ledd2 
    
    return(l,ledd2,eddrat)
    

## input (IN LOG10) ummdot:mass times accretion rate (M0**2 /yr) with error ummdotsig
##       um   : BH mass (M0)  with error umsig
## if umlog or umsiglog < 0 then code assumes input is in linear NOT log units
## if ummdotversion = 1, code assumes ummdotlog and ummdotsiglog are in log units,
##... if -1 then ummdotlog and ummdotsiglog in linear units
##... if +2 the ummdotlog is actually umdotlog and ummdotsiglog is umdotsiglog
##... if -2 then ummdotlog is actually umdot and umdotsig
## output ledd (Eddington luminosity in ergs s_1)
##        eddrat(Eddington ratio l/ledd)
##        uncertainty in eddington ratio

def edd2(umlog,umsiglog,ummdotlog,ummdotsiglog,eta,ummdotversion=1):
#    ledd=4*pi*G*um*msun*mp/thomson_e *watt2erg * c
    
     

    
    if (umlog > 0):
     um = 10**umlog*(1. + 0.5*np.log(10)**2 * umsiglog*umsiglog)
    else:
     um = np.abs(umlog)
      
    if (umsiglog > 0):
     umsig = um*np.log(10.)*umsiglog
    else:
     umsig = 1.*umsiglog
     
    
    if (ummdotversion == 1):
     ummdot = 10**ummdotlog*(1. + 0.5*np.log(10)**2 * ummdotsiglog*ummdotsiglog)
     ummdotsig = ummdot*np.log(10.)*ummdotsiglog
     umdot=ummdot/um
     sigumdot=umdot*np.sqrt( (ummdotsig/ummdot)**2 + (umsig/um)**2 )
     
    elif (ummdotversion == -1):
     ummdot = 1.*ummdotlog
     ummdotsig = 1.*ummdotsiglog
     umdot=ummdot/um
     sigumdot=umdot*np.sqrt( (ummdotsig/ummdot)**2 + (umsig/um)**2 )
    elif (ummdotversion == 2):
     umdot = 10**ummdotlog*(1. + 0.5*np.log(10)**2 * ummdotsiglog*ummdotsiglog)
     sigumdot = umdot*np.log(10.)*ummdotsiglog
    elif (ummdotversion == -2):
     umdot = ummdotlog
     sigumdot = ummdotsiglog


    
    
    ledd2=1.26e31*um*watt2erg
    l=eta*umdot/year*msun*c**2 *watt2erg
    
    sigl=l/umdot*sigumdot
    sigledd2=ledd2/um*umsig
    
    
    
    eddrat= l/ledd2 
    sigeddrat=eddrat* np.sqrt( (sigl/l)**2 + (sigledd2/ledd2)**2 )
    
    
    
    umdotlog = ummdotlog - umlog
    umdotsig_log = np.sqrt(umsiglog*umsiglog + ummdotsiglog*ummdotsiglog)
    return(l,ledd2,eddrat,sigeddrat,umdot,sigumdot,umdotlog,umdotsig_log)
    
    
      
    
# input Mass, efficiency 
# output Eddington ratio
def eddrate(um,eta):
    eddmdot = 1.26e31 * um /eta / (c*c) # in Watts
    
    eddumdot = eddmdot * year / msun
    
    return(eddumdot)


#routine to calculate the mdot for a given M and eddington ratio !!
#inp um (black hole mass in solar masses)
#er eddington ratio (0-1)
#eta black hole accretion efficiency typically 0.1
#op  umdot solar masses per year

def ermin_mdotout(um,er,eta=0.1):
 umdot = er/460.7e6 / eta * um
 return(umdot)
