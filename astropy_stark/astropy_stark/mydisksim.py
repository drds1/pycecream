#calculates the black body spectrum of an accretion disc
#by default at 70 Mpc in mJy (if 1 else if mjy = 0 use erg/s/ang/cm^2)

#INPUT wav[...nwav] list of wavelengths to evaluate spectrum
#..... embh,emdot (the black hole mass and acretion rate in solar masses and solarmasses per year)
#..... deginc the inclination of the disc, 0 = face-on
#..... dl (optional) the luminosity distance to the disc in Mpc (default is 70)
#..... radlosch, radhisch, ngrid control the resolution of the disc grid (only change if having problems)
#..... mjy = 1 if want spectrum in mjy else get in erg/s/ang/cm^2
#
#OUTPUT mds[...nwav] spectrum
  
import numpy as np
import matplotlib.pylab as plt
import astropy_stark.mytemp0 as mt0
import astropy_stark.myplank as mp
twopi = 2*np.pi
deg2rad = twopi/360.




def mds(wav,embh,emdot,degi, dl = 70., radlosch=3.0,radhisch=10000., ngrid = 1000, mjy = 1):
 r0 = 1.0
 sb      = 5.670367e-8
 gnewt   = 6.67408e-11
 msun    = 1.98855e30
 secyr   = 31557600.0
 ld      = 2.59020683712e+13
 c       = 299792458.0
 rinld   = radlosch * 2*gnewt*msun*embh/c/c/ld
 if (rinld > r0):
  r0 = np.ceil(rinld/r0)
  #print('mydisksim.py: black hole mass too big for inner radius of 1 light day... expanding to',r0)

 cosi  = np.cos(deg2rad*degi)
 ldMpc = 1191286169.529
 ldMpc2 = ldMpc*ldMpc
 rsch   = 1.15821e-10*embh #
 radlo  = radlosch*rsch
 radhi  = radhisch*rsch
 rgrid  = np.logspace(np.log10(radlo),np.log10(radhi),ngrid)
 drgrid = rgrid[1:] - rgrid[:-1]
 drgrid = np.append(drgrid,drgrid[-1])
 #calculate appropriate temperatures for this embh and emdot
 t0v = mt0.tv0(embh,emdot,r0=r0,rinsch=radlosch)
 t0i = mt0.ti0(embh,emdot,r0=r0,hxsch=3,eta = 0.1)
 tr  = mt0.tr(rgrid,t0v,t0i,embh, r0=r0,alpha_visc=-0.75,alpha_irad=-0.75,rinsch = radlosch)
 
 
 #calculate the planck function for this radial grid
 fnu = []
 for wavnow in wav:
  
  if (wavnow == 0):
   fnu.append(0)
  else:
   bnu = mp.bnuvec(wavnow,tr,mjy=mjy)
   weight = rgrid * drgrid * bnu
   fnu.append( np.sum(weight)*twopi/dl/dl/ldMpc2 * cosi )
 
 return(fnu) 
 
 
 
 
#test the code
#
#dl = 70.0
#embh = 3.e9
#emdot = 100.
#wav = np.logspace(2.,6,1000)
#nrgrid = 1000
#degi = 0.0
#radlosch = 3.0
#radhisch = 10000.0 
#
#
#fnu = mds(wav,embh,emdot,degi, dl = 70., radlosch=3.0,radhisch=10000., ngrid = 1000, mjy = 1) 
##test diag
##for i in range(nrgrid):
## print rgrid[i], drgrid[i], bnu[i] 
# 
##test plot
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.plot(wav, fnu)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#plt.savefig('mydisksim_plot.pdf')
#
#
#
##test tr
##fig = plt.figure()
##ax1 = fig.add_subplot(111)
##ax1.plot(rgrid, tr)
##ax1.set_xscale('log')
##ax1.set_yscale('log')
##plt.savefig('mydisksimtr_plot.pdf')
##




 

