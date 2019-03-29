#big fix nov 1st 2017 for very large black hole masses (schwarzchild radius bigger than 1 ld - reference radius and so 1 - rin/rld is negative. This screws up temperature radius law. Fix by setting temperature at all inner radii to zero and only using the scaling for radi above this

import numpy as np


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!define T0 in T = T0 (r/r0)**alpha#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def temp0(dotmmlog, emlog, sigdotmmlog =0, sigemlog = 0.0, alpha_in = 0.75, sig_alphain=0.0,eta=0.1, alb=0, hxrs = 3,r0ld = 1):


 gnewt = 6.673e-11
 c     = 2.9979e8
 sig   = 5.67e-8
 emsun   = 2.e30 
 ld    = 2.59e13
 em    = 10**emlog
 sigem = em*np.log(10)*sigemlog
 
 
 rs    = 2*gnewt*em*emsun/c**2

 hx    = hxrs*rs
 
 if (r0ld < 0):
  r0    = np.abs(r0ld*rs)
 else: 
  r0    = r0ld*ld


 emsunyr = emsun/(3600.*24*365) 
 
 dotmm = 10**dotmmlog
 sigdotmm  = dotmm*np.log(10)*sigdotmmlog

 emdot = dotmm/em
 sigemdot = emdot*np.sqrt( (sigdotmm/dotmm)**2 + (sigem/em)**2 )
  
 
 a         = 3*gnewt*emsun*emsunyr/(8*np.pi*sig)
 aextra    = r0**(4*alpha_in)
 sigaextra = aextra*np.log(r0)*4*sig_alphain
 anew     = a/aextra
 siganew  = anew / aextra * sigaextra
 

 
 b        = anew*dotmm
 sigb     = b*np.sqrt( (siganew/anew)**2 + (sigdotmm/dotmm)**2)   
 
 #for i in range(np.shape(b)[0]):
 # print i,a,b[i],sigb[i]
 
 #print 'diagnostics'
 #print hx,(1.-a),eta,emsunyr,c**2,sig,r0
 #raw_input()
 
 c         = hx*(1.-alb)*eta*emsunyr*c**2/(4*np.pi*sig)
 cextra    = 1.*aextra
 sigcextra = 1.*sigaextra 
 cnew      = c/cextra
 sigcnew  = cnew /cextra * sigcextra
 

 d        = cnew*emdot
 sigd     = d * np.sqrt( (sigcnew/cnew)**2 + (sigemdot/emdot)**2 )

 #for i in range(np.shape(b)[0]):
 # print i,c,d[i],sigd[i]
 #raw_input() 
  
 baddc     = b+d 
 sig_baddc = np.sqrt(sigb*sigb + sigd*sigd)  
 
 temp0     = baddc**0.25
 sigtemp0  = temp0/baddc*0.25*sig_baddc
 #print sigtemp0,temp0, sig_baddc, baddc,'dsfdsfs', sigb,b, sigd,d, 'dfs', sigcnew,cnew, sigemdot,emdot
 
 
 return(temp0,sigtemp0,b**0.25,d**0.25)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!








#function to return the radius for a specific temperature r = (T1/T)^[1/alpha] r1
#can include errors on T1 and alpha
def rsub(Tsub, r1, T1, alpha, sigT1= 0.0, sigalpha = 0.0):

 alpha_1     = 1./alpha
 alpha_1_sig = alpha_1 * sigalpha / alpha 
 
 T1_T = T1/Tsub
 
 T1_T_alpha_1     = T1_T**alpha_1
 
 T1_T_alpha_1_sig = T1_T_alpha_1 * np.sqrt( (alpha_1/T1_T * sigT1/Tsub)**2 + (np.log(T1_T)*alpha_1_sig)**2 )
 
 
 a = T1_T_alpha_1 * r1
 siga = T1_T_alpha_1_sig * r1
 
 return(a,siga)




def tv0(em,emdot,r0=1.0,rinsch=3.0):
 sb      = 5.670367e-8
 gnewt   = 6.67408e-11
 msun    = 1.98855e30
 secyr   = 31557600.0
 ld      = 2.59020683712e+13
 c       = 299792458.0
 rinld   = rinsch * 2*gnewt*msun*em/c/c/ld
 
 if (rinld > r0):
  r0in = np.ceil(rinld/r0)*r0
  print('mytemp0.py inner radius bigger than reference radius (Mbh too big)  changing reference radius')
 else:
  r0in = r0
 
 
 
 tv0out4 = em*emdot*3*gnewt*msun*msun/secyr/8/np.pi/sb/(r0*ld)**3 * (1. - np.sqrt(rinld/r0in))
 
 tv0out2 = np.sqrt(tv0out4)
 tv0out  = np.sqrt(tv0out2)
 return(tv0out)
 
 
 
def ti0(em,emdot,r0=1,hxsch=3,eta = 0.1):
 sb      = 5.670367e-8
 gnewt   = 6.67408e-11
 msun    = 1.98855e30
 secyr   = 31557600.0
 ld      = 2.59020683712e+13
 c       = 299792458.0
 hxld    = hxsch * 2*gnewt*msun*em/c/c/ld
 
 d2_ld    = hxld*hxld + r0*r0
 d1ld     = np.sqrt(d2_ld)
 d3ld     = d2_ld*d1ld
 ti0out4   = eta*msun/secyr*emdot*c*c/8/np.pi/sb * hxld/d3ld / ld/ ld
 ti0out2   = np.sqrt(ti0out4)
 ti0out    = np.sqrt(ti0out2)
 return(ti0out) 
 
 
 
 
 
#input grid r[....nr]] in light days, output t[r]
def tr(r,t0v,t0i,embh, r0=1,alpha_visc=-0.75,alpha_irad=-0.75,rinsch = 3):

 t0i2 = t0i*t0i
 t0i4 = t0i2*t0i2
 
 t0v2 = t0v*t0v
 t0v4 = t0v2*t0v2
 
 av4 = alpha_visc*4
 
 ai4 = alpha_irad*4
 
 rsch = 1.15821e-10*embh 
 rinld = rinsch*rsch
  
  
 if (rinld > r0):
  r0in = np.ceil(rinld/r0)*r0
  print('mytemp0.py inner radius bigger than reference radius (Mbh too big)  changing reference radius')
 else:
  r0in = r0
  
 #print r.shape
 #print np.mean(r)
 #print t0v4*(r)**av4
 
 tv4 = t0v4*(r/r0)**av4 * (1 - np.sqrt(rinld/r)) / (1 - np.sqrt(rinld/r0in))
 ti4 = t0i4*(r/r0)**ai4
 
 tout2 = np.sqrt(tv4 + ti4)
 tout  = np.sqrt(tout2)
 return(tout)

 
 




