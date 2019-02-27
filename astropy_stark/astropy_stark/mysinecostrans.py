import numpy as np
from mylcgen import *
from myrandom import *
import matplotlib.pylab as plt



#function to fit the sine or cosine transformastion to an input light curve
#parameters: 
#time,y,sig the data arrays
# freq if 1 then use flo and fhi to obtain a suitable array of frequencies to calc the transform
# flo,fhi if -ve (and freq == 1) then go to abs(fhi) x the nyquist frequency for high frequency and abs(flo) times the repeat period for low frequency, df = low frequency
# dtres if -ve then go to abs(dtres) x the average time resolution
# costrans if 1 calcuate the cosine transform else (if 0) do the sine transformation

def sct(time,y,sig,freq=np.zeros(0),flo=-0.5,fhi=-2,dtres=-10,costrans=1,ploton = 1):
 
 #preamble !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 twopi = 2*np.pi
 tlen  = time[-1] - time[0]
 dt    = np.median(time[1:] - time[:-1])
 ndat  = np.shape(y)[0]
 tlo   = np.min(time)
 thi   = np.max(time)
 
 if (dtres < 0):
  dtresplot = dt/np.abs(dtres)
 else: 
  dtresplot = dtres
  
 #generate angular frequencies
 if (np.shape(freq)[0] == 0):
  wlo = twopi*flo
  whi = twopi*fhi
  if (flo <0):
   wlo = twopi/tlen*np.abs(flo)
  if (fhi < 0):
   whi = twopi/dt * np.abs(fhi)
  w = np.arange(wlo,whi+wlo,wlo)
 else:
  w = twopi*freq
 
 #end of preamble!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (costrans == 1):
  phase = 0
 elif (costrans == 0):
  phase = -np.pi/2
 
 
 #now compute the hessian matrix
 nw = np.shape(w)[0]
 np2 = 2*nw
 
 cvec = np.ones(nw)
 hes = np.ones((nw,nw))
 sig2 = sig*sig
 #save the sine and cosine value in  2d array
 cwt = np.ones((ndat,nw))
 swt = np.ones((ndat,nw))
 for iw in range(nw):
  cwt[:,iw] = np.cos(w[iw]*time +phase)
 
 for ikp in range(nw):
  wp = w[ikp]
  
  wp2 = wp*wp
  cvec[ikp]   = np.sum(y*cwt[:,ikp]/sig2)
  #print ikp,nw
  for ik in range(nw):
   #cosine column
   hnow_c = np.sum(cwt[:,ik]*cwt[:,ikp]/sig2)
   if (ik == ikp):
    hnow_c = hnow_c + wp2
   hes[ikp,ik]   = hnow_c
 
    
 
       
 #invert the matrix to get the covariance matrix
 cov = np.linalg.inv(hes)
 
 
 #find the parameters
 parm  = cov.dot(cvec)
 ckout = parm
 
 
 
 
 
 
 
 
 
 #plot the input and output light curves to see if they match
 fout = w/twopi
 
 tplot = np.arange(tlo,thi+dtresplot,dtresplot)
 xplot = []
 for tnow in tplot:
  xplot.append( np.sum(ckout*np.cos(w*tnow + phase) ) )
 
 if (ploton == 1):
  fig = plt.figure()
  ax1 = fig.add_subplot(211) 
  ax1.errorbar(time,y,sig,ls='',color='k')
  if (costrans==1):
   ax1.plot(tplot,xplot,color='b',label='cosine fit')
  else:
   ax1.plot(tplot,xplot,color='b',label='sine fit')
  ax1.set_xlabel('Time (days)')
  ax1.set_ylabel('flux')
  plt.legend()
  
  ax2 = fig.add_subplot(212) 
  ax2.plot(fout,ckout**2,ls='',marker='o')
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  ax2.set_xlabel('Frequency (cyc/day)')
  ax2.set_ylabel(r'power ($S_k^2 + C_k^2$)')
  plt.savefig('test_r_hesfit_sine_cosine.pdf')
 
 
 
 return(fout, parm, tplot, xplot)
 