#fits random walk model to set of light curves
#INPUTS:      
#timall      : [timeall[t1[1..nt1],t2[1..n2],...tnlc[1..ntnlc]],
#yall        : [y1[1..nt1],y2[1..n2],...ynlc[1..ntnlc]],
#sigall      : [sig1[1..nt1],sig2[1..n2],...signlc[1..ntnlc]]
#floin,fhiin : the lower and upper frequencies of the random walk
#ploton      : make diagnostic plot of model fit over data
# dtresin    : how hi res do you want the model (cadence days) if -1, just evaluate on same time grid as data
# nits       : get model errors using monte carlo (set hi (1000) if you need the error envelopes around the model fit or 1 to skip if you dont need the error envelopes) 

#OUTPUTS:
#:parmout      : the output Sk and Ck parameters for the model (f(t) = sum_nw S_k sin(wk t) + C_k cos(wk t))
#:freq         : The frequencies at which the model is evaluated 1...nw
#tplot,xplot   : the high resolution model (resolution specified by dtresin argument)
#covout        : covout[:nparameters,:nparameters] the output covariance matrix for the model fit

import numpy as np
from mylcgen import *
import matplotlib.pylab as plt

twopi = 2*np.pi

#if dtresin = -1 then just use model evaluated at data point time 
def fitrw(timeall,yall,sigall,floin=-1,fhiin=-1,ploton=0,dtresin=-1,nits = 1):
 
 parmout = []
 covout = []
 
 nlc = len(yall)
 ndatall    = []
 dtall   = []
 tmins = []
 tmaxs = []
 for i in range(nlc):
  ndat = np.shape(timeall[i])[0]
  ndatall.append(ndat)
  tmins.append(np.min(timeall[i]))
  tmaxs.append(np.max(timeall[i]))
  dtall.append(np.median(timeall[i][1:] - timeall[i][:-1]))
 
 
 
 
 tlo = np.min(tmins)
 thi = np.max(tmaxs)
 dt  = np.min(dtall)
 
 #input some minimum and maximum frequency range for the Fourier terms
 if (floin < 0):
  flo = 0.5/(thi-tlo)
 else:
  flo = floin
 
 if (fhiin < 0):
  fhi = 2.0/dt
 else:
  fhi = fhiin
 
 
 #angular frequency
 w = np.arange(flo,fhi+flo,flo)*2*np.pi
 #print flo,fhi,fhiin, 0.5/(thi-tlo)
 #w = np.arange(0,fhi+flo,5*flo)*2*np.pi
 #w[0] = flo/5000
 bpl = [0.5,-2,-2]
 
 #prior
 if (bpl == []):
  w2 = w**2
 else:
  bplin = list(bpl)
  bplin[0] = twopi*bpl[0]
  w2 = (1. + (w/bplin[0])**(bplin[1]-bplin[2])) / (w/bplin[0])**bplin[1]
 
 #Mistake in algebra somewhere means we should square root the prior matrix (not sure why)
 w2 = np.sqrt(w2)
 
 #now compute the hessian matrix
 nw = np.shape(w)[0]
 np2 = 2*nw
 
 
 
 #the following has to be donw for as many light curves as you have
 for ilc in range(nlc):
  sig  = sigall[ilc]
  y    = yall[ilc]
  time = timeall[ilc]
  ndat = ndatall[ilc]
 
  cvec = np.ones(np2+1)
  hes = np.ones((np2+1,np2+1))
  sig2 = sig*sig
  #save the sine and cosine value in  2d array
  cwt = np.ones((ndat,nw))
  swt = np.ones((ndat,nw))
  for iw in range(nw):
   cwt[:,iw] = np.cos(w[iw]*time)
   swt[:,iw] = np.sin(w[iw]*time)
  
  
  
  
  #offset
  #cos
  iwp = 0
  hes[0,0] = np.sum(1./sig2)
  cvec[0] = np.sum(y/sig2)
  for ik in range(1,np2+1,2):
   hes[0,ik] = np.sum(cwt[:,iwp]/sig2)
   hes[0,ik+1] = np.sum(swt[:,iwp]/sig2)
   iwp = iwp + 1
  
  iwp = 0
  for ikp in range(1,np2+1,2):
   wp = w[iwp]
   
   wp2 = wp*wp
   cvec[ikp]   = np.sum(y*cwt[:,iwp]/sig2)
   cvec[ikp+1] = np.sum(y*swt[:,iwp]/sig2)
   #print iwp,nw
   iw = 0
   
   #offset
   hes[ikp,0] = np.sum(cwt[:,iwp]/sig2)
   hes[ikp+1,0] = np.sum(swt[:,iwp]/sig2)
   
   for ik in range(1,np2+1,2):
    #cosine column
    hnow_c = np.sum(cwt[:,iw]*cwt[:,iwp]/sig2)
    hnow_s = np.sum(swt[:,iw]*cwt[:,iwp]/sig2)
    if (ik == ikp):
     hnow_c = hnow_c + w2[iwp]
    hes[ikp,ik]   = hnow_c
    hes[ikp,ik+1] = hnow_s 
     
    #sine column
    hnow_c = np.sum(cwt[:,iw]*swt[:,iwp]/sig2)
    hnow_s = np.sum(swt[:,iw]*swt[:,iwp]/sig2)
    if (ik == ikp):
     hnow_s = hnow_s + w2[iwp]
    hes[ikp+1,ik]   = hnow_c
    hes[ikp+1,ik+1] = hnow_s 
    iw = iw + 1
    
   iwp = iwp + 1
        
  #invert the matrix to get the covariance matrix
  cov = np.linalg.inv(hes)
  
  
  #find the parameters
  parm  = cov.dot(cvec)
  skout = parm[2::2]
  ckout = parm[1::2]
  osout = parm[0]
  
  
  
  
  var = np.diag(cov)
  sd  = np.sqrt(var)#/1.6
  sdos  = sd[0]
  sdcos = sd[1::2]
  sdsin = sd[2::2]
  
  
  
  parmout.append(parm)
  covout.append(cov)

  freq = w/2/np.pi


  t_ext  = 0.0#0.2*(thi - tlo)

  if (dtresin == 0):
   dtres = dt/10
   tplot = np.arange(tlo-t_ext,thi+t_ext+dtres,dtres)
  elif (dtresin == -1):
   tplot = 1.*timeall[ilc]
  elif (isinstance(dtresin, np.ndarray) == 1):
   tplot = 1.*dtresin
   dtres = tplot[1] - tplot[0]
  else:
   dtres = dtresin
   tplot = np.arange(tlo-t_ext,thi+t_ext+dtres,dtres)
   
  
  
  ntplot = np.shape(tplot)[0]
  xplot = 0*tplot
  #if (ploton == 1):
  #!!!!!!!!!!!!!!!!!!!!!!! plot the results

  
  xplot = []
  
  #resample using errors to get error snake
  #nits = 10
  xplotsave = np.zeros((ntplot,nits))
  xplot0_save = np.zeros((ntplot,nits))
  pspecsave = np.zeros((nw,nits))
  
  
  print 'calculating error sake...'
  for i in range(nits):
   #print i,nits
   #cknew = normdis(1,ckout,sdcos)
   #sknew = normdis(1,skout,sdsin)
   #osnew = normdis(1,osout,sdos)
   #need to use multivariate gaussian to sample correlated probability space
   if (nits == 1):
    parmnew = 1.*parm
   else:
    parmnew = np.random.multivariate_normal(parm,cov)
   osnew = parmnew[0]
   cknew = parmnew[1::2]
   sknew = parmnew[2::2]
   
   #cknew[0] = ckout[0]
   #sknew[0] = skout[0]
   pspecsave[:,i] = sknew*sknew + cknew*cknew
   
   xplot0_save[:,i] = osnew# + cknew[0]*np.cos(w[0]*tplot) + sknew[0]*np.sin(w[0]*tplot)
   for it in range(ntplot):
    tnow = tplot[it]
    xplotsave[it,i] =  osnew + np.sum(cknew*np.cos(w*tnow) + sknew*np.sin(w*tnow))
    
  xplot = np.mean(xplotsave,axis=1)
  xplot0_sd = np.ones(ntplot)*sdos
  
  xplotsd = np.std(xplotsave,axis=1)# - xplot0_sd
  
  
  #plot errors for power spectrum also
  pspec = ckout**2 + skout**2
  pspecsd = np.std(pspecsave,axis=1)
  
  #plot the input and output light curves to see if they match
  xplot0 = np.ones(ntplot)*osout#ckout[0]*np.cos(w[0]*tplot) + skout[0]*np.sin(w[0]*tplot)
  #for tnow in tplot:
  # xplot.append( np.sum(ckout*np.cos(w*tnow) + skout*np.sin(w*tnow)) )
  if (ploton == 1):
   fig = plt.figure()
   ax1 = fig.add_subplot(211) 
   ax1.errorbar(time,y,sig,ls='',color='k')
   ax1.plot(tplot,xplot0,color='b')
   ax1.plot(tplot,xplot0-xplot0_sd,color='b')
   ax1.plot(tplot,xplot0+xplot0_sd,color='b')
   
   
   
   #calculate the reduced chi squared
   #print 'chisq info', mrc.myredcisq(time,y,sig,tplot,xplot,npar = 0)
   
   
   
   #for i in range(nits):
   # ax1.plot(tplot,xplotsave[:,i],color='b',alpha=0.06)
   ax1.fill_between(tplot,xplot-xplotsd,xplot+xplotsd,alpha=0.2,color='k',label=None)
   
   ax1.plot(tplot,xplot,color='b')
   ax1.set_xlabel('Time (days)')
   ax1.set_ylabel('flux')
   
   
   ax2 = fig.add_subplot(212) 
   
   ax2.plot(freq,pspec,ls='',marker='o',color='k')
   ax2.fill_between(freq,pspec-pspecsd,pspec+pspecsd,color='k',alpha=0.2)
   #ax2.plot(freq,sdsin,ls='',color='g',marker='o')
   #ax2.plot(freq,sdcos,ls='',color='b',marker='o')
   ax2.set_xscale('log')
   ax2.set_yscale('log')
   ax2.set_xlabel('Frequency (cyc/day)')
   ax2.set_ylabel(r'power ($S_k^2 + C_k^2$)')
   ax2.grid(True)
   plt.savefig('test_r_hesfit_'+np.str(ilc)+'.pdf')
 
  
 return(parmout,covout,freq,tplot,xplot,xplotsd)












