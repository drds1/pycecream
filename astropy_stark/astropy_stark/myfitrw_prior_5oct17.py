import numpy as np
from mylcgen import *
from myrandom import *
import matplotlib.pylab as plt
import myresample as mrs
import myredcisq as mrc
twopi = 2*np.pi

 
def fitrw(timeall,yall,sigall,floin=-1,fhiin=-1,ploton=0,nits_montecarlo=1000):
 
 #stick all light curves together
 nlc = len(timeall)
 lc  = []
 for ilc in range(nlc):
  ntnow = np.shape(timeall[ilc])
  lct = np.ones(ntnow)*ilc
  lc.append(lct)
 lc = np.concatenate(lc)
 time = np.concatenate(timeall)
 y    = np.concatenate(yall)
 sig  = np.concatenate(sig)
 
 ndat = np.shape(y)[0]
 
 parmout = []
 covout = []
 
 
 #ndatall    = []
 #dtall   = []
 #tmins = []
 #tmaxs = []
 #for i in range(nlc):
 # ndat = np.shape(timeall[i])[0]
 # ndatall.append(ndat)
 # tmins.append(np.min(timeall[i]))
 # tmaxs.append(np.max(timeall[i]))
 # dtall.append(np.median(timeall[i][1:] - timeall[i][:-1]))
 #
 #
 #
 #
 tlo = np.min(time)
 thi = np.max(time)
 dt  = np.median(time[1:]-time[:-1])
 
 #input some minimum and maximum frequency range for the Fourier terms
 if (floin < 0):
  flo = 0.5/(thi-tlo)
 else:
  flo = floin
 
 if (fhiin < 0):
  fhi = 0.5/dt
 else:
  fhi = fhiin
 
 
 #angular frequency
 w = np.arange(0.0,fhi+flo,flo)*2*np.pi
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
 nw2 = 2*nw
 
 #sig  = sigall[ilc]
 #y    = yall[ilc]
 #time = timeall[ilc]
 #ndat = ndatall[ilc]

 cvec = np.ones(nw2+1)
 hes = np.ones((nw2+1,nw2+1))
 sig2 = sig*sig
 #save the sine and cosine value in  2d array
 cwt = np.ones((ndat,nw))
 swt = np.ones((ndat,nw))
 for iw in range(nw):
  cwt[:,iw] = np.cos(w[iw]*time)
  swt[:,iw] = np.sin(w[iw]*time)
 
 
 npar = nw2+nlc
 for ik in range(nw):
  hes[ik,:nw] = np.sum( (cwt[:,ik]/sig2*cwt[:,:]).T )
  hes[ik,nw:nw2] = np.sum( (cwt[:,ik]/sig2*swt[:,:]).T )
  hes[ik,nw2:npar] = #DO THIS the extra bit for the extra offsets
 
 for ik in range(nw,nw2,1):
  hes[ik,:nw]      = np.sum( (swt[:,ik]/sig2*cwt[:,:]).T )
  hes[ik,nw:nw2]   = np.sum( (swt[:,ik]/sig2*swt[:,:]).T )
  hes[ik,nw2:npar] = #DO THIS the extra bit for the extra offsets 
 
 for ik in range(nw2,npar,1):
  hes[ik,:nw]      = 
  hes[ik,nw:nw2]   = 
  hes[ik,nw2:npar] =
 
 
 
 
 
  
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



  dtres = dt/10
  t_ext  = 0.0#0.2*(thi - tlo)
  tplot = np.arange(tlo-t_ext,thi+t_ext+dtres,dtres)
  ntplot = np.shape(tplot)[0]
  xplot = 0*tplot
  if (ploton == 1):
   #!!!!!!!!!!!!!!!!!!!!!!! plot the results

   
   xplot = []
   
   #resample using errors to get error snake
   nits = nits_montecarlo
   xplotsave = np.zeros((ntplot,nits))
   xplot0_save = np.zeros((ntplot,nits))
   pspecsave = np.zeros((nw,nits))
   
   
   print 'calculating error sake...'
   for i in range(nits):
    print i,nits
    #cknew = normdis(1,ckout,sdcos)
    #sknew = normdis(1,skout,sdsin)
    #osnew = normdis(1,osout,sdos)
    #need to use multivariate gaussian to sample correlated probability space
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
   
   fig = plt.figure()
   ax1 = fig.add_subplot(211) 
   ax1.errorbar(time,y,sig,ls='',color='k')
   ax1.plot(tplot,xplot0,color='b')
   ax1.plot(tplot,xplot0-xplot0_sd,color='b')
   ax1.plot(tplot,xplot0+xplot0_sd,color='b')
   
   
   
   #calculate the reduced chi squared
   print 'chisq info', mrc.myredcisq(time,y,sig,tplot,xplot,npar = 0)
   
   
   
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
 
  
 return(parmout,covout,freq,tplot,xplot)












#!!!!!!!!!!!!!!!!!!! test the code using this program



timeall = []
sigall  = []
yall    = []

nlc = 1



for i in range(nlc):
 #generate test light curve and add noise
 datpre      = mylcgen(tlo=0,thi=100,dt=0.03,iseed=132423)
 npre     = np.shape(datpre[:,0])[0]
 datmean  = np.std(datpre[:,1])
 snow =  np.ones(npre)/10*datmean 
 dat = mrs.myresample(dir='',fname=[''],dtave=2.0,sampmin=0.8,sampcode=3,datin=np.array((datpre[:,0],datpre[:,1],snow)).T)
 ndat = np.shape(dat[:,0])[0]
 sig = dat[:,2]
 for i in range(ndat):
  dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
 sigall.append( sig )
 yall.append( dat[:,1] +50 )
 timeall.append( dat[:,0] )
 
 
 
a = fitrw(timeall,yall,sigall,ploton=1) 
 


#
#    