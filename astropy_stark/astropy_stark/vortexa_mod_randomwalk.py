import numpy as np
import matplotlib.pylab as plt

#import myredcisq as mrc


twopi = 2*np.pi


#improved speed and added option to fit extra frequency components on top of drw
def fitrw(timeall,yall_in,sigall_in,floin=-1,fhiin=-1,plot_tit='fig_myrwfit',dtresin=-1,nits = 1,tplotlims=[],extra_f=[],
p0=-1,bpl = [0.5,2,2],normalise = 1): 

 if (plot_tit != ''):
  ploton = 1
 else:
  ploton = 0
  
  
 w0 = bpl[0]*twopi

 parmout = []
 covout = []
 
 nlc = len(yall_in)
 ndatall    = []
 dtall   = []
 tmins = []
 tmaxs = []
 yall = []
 sigall = []
 
 sdsave = []
 meansave = []
 for i in range(nlc):
  ndat = np.shape(timeall[i])[0]
  ndatall.append(ndat)
  tmins.append(np.min(timeall[i]))
  tmaxs.append(np.max(timeall[i]))
  dtall.append(np.median(timeall[i][1:] - timeall[i][:-1]))
  
  #normalise data set if option switched on
  yn = yall_in[i]
  sn = sigall_in[i]
  if (normalise == 1):
   sd = np.std(yn)
   mean = np.mean(yn)
   ynew = (yn - mean)/sd
   snew = sn/sd
   yall.append(ynew)
   sigall.append(snew)
   sdsave.append(sd)
   meansave.append(mean)
  else:
   yall.append(yn)
   sigall.append(sn)
   
   
   
 
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
 
 dw = np.mean(w[1:] - w[:-1])
 #print flo,fhi,fhiin, 0.5/(thi-tlo)
 #w = np.arange(0,fhi+flo,5*flo)*2*np.pi
 #w[0] = flo/5000
 
 
 #prior
 if (bpl == []):
  w2 = w**2
  alpha = -2.
 else:
  bplin = list(bpl)
  #bplin[0] = twopi*bpl[0]
  #w2 = (1. + (w/(twopi*bplin[0]))**(bplin[1]-bplin[2])) / (w/(twopi*bplin[0]))**bplin[1]
 bplin = list(bpl)
 w2 = (1. + (w/(twopi*bplin[0]))**(bplin[1]-bplin[2])) / (w/(twopi*bplin[0]))**bplin[1]
 
 w_hi = w[-1]
 w_lo = w[0]
 w2_integrate_hi = w_hi*( (w_hi/(twopi*bplin[0]))**(-bplin[1])/(1.-bplin[1]) + (w_hi/(twopi*bplin[0]))**(-bplin[2])/(1.-bplin[2]) )
 w2_integrate_lo = w_lo*( (w_lo/(twopi*bplin[0]))**(-bplin[1])/(1.-bplin[1]) + (w_lo/(twopi*bplin[0]))**(-bplin[2])/(1.-bplin[2]) )
 w2_integrate = w2_integrate_hi - w2_integrate_lo
 
 
 
 #Mistake in algebra somewhere means we should square root the prior matrix (not sure why)
 w2 = w2
 
 #now compute the hessian matrix
 wold = np.array(w)
 nwold = np.shape(w)[0]
 np2old = 2*nwold
 
 
 #print 'this bit here'
 #the following has to be donw for as many light curves as you have
 
 nextra = len(extra_f)
 w_extra = np.array([ef*twopi for ef in extra_f])
 w = np.append(w,w_extra)
 nw = np.shape(w)[0]
 np2 = 2*nw
 for ilc in range(nlc):
  sig  = sigall[ilc]
  y    = yall[ilc]
  time = timeall[ilc]
  ndat = ndatall[ilc]
 
  #if p0 = -1 then use theoretical statistics and integration to set p0 
  #for random walk time series, rms = p0(whi^(alpha-1) - wlo(alpha-1)) / ((alpha+1) w0^alpha)
  rms = np.std(y)

  
  p0_out  = 2*rms/dw/w2_integrate#rms* ((alpha+1)*wold[0]**alpha / ( wold[-1]**(alpha+1) - wold[0]**(alpha+1)))
  
  
  #parameters for prior
  prior = 1.0#( wold[-1]**(alpha+1) - wold[0]**(alpha+1))/(alpha+1)/dw
  
  
  
  cvec = np.ones(np2+1)
  hes = np.ones((np2+1,np2+1))
  sig2 = sig*sig
  #save the sine and cosine value in  2d array
  cwt = np.ones((ndat,nw))
  swt = np.ones((ndat,nw))
  for iw in range(nw):
   cwt[:,iw] = np.cos(w[iw]*time)
   swt[:,iw] = np.sin(w[iw]*time)
  
  
  hnow_cc = np.tensordot(cwt.T/sig2,cwt,axes=1)#/sig2)
  hnow_sc = np.tensordot(swt.T/sig2,cwt,axes=1)#/sig2)
  
  hnow_cs = np.tensordot(cwt.T/sig2,swt,axes=1)#np.sum(cwt[:,iw]*swt[:,iwp]/sig2)
  hnow_ss = np.tensordot(swt.T/sig2,swt,axes=1)#np.sum(swt[:,iw]*swt[:,iwp]/sig2)
 
  
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
  
  #print 'calculating cov matrix'
  
  
  
  
  for ikp in range(1,np2+1,2):
   wp = w[iwp]
   wp2 = wp*wp
   cvec[ikp]   = np.sum(y*cwt[:,iwp]/sig2)
   cvec[ikp+1] = np.sum(y*swt[:,iwp]/sig2)
   iw = 0
   #offset
   hes[ikp,0] = np.sum(cwt[:,iwp]/sig2)
   hes[ikp+1,0] = np.sum(swt[:,iwp]/sig2) 
   for ik in range(1,np2+1,2):
    #cosine column
    hes[ikp,ik]   = hnow_cc[iw,iwp]
    hes[ikp,ik+1] = hnow_sc[iw,iwp]
    if ((ik == ikp) & (ik <= np2old) ):
     hes[ikp,ik] = hnow_cc[iw,iwp] + w2_integrate/(2*rms*w2[iwp])#1./(p0*dw*w2[iwp])#*prior
    #sine column
    hnow_c = hnow_cs[iw,iwp]#np.sum(cwt[:,iw]*swt[:,iwp]/sig2)
    hnow_s = hnow_ss[iw,iwp]#np.sum(swt[:,iw]*swt[:,iwp]/sig2)
    hes[ikp+1,ik]   = hnow_cs[iw,iwp]
    hes[ikp+1,ik+1] = hnow_ss[iw,iwp]
    if ((ik == ikp) & (ik <= np2old)):
     hes[ikp+1,ik+1] = hnow_ss[iw,iwp] + w2_integrate/(2*rms*w2[iwp])
    
    iw = iw + 1  
   iwp = iwp + 1
        
  print('covariance')
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
 
 
  if (isinstance(dtresin, np.ndarray) == 1):
   tplot = 1.*dtresin
   dtres = tplot[1] - tplot[0] 
  elif (dtresin == 0):
   dtres = dt/10
   tplot = np.arange(tlo-t_ext,thi+t_ext+dtres,dtres)
  elif (dtresin == -1):
   tplot = 1.*timeall[ilc]
  else:
   dtres = dtresin
   tplot = np.arange(tlo-t_ext,thi+t_ext+dtres,dtres)
   
  
  if (tplotlims != []):
   tplot = np.arange(tplotlims[0],tplotlims[1],tplotlims[2])
   
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
  
  
  print('calculating error sake...')
  #draw paraeters from multivariate gaussian
  
  #need to limit the size of the parameters to sample 
  # with multivariate gaussian otherwise will take forever
  #include max of 2000 plus the extra frequencies
  ipmax = min(500,np2+1)
  #sensibly sort the covariance matrix to identify the most correlated parameteres
  #for multivariate gaussian stepping
  cor = np.abs(np.corrcoef(cov))
  corabs_sum = np.sum(np.abs(cor),axis=0)/(np2+1)
  ipinc = np.argsort(corabs_sum)[-1::-1]
  #print corabs_sum[ipinc][:10]
  #print 'average correlation coeff (abs) for mvg'
  #print 'corabssum[:10]'
  #print corabs_sum[ipinc][:min(10,np2+1)]
  #print ''
  #
  #print 'corabssum[-10:]'
  #print corabs_sum[ipinc][max(0,np2+1-10):]
  #print ''
  #
  #print 'corabssum chop'
  #print corabs_sum[ipinc][ipmax-1]
  #print ''
  
  ipinc = ipinc[:ipmax]
  #idsort = np.dstack(np.unravel_index(np.argsort(cor.ravel()), (np2old+1, np2old+1)))
  
  
  #ipinc = np.append(np.arange(ipmax),np.arange(np2old+1,np2+1))
  ipdel = np.array([i for i in range(np2+1) if i not in ipinc])
  #ipdel = np.arange(ipmax,np2old+1,1)
  psd = np.sqrt(np.diag(cov))
  parmnew = 1.*parm
 
  
  
  cwhires = np.ones((ntplot,nw))
  swhires = np.ones((ntplot,nw))
  for iw in range(nw):
   cwhires[:,iw] = np.cos(w[iw]*tplot)
   swhires[:,iw] = np.sin(w[iw]*tplot)
  
  
  if (np.sum(sig) == ndat):
   print('uniform error bars - auto scaling')
   sknew = parm[1::2]
   cknew = parm[2::2]
   osnew = parm[0]
   xptemp = osnew + np.dot(swhires,sknew) + np.dot(cwhires,cknew)
   xitp = np.interp(time,tplot,xptemp)
   var = np.sum((xitp - y)**2)/ndat
   cov = cov * var
 
  if (nits > 1):
   covnew = np.delete(cov,ipdel,axis=0)
   covnew = np.delete(covnew,ipdel,axis=1)
   print(np.shape(covnew),np.shape(parm),np.shape(parm[ipinc]))
   parmnow = np.random.multivariate_normal(parm[ipinc],covnew,size=nits)
   
  
  
  for i in range(nits):
   #print('errorbars',i,nits)
   
   #cknew = normdis(1,ckout,sdcos)
   #sknew = normdis(1,skout,sdsin)
   #osnew = normdis(1,osout,sdos)
   #need to use multivariate gaussian to sample correlated probability space
   
   if (nits > 1):
    parmnew = parm + np.random.randn(np2+1)*psd
    parmnew[ipinc] = parmnow[i,:]
   
   osnew = parm[0]#parmnew[0]
   cknew = parmnew[1::2]
   sknew = parmnew[2::2]
   
   #cknew[0] = ckout[0]
   #sknew[0] = skout[0]
   #print('osnew',osnew,parm[0])
   pspecsave[:,i] = sknew*sknew + cknew*cknew
   
   xplot0_save[:,i] = osnew# + cknew[0]*np.cos(w[0]*tplot) + sknew[0]*np.sin(w[0]*tplot)
   
   nk = np.shape(cknew)[0]
   #cktile = np.tile(cknew,ntplot).reshape(nk,ntplot).T
   #sktile = np.tile(sknew,ntplot).reshape(nk,ntplot).T
   #print(np.shape(cktile),np.shape(swhires))
   xplotsave[:,i] = osnew + np.dot(swhires,sknew) + np.dot(cwhires,cknew)
   
   
   #print(np.shape(np.sum(sktile*swhires,axis=1)))
   
   
   
   #a = []
   #for it in range(ntplot):
   # tnow = tplot[it]
   # a.append( osnew + np.sum(cknew*cwhires[w*tnow) + sknew*np.sin(w*tnow)) )
    #print(it,tnow,a[it],xplotsave[it,i])
   #input()
    
  if (normalise == 1):
   #xplot = xplot*sdsave[0] + meansave[0]
   #xplotsd = xplotsd*sdsave[0]
   #xplotlo = xplotlo*sdsave[0] + meansave[0]
   #xplothi = xplothi*sdsave[0] + meansave[0]
   xplotsave = xplotsave*sdsave[0] + meansave[0]
   y = y*sdsave[0] + meansave[0]
   sig = sig*sdsave[0]
  
      
    
  #xplot = np.mean(xplotsave,axis=1)
  xplot0_sd = np.ones(ntplot)*sdos
  
  
  
  xplotlo,xplot,xplothi = np.percentile(xplotsave,[15.865,50,84.135],axis=1)      
  xplotsd = (xplothi-xplotlo)/2
  
  #xplotsd = np.std(xplotsave,axis=1)# - xplot0_sd
  
  
  #plot errors for power spectrum also
  pspec = ckout**2 + skout**2
  pspeclo,pspec,pspechi = np.percentile(pspecsave,[15.865,50,84.135],axis=1) 
  #pspecsd = np.std(pspecsave,axis=1)
  
  #plot the input and output light curves to see if they match
  xplot0 = np.ones(ntplot)*osout#ckout[0]*np.cos(w[0]*tplot) + skout[0]*np.sin(w[0]*tplot)
  #for tnow in tplot:
  # xplot.append( np.sum(ckout*np.cos(w*tnow) + skout*np.sin(w*tnow)) )
  if (ploton == 1):
   fig = plt.figure()
   ax1 = fig.add_subplot(311) 
   ax1.errorbar(time,y,sig,ls='',color='k')
   ax1.plot(tplot,xplot0,color='b',ls='--')
   ax1.plot(tplot,xplot0-xplot0_sd,color='b',ls='--')
   ax1.plot(tplot,xplot0+xplot0_sd,color='b',ls='--')
   
   
   #plot data separately
   axd = fig.add_subplot(312)
   axd.scatter(time,y,color='k')
   print('errorbars')
   print(sig[:10])
   print(y[:10])
   #calculate the reduced chi squared
   #print 'chisq info', mrc.myredcisq(time,y,sig,tplot,xplot,npar = 0)
   
   
   
   #for i in range(nits):
   # ax1.plot(tplot,xplotsave[:,i],color='b',alpha=0.06)
   ax1.fill_between(tplot,xplotlo,xplothi,alpha=0.2,color='k',label=None)
   
   ax1.plot(tplot,xplot,color='b')
   ax1.set_xlabel('Time (days)')
   ax1.set_ylabel('flux')
   ylim = [np.min(y),np.max(y)]
   ax1.set_ylim(ylim)
   ax2 = fig.add_subplot(313) 
   
   ax2.plot(freq,pspec,ls='',marker='o',color='k')
   ax2.fill_between(freq,pspeclo,pspechi,color='k',alpha=0.2)
   #ax2.plot(freq,sdsin,ls='',color='g',marker='o')
   #ax2.plot(freq,sdcos,ls='',color='b',marker='o')
   ax2.set_xscale('log')
   ax2.set_yscale('log')
   ax2.set_xlabel('Frequency (cyc/day)')
   ax2.set_ylabel(r'power ($S_k^2 + C_k^2$)')
   ax2.grid(True)
   
   if (plot_tit == 'show'):
    plt.show()
   else:
    plt.savefig(plot_tit+'_'+np.str(ilc)+'.png')
   
 sig2_prior = (2*rms*w2)/w2_integrate
 

 return(parmout,covout,freq,tplot,xplot,xplotsd,xplotlo,xplothi,p0_out,w0,dw,sig2_prior,xplotsave)










#####!!!!!!!!!!!!!!!!!!! test the code using this program
####
######
#from mylcgen import *
#from myrandom import *
#import myresample as mrs
#timeall = []
#sigall  = []
#yall    = []
#
#nlc = 1
#
#tlo = 0.0
#thi = 100
#dtave = 1.0
#dtres = 0.1
#
#for i in range(nlc):
# #generate test light curve and add noise
# datpre      = mylcgen(tlo=tlo,thi=thi,dt=dtres,iseed=132423)
# npre     = np.shape(datpre[:,0])[0]
# datmean  = np.std(datpre[:,1])
# snow =  np.ones(npre)/10*datmean 
# dat = mrs.myresample(dir='',fname=[''],dtave=dtave,sampmin=0.5,sampcode=3,datin=np.array((datpre[:,0],datpre[:,1],snow)).T)
# ndat = np.shape(dat[:,0])[0]
# sig = dat[:,2]
# for i in range(ndat):
#  dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
# sigall.append( sig )
# yall.append( dat[:,1] +50 )
# timeall.append( dat[:,0] )
# 
# 
# #dat =mylcgen(tlo=tlo,thi=thi,dt=dtave,iseed=132423)
# #nd = np.shape(dat[:,0])[0]
# #timeall = [dat[:,0]]
# #yall = [dat[:,1]]
# #sigall = [np.ones(nd)/10*datmean ]
# 
#print('going in')
#
#forecast = 10.0
#dtres = np.arange(tlo-50.0,thi+forecast,dtres)
#a = fitrw(timeall,yall,sigall,plot_tit='show',floin=0.1/(thi-tlo),fhiin = 1./dtave,dtresin = dtres,
#bpl=f[0.1,2.0,2.0],
#nits = 1,extra_f=[]) 
#
#
#
#
###plot covariances
##cov = a[1][0]
##sig = np.diag(cov)[1:]
##nsig = np.shape(sig)[0]
##vars = np.array([sig[i]+sig[i+1] for i in range(0,nsig,2)])
##
##
##ncov = np.shape(cov)[0]
##freq = a[2]
##plt.clf()
##fig = plt.figure()
##ax1 = fig.add_subplot(111)
##
##ax1.plot(freq,vars)
##ax1.set_xscale('log')
##ax1.set_yscale('log')
##
#plt.show()






   