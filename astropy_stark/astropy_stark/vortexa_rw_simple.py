import numpy as np



def rw(ti,yi,si=0,tgi = -1,fbreak=-1,custom_freqs=[]):

 print(tgi[1:10])
 print(custom_freqs[:10])
 #initialise starting inputs and normalise
 twopi = 2*np.pi
 t = np.array(ti)
 tlo = np.min(t)
 thi = np.max(t)
 dt  = t[1] - t[0]
 ymean = np.mean(yi)
 ystd  = np.std(yi)
 y = (yi - ymean)/ystd
 ny = len(y)
 rms = np.std(y)
 if (type(si) == np.ndarray):
  s = si/ystd
 else:
  s = np.ones(ny)
 if (type(tgi) == np.ndarray): 
  tgrid = np.array(tgi)
 else:
  tgrid = np.arange(tlo,thi+dt,dt)
 ntgrid = len(tgrid)
 dtm = tgrid[1]-tgrid[0]
 tglo,tghi = tgrid[0],tgrid[-1]
 
 #initialise frequencies  
 fnyq = 0.5/dtm
 flo = 1.0/(tghi - tglo)
 fhi = fnyq 
 if (fbreak == -1):
  w0  = min(10*flo,fhi/2)*twopi
 else:
  w0 = twopi*fbreak
  
 if (custom_freqs == []): 
  f = np.arange(flo,fhi+flo,flo)
 else:
  f = custom_freqs
 
 w =f*twopi
 nw = np.shape(w)[0]
 np2 = 2*nw
 #define power spectrum prior
 w2 = (w/w0)**(-2)
 w2_integrate = w0**2 * (1/w[0] - 1/w[-1])
 #w2 = 1./(1+(w/w0)**2)
 #w2_integrate = w0*(1./np.tan(w[-1]/w0) - 1./np.tan(w[0]/w0))

 
 #bplin=[fbreak,2,2]
 #w2 = (1. + (w/(twopi*bplin[0]))**(bplin[1]-bplin[2])) / (w/(twopi*bplin[0]))**bplin[1]
 #
 #w_hi = w[-1]
 #w_lo = w[0]
 #w2_integrate_hi = w_hi*( (w_hi/(twopi*bplin[0]))**(-bplin[1])/(1.-bplin[1]) + (w_hi/(twopi*bplin[0]))**(-bplin[2])/(1.-bplin[2]) )
 #w2_integrate_lo = w_lo*( (w_lo/(twopi*bplin[0]))**(-bplin[1])/(1.-bplin[1]) + (w_lo/(twopi*bplin[0]))**(-bplin[2])/(1.-bplin[2]) )
 #w2_integrate = w2_integrate_hi - w2_integrate_lo
 
 
 
 
 
 
 
 #save the sine and cosine value in  2d array to speed computation
 cwt = np.ones((ny,nw))
 swt = np.ones((ny,nw))
 for iw in range(nw):
  cwt[:,iw] = np.cos(w[iw]*t)
  swt[:,iw] = np.sin(w[iw]*t)
 
 
 #perform hessian fit
 cvec = np.ones(np2)
 hes = np.ones((np2,np2))
 s2 = s*s
 cwtT = cwt.T
 swtT = swt.T
 
 hnow_cc = np.tensordot(cwtT/s2,cwt,axes=1)
 hnow_sc = np.tensordot(swtT/s2,cwt,axes=1)
 hnow_cs = np.tensordot(cwtT/s2,swt,axes=1)
 hnow_ss = np.tensordot(swtT/s2,swt,axes=1)

  
  
 import time
 t1 = time.time()
 idk = np.arange(0,np2,2)
 idk1 = np.arange(1,np2,2)
 hes[::2,::2] = hnow_cc[:,:].T
 hes[::2,1::2] = hnow_sc[:,:].T
 hes[1::2,::2] = hnow_cs[:,:].T
 hes[1::2,1::2] = hnow_ss[:,:].T 
 ##add on prior
 x = np.sqrt(w2_integrate)/(2*rms*w2)

 
 
 
 idk = np.arange(0,np2,2)
 idk1 = np.arange(1,np2,2)
 hes[::2,::2] = hnow_cc[:,:].T
 hes[::2,1::2] = hnow_sc[:,:].T
 hes[1::2,::2] = hnow_cs[:,:].T
 hes[1::2,1::2] = hnow_ss[:,:].T 
 
 ##add on prior unless fitting custom frequencies 
 if (custom_freqs == []):
  x = np.sqrt(w2_integrate)/(2*rms*w2)
 else:
  x = 0
  
 hes[idk,idk] = hes[idk,idk] + x
 hes[idk1,idk1] = hes[idk1,idk1] + x
 ys2 = y/s2
 cvec[idk] = np.dot(cwtT,ys2)
 cvec[idk1] = np.dot(swtT,ys2)

 t2 = time.time()
 
 print('time1',t2-t1)

 
 #t1 = time.time()
 #iwp = 0
 #for ikp in range(0,np2,2):
 # wp = w[iwp]
 # wp2 = wp*wp
 # cvec[ikp]   = np.sum(y*cwt[:,iwp]/s2)
 # cvec[ikp+1] = np.sum(y*swt[:,iwp]/s2)
 # iw = 0
 # for ik in range(0,np2,2):
 #  #cosine column
 #  hes[ikp,ik]   = hnow_cc[iw,iwp]
 #  hes[ikp,ik+1] = hnow_sc[iw,iwp]
 #  #sine column
 #  hes[ikp+1,ik]   = hnow_cs[iw,iwp]
 #  hes[ikp+1,ik+1] = hnow_ss[iw,iwp]
 #  #priors
 #  if (ik == ikp):
 #   hes[ikp,ik] = hnow_cc[iw,iwp] +  np.sqrt(w2_integrate)/(2*rms*w2[iwp]) 
 #   hes[ikp+1,ik+1] = hnow_ss[iw,iwp] + np.sqrt(w2_integrate)/(2*rms*w2[iwp]) 
 #  iw = iw + 1  
 # iwp = iwp + 1
 #t2 = time.time()
 #print('time2',t2-t1)
 
 #print('hes')
 #for i in range(5):
 # print(hes[i,:5])
 # 
 # 
 #print('hes2')
 #for i in range(5):
 # print(hes2[i,:5])
 #input()
 
 #invert the matrix to get the covariance matrix
 cov = np.linalg.inv(hes)
 

 
 #find the parameters
 parm  = cov.dot(cvec)
 skout = parm[1::2]
 ckout = parm[0::2]
 ymodout = np.dot(swt,skout) + np.dot(cwt,ckout)
 var = np.sum((ymodout - y)**2)/(ny - 1)
 
 if (type(si) != np.ndarray):
  cov = cov * var

 
 #monte carlo sample the fit to estimate uncertaintie
 nits = 2000
 cwhires = np.ones((ntgrid,nw))
 swhires = np.ones((ntgrid,nw))
 for iw in range(nw):
  cwhires[:,iw] = np.cos(w[iw]*tgrid)
  swhires[:,iw] = np.sin(w[iw]*tgrid)
 parmnow = np.random.multivariate_normal(parm,cov,size=nits)
 ygridsave = np.zeros((ntgrid,nits))
 for i in range(nits):
  parmnew = parmnow[i,:]
  cknew = parmnew[0::2]
  sknew = parmnew[1::2]
  ygridsave[:,i] = np.dot(swhires,sknew) + np.dot(cwhires,cknew)
 ygridsave = ygridsave*ystd + ymean 
 ygridlo,ygrid,ygridhi = np.percentile(ygridsave,[15.865,50,84.135],axis=1)      
 ygridsd = (ygridhi-ygridlo)/2
 
 ygridop = np.array([ygridlo,ygrid,ygridhi]).T
 
 
 #fit stat
 ymod_itp = np.interp(ti,tgrid,ygrid)
 rcoef = np.corrcoef(ymod_itp,yi)[0,1]
 
 #import matplotlib.pylab as plt
 #plt.clf()
 #fig = plt.figure()
 #ax1 = fig.add_subplot(111)
 #ax1.scatter(t,yi,color='r')
 #ax1.plot(tgrid,ygridop[:,1],color='b')
 #ax1.fill_between(tgrid,ygridop[:,0],ygridop[:,2],alpha=0.3,color='b')
 #ax1.set_xlabel('time')
 #plt.show()
 #plt.clf()
 #np.savetxt('rw_mod_bad_series.dat',np.array([t,s,y]).T)

 return(tgrid,ygridop,f,ckout,skout,rcoef)
 
 
 
 





#test

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
#y = yall[0]
#t = timeall[0]
#s = sigall[0]
#
##dat =mylcgen(tlo=tlo,thi=thi,dt=dtave,iseed=132423)
##nd = np.shape(dat[:,0])[0]
##timeall = [dat[:,0]]
##yall = [dat[:,1]]
##sigall = [np.ones(nd)/10*datmean ]
#
#print('going in')
#forecast = 0
#dtm = 0.01
#tmod = np.arange(tlo-forecast,thi+forecast,dtm)
#
##dat = np.loadtxt('/Users/david/github_datascience/vortexa/analytics/arbitrage_flow_study/rw_mod_bad_series.dat')
##t = dat[:,0]
##y = dat[:,1]
##s = dat[:,2]
##tlo = np.min(t)
##thi = np.max(t)
##tmod = np.arange(tlo,thi,10.)
#fc = np.arange(1./100,1./4,1./100)
##fc = np.array([1./10,1./20.,1./40,1./1000,1./30,1./5])
#
#tg,xg,w,ckout,skout,rcoef = rw(t,y,si=s,tgi = tmod,custom_freqs=fc)
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax1.errorbar(t,y,s)
#ax1.plot(tg,xg[:,1])
##ax1.fill_between(tg,xg[:,0],xg[:,2],alpha=0.5)
#
#ax1 = fig.add_subplot(212)
#ax1.plot(w,ckout**2+skout**2)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#
#plt.show()
###
##