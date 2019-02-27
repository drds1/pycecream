import numpy as np
import mylcgen as mlc
import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import scipy.signal as ss







def ccf_simp(lc1_in,lc2_in,dt=1.0,r_crit = 0.8,laglim='auto'):

 lc1 = np.array(lc1_in,dtype=float)
 lc2 = np.array(lc2_in,dtype=float)
 
 #start of function
 #dt = 0.2 #time grid resolution
 
 #only select overlapping periods to take part in time series analysis
 tmin = max(np.min(lc1[:,0]),np.min(lc2[:,0]))#min(np.min(lc1[:,0]),np.min(lc2[:,0]))
 tmax = min(np.max(lc1[:,0]),np.max(lc2[:,0]))#max(np.max(lc1[:,0]),np.max(lc2[:,0]))
 
 
 #interpolate onto regular time grid
 tg  = np.arange(tmin,tmax+dt,dt)
 yg1 = np.interp(tg,lc1[:,0],lc1[:,1])
 yg2 = np.interp(tg,lc2[:,0],lc2[:,1])
 
 yg1_mean = np.mean(yg1)
 yg1_sd   = np.std(yg1)
 yg1 = (yg1 - yg1_mean)/yg1_sd

 yg2_mean = np.mean(yg2)
 yg2_sd   = np.std(yg2) 
 yg2 = (yg2 - yg2_mean)/yg2_sd
 
 
 #compute ccf
 ccf = ss.correlate(yg1,yg2)
 nccf = np.shape(ccf)[0]
 tccf = (np.arange(nccf) - nccf/2)*dt
 
 if (laglim != 'auto'):
  idccf = np.where((tccf > laglim[0]) & (tccf < laglim[1]))[0]
  ccf = ccf[idccf]
  tccf = tccf[idccf]
  nccf = np.shape(ccf)[0]
  #for it in range(nccf):
  # print(it,tccf[it], ccf[it])
  #print(np.argmax(ccf))
 
 
 ccf_max = np.max(ccf)
 
 #Peterson, Edelson 2018 definition of ccf inclusion threshold (>0.8 ccfmax)
 ccf_crit = r_crit*ccf_max
 idinc = np.where(ccf > ccf_crit)[0] 
 
 #weighted mean lag
 #lagmean = np.sum(ccf[idinc]*tccf[idinc])/np.sum(ccf[idinc])
 #laglo = tccf[idinc[0]]
 #laghi = tccf[idinc[-1]]
 
 
 #lagpeak
 idmax = np.argmax(ccf)
 
 try:
  idhi = idmax + np.where(ccf[idmax:] < ccf_crit)[0][0]
 except:
  print('ccf upper lag limit goes off the end of the grid. Re run with higher upper laglim')
  idhi = nccf-1
 try:
  idlo = np.where(ccf[:idmax] < ccf_crit)[0][-1]
 except:
  idlo = 0
  print('ccf lower lag limit goes off the lower end of the grid. Re run with lower lower laglim')
 
 lagpeak = tccf[idmax]
 laglo   = tccf[idlo]
 laghi   = tccf[idhi]
 
 print ('ccf low, peak, hi ',laglo,lagpeak,laghi)
 return(tccf,ccf,laglo,lagpeak,laghi,idlo,idmax,idhi)







#frrss ccf sampler

def ccf_frrss(lc1,lc2,dt=1.0,fraction_rss=0.8,nsim = 500,
centroid_frac = 0.8,flux_randomisation=0):
#
 n1 = np.shape(lc1)[0]
 n2 = np.shape(lc2)[0]
 
 n1pick = np.int(fraction_rss*n1)
 n2pick = np.int(fraction_rss*n2)

 #only select overlapping periods to take part in time series analysis
 tmin = max(np.min(lc1[:,0]),np.min(lc2[:,0]))#min(np.min(lc1[:,0]),np.min(lc2[:,0]))
 tmax = min(np.max(lc1[:,0]),np.max(lc2[:,0]))#max(np.max(lc1[:,0]),np.max(lc2[:,0]))
 tgrid = np.arange(tmin,tmax+dt,dt)

 print ('tmin max',tmin,tmax,dt)
 lagpeak = []
 ccfpeak = []
 ccf_save = []
 lagcent = []
 for i in range(nsim):
  id1 = np.sort(np.random.choice(np.arange(n1),size=n1pick,replace=True))
  id2 = np.sort(np.random.choice(np.arange(n2),size=n2pick,replace=True))
  #print(id1)
  
  if (flux_randomisation==1):
   noise1 = np.random.randn(n1pick)*lc1[id1,2]
   noise2 = np.random.randn(n2pick)*lc2[id2,2]
  else:
   noise1 = np.zeros(n1pick)
   noise2 = np.zeros(n2pick)
   
  x1 = np.interp(tgrid,lc1[id1,0],lc1[id1,1]+noise1)
  x2 = np.interp(tgrid,lc2[id2,0],lc2[id2,1]+noise2)
  x1mean = np.mean(x1)
  x1sd   = np.std(x1)
  x2mean = np.mean(x2)
  x2sd   = np.std(x2) 
  #compute ccf
  #ccf =np.correlate((x2-x2mean)/x2sd,(x1-x1mean)/x1sd,mode='full')#ss.correlate((x1-x1mean)/x1sd,(x2-x2mean)/x2sd)
  x1_new = (x1-x1mean)/x1sd
  x2_new = (x2-x2mean)/x2sd
  auto1 = np.sum(x1_new**2)#np.abs(np.correlate(x1_new,x2_new).max())
  #auto2 = np.correlate((x2-x2mean)/x2sd,(x2-x2mean)/x2sd).max()
  ccf = np.correlate(x1_new,x2_new,mode='full')/auto1
  
  if (i == 0):
   nccf = np.shape(ccf)[0]
   tccf = (np.arange(nccf) - np.floor(nccf/2))*dt
   tmax = tccf[-1]
   tmin = tccf[0]
   lolim = tmin/2
   hilim = tmax/2
   idinc = np.where((tccf < hilim) & (tccf > lolim))[0]
   lag_ccf = tccf[idinc]
   
  idpeak = np.argmax(ccf[idinc])
  ccfp = ccf[idinc][idpeak]
  tpeak = tccf[idinc][idpeak]
  ccf_save.append(ccf[idinc])
  lagpeak.append(tpeak)
  ccfpeak.append(ccfp)
  
  ccfinc = ccf[idinc]
  laginc = tccf[idinc]
  
  idcent = np.where(ccfinc > centroid_frac*ccfp)[0]
  lagcent.append( np.sum(laginc[idcent]*ccfinc[idcent])/np.sum(ccfinc[idcent]) )
  
  
  #fig = plt.figure()
  #ax1 = fig.add_subplot(211)
  #ax1.plot(tgrid,x1,label='drive',color='k')
  #ax1.plot(tgrid,x2,label='echo',color='r')
  #ax2 = fig.add_subplot(212)
  #ax2.plot(ccf,label='drive',color='k')
  #plt.savefig('test_ccf_fig'+np.str(i)+'.pdf')
  
 
  #plt.scatter(lc1[id1,0],lc1[id1,1],color='k')
  #plt.scatter(lc2[id2,0],lc2[id2,1],color='r')
  #plt.show()
 return(np.array(lag_ccf),ccf_save,np.array(lagpeak),np.array(ccfpeak),np.array(lagcent))
 
# 
#
# return(




#t  = np.arange(240)
#nt = t.size
#x =  np.zeros(nt)
#x[10:13] = 1.0
#lag = np.zeros(nt)
#lag[12:15] = 0.5
#
#lc1 = np.array([t,x]).T
#lc2 = np.array([t,lag]).T
#nsim = 500
#lag_ccf,ccf_save,lagpeak,ccfpeak = ccf_frrss(lc1,lc2,dt=0.5,fraction_rss=0.8,nsim =nsim)
#
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(411)
#ax1.plot(t,x)
#ax2 = fig.add_subplot(412)
#ax2.plot(t,lag)
#ax3 = fig.add_subplot(413)
#ax3.hist(lagpeak,bins=50)
#ax4 = fig.add_subplot(414)
##
##ccf = np.correlate((x-x.mean())/x.std(),(lag-lag.mean())/lag.std(),mode='full')
##ax4.plot(lag_ccf,ccf)
#[ax4.plot(lag_ccf,ccf_save[i],color='k',label=None) for i in range(nsim)]
#plt.show()
#







#sample code to test
#ilag = 5
#sigma = 0.3
#
##generate synthetic data
#drive = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=100,dt=1.0,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
#
#t = drive[:,0]
#x = drive[:,1]
#xsd = np.std(x)
#nx = np.shape(x)[0] 
#x = (x - np.mean(x))/xsd
#xsd = np.std(x)
#sigx = np.ones(nx)*xsd*sigma
#x = x + np.random.randn(nx)*sigx
#
#lc1 = np.array([t,x,sigx]).T
#
##define convolution kernel
#conv = np.zeros(nx)
#conv[ilag-1:ilag+1] = 1
#
#
#
#echo  = np.convolve(x, conv, mode='same')
#echo2 = np.convolve(x, conv, mode='full')[:nx]
#
#echo2 = (echo2 - np.mean(echo2))/np.std(echo2) 
#esd = np.std(echo2)
#ne2=np.shape(echo2)[0]
#sige = np.zeros(ne2)+esd*sigma
#echo2 = echo2 + np.random.randn(ne2)*esd*sigma
#lc2 = np.array([t,echo2,sige]).T








