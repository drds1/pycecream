import numpy as np
import matplotlib.pylab as plt
import mytemp0 as mt0
import scipy.stats as ss

twopi = np.pi*2
fourpi = 2*twopi
deg2rad = np.pi/180
planck = 6.626307004e-34
c      = 2.99792458e8
boltz  = 1.38064852e-23
uw     = 1.e-10


def pytfb_sub(taus,rgrid,ygrid,zx,tv4,ti4,cosang,wav):
 tauhi  = taus[-1]
 taulo  = taus[0]
 ntaus  = np.shape(taus)[0]
 psis   = np.zeros(ntaus)
 ux = planck * c / ( uw * boltz * wav)


 nrad = np.shape(rgrid)[0] 
 #each point will have a time delay and a weighting append these to a list
 #for each point in the disc
 
 #loop of azimuths
 delsave = []
 wsave   = []
 if (diagnose == 1):
  azsave = []
  rsave = []
 

 for i in range(nrad-1):
  ynow = ygrid[i]
  r    = rgrid[i]
  zsuby=zx-ynow
  dd=zsuby*zsuby+r*r
  d=sqrt(dd)

  #quick way 
  canow     = cosang[i]
  tv4now    = tv4[i]
  ti4now    = tiv[i]*canow
  
  ttot4now   = tv4now + ti4now
  
  T2         = np.sqrt(ttot4now)
  T          = np.sqrt(T2)
  T3         = T*T2



  

  radlo = r
  radhi = rgrid[i+1]
  rhalf = (radlo + radhi)/2
  ylo   = ynow
  yhi   = ygrid[i+1]
  dy    = yhi - ylo
  drad  = radhi - radlo
  
  #now azimuth grid
  azwidth = drad/radlo
  daz = azwidth
  azgrid = np.arange(0.0, twopi,azwidth)
  naz = np.shape(azgrid)[0]
  nazsub1 = naz - 1
  azgrid_s = azgrid[1:] - azgrid[:-1]
  
  #print drad, radhi, azwidth,'sdfsd',radhi*azwidth
  #print radlo, radhi, azwidth, drad, naz,'sdfsf'
  
  raz   = np.random.uniform(radlo,radhi,nazsub1)   
  yaz  = ylo + (raz - radlo)/drad*dy#*np.sqrt(raz*raz + zx2)
  az   = np.random.uniform(low=azgrid[:-1],high=azgrid[1:],size=nazsub1)#np.random.uniform(low=0,high=1,size=nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
  caz  = np.cos(az)
  tdel=( zsuby*ci - raz*caz + d ) # calculate time delay in days

  #calculate the weights
  x=ux/T        # record h c /(lamda K T)
  integrand = fourpi*x*x/(np.cosh(x)-1.)    # the xx/(cosh(x)-1) is part of dBdT !! the integrand of the transfer function£
  solid_angbit = np.abs(np.cos(np.arctan(dy/drad)))
  dsa        = rhalf*dr*daz * solid_angbit
  dTdL      = 1./T3 * canow / dd
  weight=integrand*dsa*dTdL
  
  wsave.append([weight]*nazsub1)
  delsave.append(tdl)
  if (diagnose == 1):
   azsave.append(az)
   rsave.append(raz)


 
 #introduce a lower threshold weight to abort the iterations (to save time)
 
 
 
 delsave = [item for sublist in delsave for item in sublist]
 wsave   = [item for sublist in wsave for item in sublist]
 if (diagnose == 1):
  rsave   = [item for sublist in rsave for item in sublist]
  azsave  = [item for sublist in azsave for item in sublist]



 delsave = np.array(delsave)
 wsave   = np.array(wsave)
 if (diagnose == 1):
  rsave   = np.array(rsave)
  azsave  = np.array(azsave)
 
 nds = np.shape(delsave)[0]
 
 
 
 #add smoothing function to approximate delta function
 #ntaus = np.shape(taus)[0]
 sigtau = 5*dtau/2
 sigtaulim = 3*sigtau
 for i in range(ntaus):
  taunow = taus[i]
  idxlo = taunow - sigtaulim
  idxhi = taunow + sigtaulim
  idnow = np.where((delsave > idxlo) & (delsave < idxhi))[0]
  
  a = (taunow - delsave[idnow])/sigtau
  a2 = -a*a/2
  ea2 = np.exp(a2)
  #ea2sum = np.sum(ea2) dont seem to need ea2sum but put back in later if necessary
  ea2out= np.sum(wsave[idnow]*ea2)
  psis[i] = ea2out#/ea2sum
  

 
 psis = np.nan_to_num(psis,0)  

 if (norm == 1):
  pt   = (psis - np.min(psis))
  psis = pt/np.max(pt)
 

 
 if (diagnose == 1):
   return(psis,delsave,wsave,rsave,azsave)
 else:
  return(psis)
#






#outputs are psis[ntau] to input set of time grid[ntau]
#make code to test and plot make nice plots
# 
#import timeit 
#
#wavref = 4000.
#embhref   = 1.e8
#emdotref  = 1.0
#degref    = 0.0
#
#t0v   = -1.6e4
#t0i   = -1.2e4
#taus   = np.arange(0.0,30.,0.1)
#ntau = np.shape(taus)[0]
#
#
#
#
#
#
#color=['k','r','b','cyan','purple','orange','gree','skyplue','magenta']
#
#deg = [0.0,20.0,40.0,60.0,80.0]
#fig = plt.figure(111)
#ax1 = fig.add_subplot(111)
#idx = 0
#for dnow in deg:
# psis   = pytfb_sub(taus,embhref,emdotref,wavref, dnow,t0vin = t0v, t0iin = t0i, norm = 1,quick = 1,oldsmooth=0,newsmooth=1,diagnose = 0)
# ax1.plot(taus,psis,color=color[idx],label='deg = '+np.str(dnow))
# #calculate the mean
# mean = np.sum(taus*psis)/np.sum(psis)
# ax1.plot([mean,mean],[0.0,1.0],color=color[idx],label=None)
# idx = idx + 1
#ax1.set_xlabel('delay (days)')
#ax1.set_ylabel('response (arbitrary units)')
#plt.savefig('plot_response_inc.pdf')
#











#other diagnostic tests testing smoothing function 



#start = timeit.default_timer()
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 1)
#
#ds = psis[1]
#ws = psis[-1]
#
#ida = np.argsort(ds)
#ds = ds[ida]
#ws = ws[ida]
#nws = len(ws)
#
##weights = np.zeros((ntau,nws))
#
##tauds = np.zeros((ntau,nws))
##tauds[:,:] = taus
##tauds[:,:] = taus - ds
#
##dssub = ds - taus
##weights = np.exp(ds - taus)
#
#
#psis = psis[0]
#stop = timeit.default_timer()
#
#
#
#print 'run time', stop-start
# 
#fig = plt.figure()
#ax1 = fig.add_subplot(111)

#ax1.plot(ds,ws/ws.max())


#start = timeit.default_timer()
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 0)
#
#ds = psis[1]
#ws = psis[-1]
#
#ida = np.argsort(ds)
#ds = ds[ida]
#ws = ws[ida]
#nws = len(ws)
#
##weights = np.zeros((ntau,nws))
#
##tauds = np.zeros((ntau,nws))
##tauds[:,:] = taus
##tauds[:,:] = taus - ds
#
##dssub = ds - taus
##weights = np.exp(ds - taus)
#
#
#psis = psis[0]
#stop = timeit.default_timer()
#print 'run time', stop-start
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.plot(taus,psis,label='quick off old smooth')
#ax1.plot(ds,ws/ws.max(),label='quick off no smooth')
#










#experiment 2 
#
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 0)
#ds = psis[1]
#ws = psis[-1]
#ida = np.argsort(ds)
#ds = ds[ida]
#ws = ws[ida]
#nws = len(ws)
#
#datnew = ss.binned_statistic(ds,ws,statistic='mean',bins=taus)
#a = datnew[0]
#a = np.insert(a,0,0)
#ax1.plot(taus,a)
#
#plt.savefig('pytfbplottest_quick.pdf')

#scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
#try making histogram of points and sampling mean from each bin



















#experiment 3 
#
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 0)
#ds = psis[1]
#ws = psis[-1]
#ida = np.argsort(ds)
#ds = np.array(ds[ida])
#ws = np.array(ws[ida])
#nws = len(ws)
#dtau = taus[1] - taus[0]
#sigtau = 5.*dtau
#ntaus = np.shape(taus)[0]
#psinew = np.zeros(ntaus)
#
##taunew = np.insert(np.arange(-10.0
#for i in range(ntaus):
# taunow = taus[i]
# a = (taunow - ds)/sigtau
# a2 = -a*a/2
# ea2 = np.exp(a2)
# ea2sum = np.sum(ea2)
# ea2sum = np.sum(ws*ea2/ea2sum)
# psinew[i] = ea2sum
# 
#ax1.plot(taus,psinew/psinew.max(),label='quick off new smooth')
#ax1.legend()
#plt.savefig('pytfbplottest_quick.pdf')
#
#scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
#try making histogram of points and sampling mean from each bin






##experiment 3
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#
#print 'experiment 1 old smoothing quick = 1'
#start = timeit.default_timer()
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 1,oldsmooth=1,newsmooth=0,diagnose = 1)
#stop = timeit.default_timer()
#print 'run time', stop-start
#ax1.plot(taus,psis[0]/psis[0].max(),label='experiment 1 old smoothing quick = 1')
#
#print 'experiment 2 no smoothing quick = 1'
#start = timeit.default_timer()
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 1,oldsmooth=0,newsmooth=0, diagnose = 1)
#stop = timeit.default_timer()
#print 'run time', stop-start
#ax1.plot(psis[1],psis[2]/psis[2].max(),label='experiment 2 no smoothing quick = 1')
#
#
#
#print 'experiment 3 new smoothing quick = 1'
#start = timeit.default_timer()
#psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1,quick = 1,oldsmooth=0,newsmooth=1,diagnose = 0)
#stop = timeit.default_timer()
#print 'run time', stop-start
#ax1.plot(taus,psis/psis.max(),label='experiment 3 new smoothing quick = 1')
#
#ax1.legend()
#plt.savefig('pytfbplottest_quick.pdf')
#
#
#
#
###plot the azimuth and radius grid to identify the source of hte hump
#rsave = psis[3]
#azsave = psis[4]
#fig = plt.figure()
#ax1 = fig.add_subplot(311)
#ax1.hist(rsave,bins = 200)
#ax1 = fig.add_subplot(312)
#ax1.hist(azsave,bins = 200)
#ax1 = fig.add_subplot(313)
#ax1.hist(psis[1],bins = 200)
##plt.show()
#plt.savefig('pytfbplottest_hist.pdf')
#
#
##fig = plt.figure()
##ax1 = fig.add_subplot(111)
##ax1.scatter(rsave*np.cos(azsave),rsave*np.sin(azsave))
##plt.savefig('pytfbplottest_scatter.pdf')
#












