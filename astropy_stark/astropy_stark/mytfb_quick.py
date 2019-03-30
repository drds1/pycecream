import numpy as np
import matplotlib.pylab as plt
import astropy_stark.mytemp0 as mt0
import scipy.stats as ss

twopi = np.pi*2
deg2rad = np.pi/180
planck = 6.626307004e-34
c      = 2.99792458e8
boltz  = 1.38064852e-23




def pytfb_sub(taus,embh,emdot,wavang, deginc, t0vin=-1, t0iin = -1, alpha_visc = -0.75, hxsch = 3.0, alpha_irad = -0.75, eta = 0.1, rlosch = 3.0, norm = 1, quick = 1, xstop = 15, udlnr = 0.01,thcent=1.0,thfwhm=0.2,oldsmooth = 0,newsmooth = 1,diagnose=0):
 tauhi  = taus[-1]
 taulo  = taus[0]
 ntaus  = np.shape(taus)[0]
 psis   = np.zeros(ntaus)
 
 #if you want a top hat response then its easy else do what you had before
 if (wavang < 0.0):
  idxinc = np.where((taus > thcent - thfwhm/2) & (taus < thcent + thfwhm/2))[0]
  psis[idxinc] = 1
 
 
 else:
 #either input desired reference temperature t0 at 1 light day, or calculate the
 #value appropriate for black body accretion disc given m and mdot 
  if (t0vin < 0):
   t0v   = mt0.tv0(embh,emdot)
  else:
   t0v   = t0vin
  
  if (t0iin < 0):
   t0i   = mt0.ti0(embh,emdot, eta = eta)
  else:
   t0i   = t0iin 
   
 
  #rlosch = 3.0
  #wavang = 4000.
  #hxsch  = 3.0
  
  
  nrad   = 10000
  
  #need to calculate 1/T^3 (r/r0)^alpha_irad x^2/wav^2 / (cosh(x) - 1) delta (tau - tau(r,theta)) dtau
  
  
  t0v2   = t0v*t0v
  t0v4   = t0v2*t0v2
  
  t0i2   = t0i*t0i
  t0i4   = t0i2*t0i2
  
  rsch   = 1.15821e-10*embh #scwarzchild radius in light days
  hx     = hxsch*rsch
  hx2    = hx*hx
  cosi   = np.cos(deginc*deg2rad)
  sini   = np.sin(deginc*deg2rad)
  hxci = hx*cosi
  hc_kwav = planck*c/wavang/1.e-10/boltz
  wavang2 = wavang*wavang
  rlold  = rlosch*rsch
  
  
  #print t0v,t0i
  
  #this estimate for the cutoff radius is based on the max radius of the
  #highest lag (re-arrange equation 3 in Starkey et al 2017)
  #should be ok but might exclude some significant low lags for a VERY hot, edge on disk
  rhisch = -1.* max((tauhi - hxci)/(1.+sini), tauhi) #if -ve then this is in light days
  #print rhisch
  if (rhisch < 0):
   rhild  = np.abs(rhisch)
  else:
   rhild  = rhisch*rsch
  
  
  
  
  
  av4   = alpha_visc*4 
  ai4   = alpha_irad*4 
  #now define radius grid be smart here.
  
  #use a cutoff x_stop to determine when to stop the radius grid
  dtau   = taus[1] - taus[0]
  rhilog = dtau
  rtemp = np.array([rhilog,10*rhilog]) 
  rtl   = np.log(rtemp)
  ttemp4 = t0v4*(rtemp)**av4 + t0i4*(rtemp)**ai4
  ttemp = np.sqrt(np.sqrt(ttemp4))
  ttl    = np.log(ttemp) 
  tstop  = hc_kwav/xstop
  grad = (rtl[1] - rtl[0]) / (ttl[1] - ttl[0])
  rhil = rtl[0] +  grad*( np.log(tstop) - ttl[0] )
  rhild   = np.exp(rhil)
  
  #rt4   = rtemp**4
  #ttemp4 = t0v4*(rtemp)**av4 + t0i4*(rtemp)**ai4
  #grad = (rt4[1] - rt4[0])/(ttemp4[1] - ttemp4[0])
  #t_stop4 = (hc_kwav/xstop)**4
  #rhild4 = rt4[0] + grad*(t_stop4 - ttemp4[0])
  #ttemp=ttemp4**0.25
  #xtemp = hc_kwav/ttemp
  #rhild4 = np.interp(t_stop4,ttemp4,rt4)
  #rhild = np.sqrt(np.sqrt(rhild4))
  #print rtemp
  #print ttemp
  #print xtemp
  #print xstop, tstop, rhild, grad, rtl, ttl
  #print rhild
  #print 'now defined a stopping radius e.g radius at which'
  #print 'temp give x =hc/(wav T) > xstop. xstop default is 15 now set up radius grid with this'
  #print 'upper radius BUT how can it be 2 light days because we need response at time delays',
  #print 'longer than this? FIX THIS' 
  #raw_input()
  #build grid
  rgridlog = np.exp( np.arange(np.log(rlold), np.log(rhilog), udlnr)[:-1] )#np.logspace(np.log(rlold),np.log(rhilog))
  rgridlin = np.arange(rhilog,rhild,rhilog)
  rgrid    = np.concatenate((rgridlog,rgridlin))
  nrad     = np.shape(rgrid)[0]
  dr       = rgrid[1:] - rgrid[:-1]
  
  
  
  
  #rgrid  = np.logspace(np.log10(rlold),np.log10(rhild),nrad)
  
 
  
  
  
  
  #calculate temperature at each radius grid
  tv4   = t0v4*(rgrid)**av4
  ti4   = t0i4*(rgrid)**ai4
  rir0b = ti4/t0i4 #this is just (r/r0)^alpha_irad
  ttot4 = tv4 + ti4
  ttot2 = np.sqrt(ttot4) 
  ttot  = np.sqrt(ttot2)
  ttot3 = ttot2*ttot
  
  #each point will have a time delay and a weighting append these to a list
  #for each point in the disc
  
  #loop of azimuths
  delsave = []
  wsave   = []
  if (diagnose == 1):
   azsave = []
   rsave = []
  
  if (quick == 1):
   for i in range(nrad-1):
    #quick way 
    ttotnow   = ttot[i]
    ttot3now  = ttot3[i]
    radlo = rgrid[i]
    radhi = rgrid[i+1]
    drad  = dr[i]#radhi - radlo
    
    #now azimuth grid
    azwidth = drad/radlo
    azgrid = np.arange(0.0, twopi,azwidth)
    naz = np.shape(azgrid)[0]
    nazsub1 = naz - 1
    azgrid_s = azgrid[1:] - azgrid[:-1]
    
    #print drad, radhi, azwidth,'sdfsd',radhi*azwidth
    #print radlo, radhi, azwidth, drad, naz,'sdfsf'
    
    raz   = np.random.uniform(radlo,radhi,nazsub1)   
    daz  = np.sqrt(raz*raz + hx2)
    az   = np.random.uniform(low=azgrid[:-1],high=azgrid[1:],size=nazsub1)#np.random.uniform(low=0,high=1,size=nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
    caz  = np.cos(az)
    tdl  = hxci - raz*caz*sini + daz
    x    = hc_kwav/ttotnow#hc_kwav/ttot#hc_kwav/ttotnow
    
     
    #print radlo, x
    x2   = x*x
    #the radlo * drad *azwidth is the solid angle element
    #azwidth can be left off here as it is always the same
    weight = rir0b[i]/ttot3now * x2/wavang2/(np.cosh(x) - 1) * radlo*drad*azwidth                 #rir0b/ttot3 * x2/wavang2/(np.cosh(x) - 1) * radlo*drad*azwidth      #
    wsave.append([weight]*nazsub1)
    delsave.append(tdl)
    if (diagnose == 1):
     azsave.append(az)
     rsave.append(raz)
    
  else:
   for i in range(nrad-1):
    radlo = rgrid[i]
    radhi = rgrid[i+1]
    drad  = radhi - radlo
    
    #now azimuth grid
    azwidth = drad/radhi
    azgrid = np.arange(0.0, twopi,azwidth)
    naz = np.shape(azgrid)[0]
    nazsub1 = naz - 1
    azgrid_s = azgrid[1:] - azgrid[:-1]
    
    raz   = np.random.uniform(radlo,radhi,nazsub1)
    daz  = np.sqrt(raz*raz + hx2)
    az   = np.random.uniform(low=azgrid[:-1],high=azgrid[1:],size=nazsub1)#np.random.uniform(azgrid[:-1],azgrid[1:],1)
    caz  = np.cos(az)
    tdl  = hxci - raz*caz*sini + daz
    tv4   = t0v4*(raz)**av4 #modification for inner radius to the right (negligible)* (1 - np.sqrt(rlold/raz)) / (1 - np.sqrt(rlold/1))
    ti4   = t0i4*(raz)**ai4
    rir0b = ti4/t0i4 #this is just (r/r0)^alpha_irad
    ttot4 = tv4 + ti4
    ttot2 = np.sqrt(ttot4) 
    ttot  = np.sqrt(ttot2)
    ttot3 = ttot2*ttot
    x      = hc_kwav/ttot
    x2     = x*x
    weight = rir0b/ttot3 * x2/wavang2/(np.cosh(x) - 1) * radlo*drad*azwidth
    wsave.append(weight)
    delsave.append(tdl)
    if (diagnose == 1):
     azsave.append(az)
     rsave.append(raz)
    #print 'az',az
    #print ''
    #print 'azgrid',azgrid
   #tv4   = t0v4*(raz)**av4 #modification for inner radius to the right (negligible)* (1 - np.sqrt(rlold/raz)) / (1 - np.sqrt(rlold/1))
   #ti4   = t0i4*(raz)**ai4
   #rir0b = ti4/t0i4 #this is just (r/r0)^alpha_irad
   #ttot4 = tv4 + ti4
   #ttot2 = np.sqrt(ttot4) 
   #ttot  = np.sqrt(ttot2)
   #ttot3 = ttot2*ttot
   
   
   
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
  if (oldsmooth == 1):
   sigtau2 = sigtau*sigtau
   
   #print('ok now it works but is slow')
   #print('need to find a numpy or scipy gaussian smoothing function')
   #print('ready made function wil be faster than anything I write down below')
   #print('was fine for fortran but python prefers ready made smoothing function')
   #print('google numpy or scipy gaussian filter 1d etc')
   #raw_input()
   #psis_2 = scipy.ndimage.filters.gaussian_filter1d(, sigtau, axis=-1, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0) 
   psis = np.zeros(ntaus)
   for id in range(nds):
    delnow = delsave[id]
    idlo = int( max(1,np.floor((delnow - sigtaulim - taulo)/dtau)) )  
    #print 'dfsdf',ntaus, np.ceil((delnow + sigtaulim)/dtau), np.max(ntaus,np.ceil((delnow + sigtaulim)/dtau))
    idhi = int( min(ntaus,np.ceil((delnow + sigtaulim - taulo)/dtau)) )
    idinc = np.arange(idlo,idhi,1)
    #print id,delnow,idlo,idhi,taus[idlo],taus[idhi]
    tdelsubtau = delnow - taus[idinc]
    gtemp = np.exp(-0.5*tdelsubtau*tdelsubtau/sigtau2)
    gsum = np.sum(gtemp)
    wnow  = wsave[id]
    psis[idinc] = psis[idinc] + wnow*gtemp/gsum 
  
  
  elif (newsmooth ==1):  
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
   
  #if no smoothing just sort the uneven array into ascending time order and release it (do not use psis as output in this case. BEST NOT TO USE THIS SETTING) 
  else:
   ida = np.argsort(delsave)
   delsave = delsave[ida]
   wsave = wsave[ida]
   if (diagnose == 1):
    azsave = azsave[ida]
    rsave = rsave[ida]   
  
  psis = np.nan_to_num((psis,0)) 
  psis = psis[0] 
 # for i in range(ntaus-1):
 #  taulo = taus[i]
 #  tauhi = taus[i+1]
 #  taulimlo = taulo - sigtaulim
 #  taulimhi = tauhi + sigtaulim
 #  taucent = (taulo + tauhi)/2
 #  idxsmooth = np.where((delsave >= taulimlo) & (delsave < taulimhi))[0]
 #  taunow    = delsave[idxsmooth]
 #  wnow      = wsave[idxsmooth]
 #  
 #  a = taunow - taucent
 #  b = sigtau
 #  aa_bb2 = a*a/b/b
 #  gauss     = np.exp(-aa_bb2/2)
 #  gsum      = np.sum(gauss)
 # 
 #
 #  #WONT SOME SURFACE ELEMENTS BE COUNTED TWICE AND SO BE INCLUDED IN MORE THAN ONE GAUSSIAN WEIGHTING? YES BUT I DONT THINK THIS MATTERS
 #  psis[i] = np.sum( wnow*gauss/gsum )
 #  
 #  
 #  #print taus[i], psis[i], ntaus
  #if norm == 1then normalise so max == 1
  if (norm == 1):
   #print psis[0]
   pt   = (psis - np.min(psis))
   psis = pt/np.max(pt)
 

 
 if (diagnose == 1):
   return(psis,delsave,wsave,rsave,azsave)
 else:
  return(psis)
#






#outputs are psis[ntau] to input set of time grid[ntau]
#make code to test and plot make nice plots
## 
#import timeit 
#
#wavref = 4000.
#embhref   = 1.e8
#emdotref  = 1.0
#degref    = 0.0
##
#t0v   = -1.6e4
#t0i   = -1.2e4
#taus   = np.arange(0.0,30.,0.1)
#ntau = np.shape(taus)[0]
##
##
##
##
##
##
#color=['r','b','cyan','purple','orange','gree','skyplue','magenta']
#fortrancomp = '/Users/ds207/Documents/standrews/sta/fort/fortcode/delaydist/tfbtest.dat'
#
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
#if (fortrancomp != ''):
# fc = np.loadtxt(fortrancomp)
# ax1.plot(fc[:,0],fc[:,1],linewidth=2,color='k',label='fortran comp')
# mean = np.sum(fc[:,0]*fc[:,1])/np.sum(fc[:,1])
# ax1.plot([mean,mean],[0.0,1.0],color='k',linewidth=2,label=None)
#ax1.set_xlim([0,20.])
#ax1.legend()
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












