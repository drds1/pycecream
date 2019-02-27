import numpy as np
import matplotlib.pylab as plt
import mytemp0 as mt0

twopi = np.pi*2
deg2rad = np.pi/180
planck = 6.626307004e-34
c      = 2.99792458e8
boltz  = 1.38064852e-23




def pytfb_sub(taus,embh,emdot,wavang, deginc, t0vin=-1, t0iin = -1, alpha_visc = -0.75, hxsch = 3.0, alpha_irad = -0.75, eta = 0.1, rlosch = 3.0, norm = 1, quick = 1, xstop = 15, udlnr = 0.01,thcent=1.0,thfwhm=0.2):
 tauhi  = taus[-1]
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
    az   = np.random.uniform(nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
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
    az   = np.random.uniform(nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
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
  delsave = np.array(delsave)
  wsave   = np.array(wsave)
  
  nds = np.shape(delsave)[0]
  
  
  
  #add smoothing function to approximate delta function
  #ntaus = np.shape(taus)[0]
  sigtau = 5*dtau/2
  sigtau2 = sigtau*sigtau
  sigtaulim = 5*sigtau
  
  
  print 'ok now it works but is slow'
  print 'need to find a numpy or scipy gaussian smoothing function'
  print 'ready made function wil be faster than anything I write down below'
  print 'was fine for fortran but python prefers ready made smoothing function'
  print 'google numpy or scipy gaussian filter 1d etc'
  #raw_input()
  #psis = scipy.ndimage.filters.gaussian_filter1d(, sigma, axis=-1, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)
  psis = np.zeros(ntaus)
  for id in range(nds):
   delnow = delsave[id]
   idlo = int( max(1,np.floor((delnow - sigtaulim)/dtau)) )
   
   #print 'dfsdf',ntaus, np.ceil((delnow + sigtaulim)/dtau), np.max(ntaus,np.ceil((delnow + sigtaulim)/dtau))
   idhi = int( min(ntaus,np.ceil((delnow + sigtaulim)/dtau)) )
   idinc = np.arange(idlo,idhi,1)
   #print idinc
   
   tdelsubtau = delnow - taus[idinc]
   gtemp = np.exp(-0.5*tdelsubtau*tdelsubtau/sigtau2)
   gsum = np.sum(gtemp)
   wnow  = wsave[id]
   psis[idinc] = psis[idinc] + wnow*gtemp/gsum 
  
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
   pt   = (psis - np.min(psis))
   psis = pt/np.max(pt)
 
 return(psis)
#






#outputs are psis[ntau] to input set of time grid[ntau]
#make code to test and plot 
 
import timeit 
wavang = 4000.
embh   = 1.e8
emdot  = 1.0
deginc = 0.0

t0v   = -1.6e4
t0i   = -1.2e4
taus   = np.arange(0.0,30.,0.1)



start = timeit.default_timer()
psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1)
stop = timeit.default_timer()



print 'run time', stop-start
 
fig = plt.figure()
ax1 = fig.add_subplot(111)


ax1.plot(taus,psis)
plt.savefig('pytfbplottest_quick.pdf')
















