import numpy as np
import matplotlib.pylab as plt
import mytemp0 as mt0

twopi = np.pi*2
deg2rad = np.pi/180
planck = 6.626307004e-34
c      = 2.99792458e8
boltz  = 1.38064852e-23




def pytfb_sub(taus,embh,emdot,wavang, deginc, t0vin=-1, t0iin = -1, alpha_visc = -0.75, hxsch = 3.0, alpha_irad = -0.75, eta = 0.1, rlosch = 3.0, norm = 1):
 tauhi  = taus[-1]
 ntaus  = np.shape(taus)[0]
 psis   = np.zeros(ntaus)
 
 
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
 naz    = 100
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
 
 
 nazsub1 = naz - 1
 
 #now define radius grid
 rgrid  = np.logspace(np.log10(rlold),np.log10(rhild),nrad)
 azgrid = np.linspace(0.0, twopi,naz)
 azgrid_s = azgrid[1:] - azgrid[:-1]
 #calculate temperature at each radius grid
 av4   = alpha_visc*4 
 ai4   = alpha_irad*4
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
 azwidth = twopi/nazsub1
 azsave=[]
 razsave=[]
 for i in range(nrad-1):
  
  #quick way 
  #tv4now    = tv4[i]
  #ti4now    = ti4[i]
  #ttotnow   = ttot[i]
  #ttot3now  = ttot3[i]
  
  
  radlo = rgrid[i]
  radhi = rgrid[i+1]
  drad  = radhi - radlo
  raz   = np.random.uniform(radlo,radhi,nazsub1)
  
  tv4   = t0v4*(raz)**av4 #modification for inner radius to the right (negligible)* (1 - np.sqrt(rlold/raz)) / (1 - np.sqrt(rlold/1))
  ti4   = t0i4*(raz)**ai4
  rir0b = ti4/t0i4 #this is just (r/r0)^alpha_irad
  ttot4 = tv4 + ti4
  ttot2 = np.sqrt(ttot4) 
  ttot  = np.sqrt(ttot2)
  ttot3 = ttot2*ttot
 
  
  daz  = np.sqrt(raz*raz + hx2)
  
  
  az   = np.random.uniform(nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
  caz  = np.cos(az)
  tdl  = hxci - raz*caz*sini + daz
  

  
  x    = hc_kwav/ttot#hc_kwav/ttotnow
  x2   = x*x
  
  #the radlo * drad *azwidth is the solid angle element
  #azwidth can be left off here as it is always the same
  weight = rir0b/ttot3 * x2/wavang2/(np.cosh(x) - 1) * radlo*drad*azwidth      #rir0b[i]/ttot3now * x2/wavang2/(np.cosh(x) - 1)
  wsave.append(weight)
  delsave.append(tdl)
  
  azsave.append(az)
  razsave.append(raz)
  #introduce a lower threshold weight to abort the iterations (to save time)
  
 
 
 delsave = [item for sublist in delsave for item in sublist]
 wsave   = [item for sublist in wsave for item in sublist]
 azsave = [item for sublist in azsave for item in sublist]
 razsave   = [item for sublist in razsave for item in sublist]
 
 
 razsave = np.array(razsave)
 azsave = np.array(azsave)
 delsave = np.array(delsave)
 wsave   = np.array(wsave)
 nsave =np.shape(delsave)[0]
 dat = np.zeros((nsave,3))
 dat[:,0] = razsave
 dat[:,1] = azsave
 dat[:,2] = delsave
 np.savetxt('shittyprob.txt',dat)
 
 
 
 #save the results
 
 
 
 
 #add smoothing function to approximate delta function
 #ntaus = np.shape(taus)[0]
 sigtau = taus[1] - taus[0]
 sigtaulim = 2.5*sigtau
 
 
 for i in range(ntaus-1):
  taulo = taus[i]
  tauhi = taus[i+1]
  taulimlo = taulo - sigtaulim
  taulimhi = tauhi + sigtaulim
  taucent = (taulo + tauhi)/2
  idxsmooth = np.where((delsave >= taulimlo) & (delsave < taulimhi))[0]
  taunow    = delsave[idxsmooth]
  wnow      = wsave[idxsmooth]
  
  a = taunow - taucent
  b = sigtau
  aa_bb2 = a*a/b/b
  gauss     = np.exp(-aa_bb2/2)
  gsum      = np.sum(gauss)
 

  #WONT SOME SURFACE ELEMENTS BE COUNTED TWICE AND SO BE INCLUDED IN MORE THAN ONE GAUSSIAN WEIGHTING? YES BUT I DONT THINK THIS MATTERS
  psis[i] = np.sum( wnow*gauss/gsum )
  
  
  #print taus[i], psis[i], ntaus
 #if norm == 1then normalise so max == 1
 if (norm == 1):
  pt   = (psis - np.min(psis))
  psis = pt/np.max(pt)
 return(psis)







#outputs are psis[ntau] to input set of time grid[ntau]
#make code to test and plot 
 
import timeit 
wavang = 4000.
embh   = 1.e8
emdot  = 1.0
deginc = 0.0

t0v   = -1.6e4
t0i   = -1.2e4
taus   = np.arange(0.0,15.,0.05)



start = timeit.default_timer()
psis   = pytfb_sub(taus,embh,emdot,wavang, deginc,t0vin = t0v, t0iin = t0i, norm = 1)
stop = timeit.default_timer()



print 'run time', stop-start
 
fig = plt.figure()
ax1 = fig.add_subplot(111)


ax1.plot(taus,psis)
plt.savefig('pytfbplottest_2.pdf')
















