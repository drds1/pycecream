import numpy as np
import astropy


#function to perform running sigma clip on inpu light curve

def myrunsigclip(t,x,width,sd=4,iters=None):

 xchop = np.zeros((0))
 tchop = np.zeros((0))
 tmin  = t.min()
 tmax  = t.max()
 
 tlo = 1.*tmin
 thi = tlo + width
 tloop=[tlo]
 while (thi <= tmax):
  idx = (t > tlo) & (t <= thi)
  tnow = t[idx]
  xnow = x[idx]
  
  xnow_clip = astropy.stats.funcs.sigma_clip(xnow, sig=sd, iters=iters)
  idxchop = xnow_clip[1]
  xnow_clip = xnow_clip[0]

  
  tchop=np.append(tchop,tnow[idxchop])
  xchop=np.append(xchop,xnow_clip)
  tlo = tlo + width
  thi = thi + width
  tloop.append(tlo)
 return(tchop,xchop,tloop)
  
  
#function to perform running sigma clip on input light curve (no segments re-evauate for each point)

def myrunsigclip2(t,x,width,sd=4,iters=None):
 fracl=[]
 nx = x.shape[0]
 xtally = np.zeros(nx)
 xpicked = np.zeros(nx)
 idx = np.arange(nx)

 xchop = np.zeros((0))
 tchop = np.zeros((0))
 tmin  = t.min()
 tmax  = t.max()
 
 ic = 0
 for tn in t:
 
  tlo = max(tmin,tn-width/2)
  thi = min(tmax,tn+width/2)

  tloop=[tlo]
  #idxnow = idx[(t > tlo) & (t <= thi)]
  idxnow = idx[(t > tlo) & (t <= thi) & (t != tn)]
  xpicked[idxnow] = xpicked[idxnow] + 1#[xa + 1 for xa in xpicked[idxnow]]
  tnow = t[idxnow]
  xnow = x[idxnow]
  
  xnow_clip = astropy.stats.funcs.sigma_clip(xnow, sig=sd, iters=iters)
  idxchop = idxnow[xnow_clip[1]]#xnow_clip[1]
  idxn    = idxnow[np.where(tnow == tn)[0]]
  xnow_clip = xnow_clip[0]
  xtally[idxchop]= xtally[idxchop] + 1
  #xtally[idxn] = xtally[idxn] - 1 #do not count the current value as a success
  
  #tests
  #nxn = np.shape(idxnow)[0]
  #nxc = np.shape(idxchop)[0]
  #print idxnow, nxn
  #print idxchop,nxc
  #print xtally[idxnow],idxn
  #raw_input()
  
  #for i in range(nxn):
  # print i, idxnow
  
  ic = ic + 1
  frac = int(100.*ic/nx)
  if (np.mod(frac,10) == 0):
   test = frac in fracl
   if (test == 0):
    print 'myrunsigclip2',frac,'pc complete'
    fracl.append(frac)

   
 #idxreturn = idx[np.where(np.array(xtally > 0))[0]]
 idxreturn = idx[np.where(np.array(xtally == xpicked))[0]]
 tchop =  t[idxreturn]
 xchop = x[idxreturn]
 
 
 return(tchop,xchop,xpicked,xtally)
 

 