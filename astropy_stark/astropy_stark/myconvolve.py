import numpy as np
#convolve xd by tf 1d


def myconvolve(xd,tf):

 nmodtf = np.shape(tf)[0]
 nd     = np.shape(xd)[0]
 #manual convolution (i dont know how to use the python one!)
 modf= []
 for i in range(nd):
  xnow = xd[i]
  idxlo = max(0,i-nmodtf)
  idxhi = i
  idxd = np.arange(idxlo,idxhi,1)[::-1]
  idxtf = np.arange(idxhi  - idxlo)
  try:
   modf.append( np.sum(xd[idxd]*tf[idxtf]) )
  except:
   modf.append(0)
 modf = np.array(modf)
 return(modf)
 
 

def mc3(td,xd,tau,tf):



 nmodtf = np.shape(tf)[0]
 nd     = np.shape(xd)[0]
 #manual convolution (i dont know how to use the python one!)
 modf= []
 
 taulo = np.min(tau)
 tauhi = np.max(tau)
 dtau  = np.mean(tau[1:]-tau[:-1])
 ntau  = np.shape(tau)[0]

 dtgrid = np.mean(td[1:]-td[:-1])
 #check if driver and tau grid are on same scale, interpolate tf if not
 
 if (dtau != dtgrid):
  #print('myconvolve response and driver not on same time grid, interpolating response function.',dtgrid,dtau)
  tauint = np.arange(taulo,tauhi+dtgrid,dtgrid)
  tfint  = np.interp(tauint,tau,tf) 
 else:
  tauint = tau
  tfint = tf
 
 itflo = int(np.floor(taulo/dtgrid))
 itfhi = int(np.floor(tauhi/dtgrid))
 
 
 for i in range(nd):
  xnow = xd[i]
  idxlo = max(0,i-itfhi)
  idxhi = min(nd-1,i-itflo)
  idxd = np.arange(idxlo,idxhi,1)[::-1]
  idxtf = np.arange(idxhi  - idxlo)
  #print np.shape(idxtf),np.shape(idxd),np.sum(xd[idxd]*tf[idxtf])
  try:
   modf.append( np.sum(xd[idxd]*tfint[idxtf]) )
  except:
   modf.append(0)
 modf = np.array(modf)
 #print modf
 #raw_input()
 return(modf)