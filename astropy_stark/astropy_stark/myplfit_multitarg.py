
#fit a power law to multiple target observations
#such that we fit y = y0_{1...N} x/x0^slope 
#i.e one slope but N normalisation offsets
#useful for multi object accretion disk lag satudies


import numpy as np



#inputs
#x[Ntarg,Ndim]
#y[Ntarg,Ndim]
#sigy[Ntag,Ndum]


#outputs
#os[Ntarg],sigos[Ntarg] the normalisations and variance
#slope, sigslope the slope and variance

def mypl_multi(x,y,sig,nits = 100):

 #now do analysis
 ntarg,ndim = np.shape(x)
   
 xlog    = np.log10(x)
 ylog    = np.log10(y)
 sigylog = sig/lag/np.log(10)
 wavl0 = np.sum(xlog/sigylog**2)/np.sum(1./sigylog**2)
 
 sy2 = sigylog**2
 
 nits = 100
 #use iterated optimal scaling to fit multi offset line
 slope = 1.0#4./3
 os = np.ones(ntarg)
 varos = np.ones(ntarg)
 
 xlwl0 = xlog-wavl0
 for it in range(nits):
  
  #update os
  for ip in range(ntarg):
   bot = np.sum(1./sy2[ip,:])
   osnew = np.sum((ylog[ip,:] - slope*xlwl0[ip,:])/sy2[ip,:])/bot
   varsnew = 1./bot
   os[ip] = osnew
   varos[ip] = varsnew
 
  ost = np.tile(os,(ndim,1)).T
  #update slope
  top = np.sum((ylog-ost)*xlwl0/sy2)
  bot = np.sum(xlwl0**2/sy2)
  slope = top/bot
  varslope = 1./bot
 
 
 
 return(os,varos,slope,varslope) 
 
 
 
