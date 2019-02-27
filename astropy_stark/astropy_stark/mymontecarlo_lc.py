import numpy as np



def mymc_lc(taunow,lc,sig,nits = 100,means=0):

 
 ndat = np.shape(lc)[0]
 lcsave = np.zeros((ndat,nits))
 
 taumeansave = []
 for it in range(nits):
  lcnow = lc + np.random.randn(ndat)*sig
  lcsave[:,it] = lcnow
  if (means == 1):
   taumean = np.sum(lcnow*taunow)/np.sum(lcnow)
   taumeansave.append(taumean)
   
 lcout = np.array(np.percentile(lcsave,[15.865,50,84.315],axis = 1)).T
 if (means == 0):
  return(lcout)
 else:
  taumeansave = np.array(taumeansave)
  means_out = np.percentile(taumeansave,[15.865,50,84.315])
  return(lcout,means_out)