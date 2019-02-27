# python function takes an input array calculates the mean, lower and upper quartiles
# ip: parm[:nparm], conlev the confidence level required i.e 0.68 for 1 sigma etc
# op: the median, lower and upper confidence regions
import numpy as np

def mymedquart(parm,conlev=0.68):
 
 parmsort = np.sort(parm)
 nparm = parmsort.shape[0]
 
 
 idxparmmed = nparm/2

# calculate the median  
 if (np.mod(nparm,2) ==0):
  parmmed = (parmsort[idxparmmed-1] + parmsort[idxparmmed])/2 
 else:
  parmmed = parmsort[idxparmmed]
#



 a = conlev*nparm/2
 idxlo = int(np.max((idxparmmed - a,0)))
 idxhi = int(np.min((idxparmmed + a,nparm-1)))
 
 return(parmmed,parmsort[idxlo], parmsort[idxhi])


# test below on 5/2/2014. seems to work
#dat = np.random.randn(1000)
#
#op = mymedquart(dat,0.68)
#
#from pylab import *
#
#hist(dat)
#
#vlines(op[1],0,1000)
#vlines(op[0], 0, 1000)
#vlines(op[2], 0, 1000)
#
#show()



