import numpy as np
from mybnu import *

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!calculate the half light radius at a given wavelength!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def rhalflam(wav,temp,rld):
 #print wav,'halflight'

 ntemp = np.shape(temp)[0]
 #print wav, 'just before'

 #try:
 bnu2  = np.vectorize(bnu)
 try:
  ilo = np.where(temp == 0)[0][-1]+1
 except:
  ilo = 0
 ilo=0
 plank = np.double(bnu2(np.double(np.array(temp[ilo:])),wav))
 #print plank, 'myhalflight', np.double(bnu2(np.array(temp[ilo:]),wav)),'ghghg',ilo
 #print 'checking', ilo, temp[ilo-1],temp[ilo:ilo+20]
 
 #raw_input()
 #except:
  #plank = bnu(temp[1:],wav)

 dr = rld[1:] - rld[:-1]
 rldnew = rld[1:]
 
 
 sa = 2*np.pi*rldnew*dr
 ar = sa*plank
 
 ltot = np.sum(ar)
 lcumsum = np.cumsum(ar)
 
 try:
  idxhalf = np.where(lcumsum >= ltot/2)[0][0]
 except:
  idxhalf = np.where(lcumsum >= ltot/2)[0]
  
 rhalf = rld[idxhalf]
 
 return(rhalf)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!


