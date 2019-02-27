import numpy as np
import myreadlines

#define an array to load in parameters and reject certain values, returning the masked
#points in an array idxmask

def myloadpar(fname,rejectinf==1,rejectnan==1,rejectother==[0]):
 try:
  dat = np.loadtxt(fname)
 except:
  dat = myreadlines.myreadlines(fname)

 nd  = dat.shape
 nalong = nd[1]
 nvert = nd[0]
 idxmask = np.zeros(nvert)
 nother = len(rejectother)
 
 for idx in range(nalong):
  x = np.median(dat[:,idx])
  
  if (rejectnan == 1):
   idx_xnan = np.isnan(abs(x))
   idxmask[idx_xnan] = 1

  if (rejectinf == 1):
   idx_xinf= abs(x) == np.inf
   idxmask[idx_xinf] = 1

  for i2 in range(nother):
   val = rejectother[i2]
   idx_num = x == val
   idxmask[idx_num] = 1
   
 datnew = dat[idxmask==0,:]
 return(datnew,idxmask)
   
