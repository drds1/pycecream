# python code to calculate the structure function (dy vs dx for pairs of points)
# Inps: x,y array
# Ops : dx, dy the sorted difference in x and y's (should be (nxy - 1)/2 values )
import numpy as np 

def mysf(x,y):

 nxy = x.shape[0]
 
 dx = np.array(())
 dy = np.array(())
 for i in range(nxy-1):
  #print i, x[i].shape[0]
  dxit = x[i+1:] - x[i]
  dyit = y[i+1:] - y[i]
  dx   = np.append(dx,dxit)
  dy   = np.append(dy,dyit)
  
  idx_sort = np.argsort(dx)
  dx   = dx[idx_sort]
  dy   = dy[idx_sort] 
 return(dx,dy)
 
 
 
 
### test it!
#
#import numpy as np
#import pylab as plt
#import mysf
#
#
#fname = 'testfile.dat'
#
#dat = np.loadtxt(dat)
#
#sf = mysf.mysf(x,y)
#dx = sf[0]
#dy = sf[1]
#
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax1.plot(x,y)
#
#ax2 = fig.add_subplot(221)
#ax2.plot(dx,dy)
#
#plt.show()
#
##
