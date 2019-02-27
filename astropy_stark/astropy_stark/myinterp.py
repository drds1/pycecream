import scipy.interpolate as si
import numpy as np
import scipy

#inputs
#x,y,z data
#xi,yi desired interpolation grid

def my2ditp(x,y,z,xi,yi,gauss=0):
 
 #do standard stuff
 xim, yim = np.meshgrid(xi, yi)
 rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
 zim = rbf(xim, yim)
 
 nia = np.shape(zim)[0]

 
 
 ri = xi**2 + yi**2
 r  = x**2 + y**2
 #now alter zi so that it is treated differently when data is extrapolated
 for ia in range(nia):
  
  xinow = xim[:,ia]
  yinow = yim[:,ia]
  rinow = np.sqrt(yinow**2 + xinow**2)
  
  xinow0 = xinow[0]
  if (xinow0 > 0):
   try:
    idx = np.where(x > xinow0)[0]
    ymax = np.max(y[idx])
    idxim = np.where(yinow > ymax)[0]
    zim[idxim,ia] = 0
   except:
    pass  
  elif (xinow0 < 0):
   try:
    idx = np.where(x < xinow0)[0]
    #print idx
    ymax = np.max(y[idx])
    idxim = np.where(yinow > ymax)[0]
    zim[idxim,ia] = 0
   except:
    pass
  
  #print xinow0
  #print yinow
  #print zim[:,ia]
  #raw_input()
 return(zim)
   