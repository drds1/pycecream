#!!! ingest input data file, convert dmag to rat flux with mean level, also offset times by
#!!! desired amount

import numpy as np

def myfluxchange(fname, fmean, dtime, dat = 0):
 
 if (type(dat) == int):
  dat = np.loadtxt(fname)
 
 nx = np.shape(dat)
 ny = nx[1]
 nx = nx[0]
 datnew = np.zeros((nx,ny)) 
 print nx, ny
  
 datnew[:,0] = dat[:,0] + dtime		#times
 
 fnew = fmean*10**(-0.4*dat[:,1])	# data points
 
 sigfnew = np.abs (fmean * np.log(10.0) * 0.4*dat[:,2])
 
 
 datnew[:,1] = fnew
 
 datnew[:,2] = sigfnew
 
 for i in range(dat[:,0].shape[0]):
  print i, dat[i,0], dat[i,1], datnew[i,1], dat[i,2]
  raw_input()
 
 np.savetxt('mod_'+fname,dat)
 
 return()

 
 