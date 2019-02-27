#!!!! Extracts shorter light curve between two times and dumps it to another text file. 
#!! See below for arguments

import numpy as np
import pickle
import glob
import os
import sys
from pylab import *

def mypicktimes(dir,lines,tloold,thiold,outfile='creamnames.dat',plotop=0,submin=0):
 #dir     = sys.argv[1] #the directory storing the light curves e.g '../fort/fortcode/fakelc/kep_18mar'
 #lines   = sys.argv[2] # a list of the files .e.g ['file1.dat','file2.dat','...'] etc
 #dtave   = sys.argv[3] #e.g 0.5 will resample with mean half day cadence
 #tloold  = sys.argv[4] #lower time
 #thiold  = sys.argv[5] #upper time
 #outfile = sys.argv[6] #name of file storing the names of the outputted light curves
 
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! arguments above
 #tloold  = 241.0
 #thiold  = 341.0
 ntnew   = thiold - tloold + 1
 toffset = 0.0
 ioffset = [2] #this is the kepler light curve where we apply the offset
 
 pwd    = os.getcwd()
 #dir    = './kep_18mar'#./fake_backups_CREAM_paper'
 #fnames = 'mypicktimesdat.dat'
 #snrfile = 'addnoise_snr.dat'
 #outfile = 'mcmcmultinames.dat'
 #dirs = 'default'
 
 #with open(fnames) as f:
 #    lines = f.read().splitlines()
 
 
 os.chdir(dir)
 
 #dirprev = glob.glob('*'+dirs+'*')
 #nprev = len(dirprev)
 #dirnew = dirs+'_'+str(nprev+1)
 #os.mkdir(dirnew)
 #pwd = os.getcwd()
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! read files lcfile and snr file
 #  
 
 
     
 nf = len(lines)
 
 fnamenew=[]
 for i in range(nf):
  dat = np.loadtxt(lines[i])
  
  
  if (i in ioffset):
   tlo = tloold
   thi = thiold
  else:
   tlo = tloold + toffset
   thi = thiold + toffset
  
  t = dat[:,0]
  print np.min(t),np.max(t),tlo,thi
  idxlo = np.where(t > tlo)[0][0]
  
  nhi   = np.where(t > thi)[0].shape[0]
  if (nhi > 0):
   idxhi = np.where(t > thi)[0][0]-1 
  else:
   idxhi = t.shape[0] - 1
  
  ntnew = idxhi - idxlo
  datnew = np.zeros((ntnew,3)) 
  datnew[:,0] = dat[idxlo:idxhi,0]
  datnew[:,1] = dat[idxlo:idxhi,1]
  datnew[:,2] = dat[idxlo:idxhi,2] 
  
  
  if (submin == 1):
   datnew_min = np.min(datnew[:,0])
   datnew[:,0] = datnew[:,0] - datnew_min
   
  if (plotop == 1):
   errorbar(datnew[:,0],datnew[:,1],datnew[:,2],ls='')
   show()
   
  
  fnamenew.append( 'length_'+lines[i] )
  
  #os.chdir(dirnew)
  np.savetxt(fnamenew[i],datnew)
  #os.chdir(pwd)
  #os.chdir(dir)
  
  
  
 #write to file
 #os.chdir(dirnew)
 of = open(outfile,'w') 
 for item in fnamenew:
   of.write("%s\n" % item)
  
 of.close()
 os.chdir(pwd)
 return()