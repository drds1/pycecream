#combine 2textfiles into 1
import myreadlines 
import numpy as np

def mycomb2to1(fname,fnameop='',cutoff=[]):
 
 #
# with file(fname) as f:
#  fname = f.read().splitlines()
 
 nf  = len(fname)
 cutoff_in = []*nf
 #print fname,len(fname)
 
 if (len(cutoff) == 0):
  cutoff_in[:] = 0
 else:
  cutoff_in[:] = cutoff[:]
 
 
 dat=[]
 ncol = []
 for i in range(nf):
  try:
   temp = np.loadtxt(fname[i])
  except:
   temp = myreadlines.myreadlines(fname[i])
  
  print 'mycomb2to1.py:', fname[i],nf, temp.shape[0]
  #only read after the cutoff
  if (cutoff_in[i] < 0):
   nrolold = temp.shape[0] 
   nro =  abs(nrolold*cutoff_in[i])
  else:
   nco = cutoff_in[i]
  
  temp = temp[nro:,:]
   
  #make sure both files have same number of columns
  ncol.append(temp.shape[1])
  ncolmin = min(ncol)
  if (all(x==ncol[0] for x in ncol) ==0):
   print 'mycomb2to1 warning: different number of colums for the two arrays you are trying to combine. I will combine to the lowest number...'  
   print 'column lengths...', ncol
  dat.append(temp)
 
 datnew=[]
 for i in range(nf):
  datnew.append(dat[i][:,:ncolmin])
   
 datop = np.concatenate((datnew[:]),axis=0)
 
 print ''
 print 'mycomb2to1.py: after combining. New umber of points=',np.shape(datop) 
 if (fnameop != ''):
  np.savetxt(fnameop,datop)
 
 return(datop)