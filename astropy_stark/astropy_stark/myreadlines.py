## read a file and dont give me crappy errors
## ouptut numpy 2d array (nrows,ncols)
#alterd code to deal with lines with uneven number of lines returns rejected points in a fail lsit

import numpy as np



def myreadlines(fname):

 f = open(fname)
 
 faillist = []
 idx = 1
 a = []
 for i in f:
  #a.append(np.float(i.split()))
  b = []
  new = i.split()
  for idx2 in range(len(new)):
   #produce exception if there is a hidden character in there
   try:
    b.append(np.float(new[idx2]))
   except:
    print 'myreadlines.py: Cannot convert to float index',idx,idx2
    #print new
    #print len(new),new[0]
    b.append(np.nan)
  a.append(b)
  idx = idx+1
  #print f.readline()
 #a = np.array(a)
 #
 
 
 nrow = len(a)
 
 #how to deal with files that have an uneven number of columns (find the mode number of columnsa and put that info in array
 #skipping elements that dont fit
 ncollist = [len(a[i]) for i in range(nrow)]
 ncol = max(set(ncollist), key=ncollist.count)
 
 #ncol = len(a[0])
 
 b = np.zeros((nrow,ncol))
 idcount = 0
 for i in range(nrow):
  try:
   b[idcount,:] = a[i][:]
   idcount = idcount + 1
  except:
   print 'something wrong with row',i,'...skipping'
   faillist.append(a[i][:])
 b = b[:idcount,:]  
   
 f.close()
 return(b,faillist)
 

