#sort telescope ids that may not be in order into order
# sort list of integers that include gaps and missing number where 1 number may be huge

import numpy as np



def mytelid_sort(telid):

 a  = list(telid)#[1,2,17,8,1,3,17,5,5,6,8,0,6]
 idasort = np.argsort(a)
 ida_sort_sort = np.argsort(idasort)
 asort = np.array(a)[idasort]
 
 b  = np.arange(np.min(asort),np.max(asort),1)
 c  = np.array(list( set(list(asort)) ^ set( list(b) ) ))
 anew = []
 
 nasort = len(asort)
 
 for i in range(nasort):
  
  anow = asort[i]
  idxless = list(np.where(c < anow)[0])
  idr = 0
 
  if (i < nasort-1):
   anext = asort[i+1]
  if ((anow != 0) & (len(idxless) > 0)):
   
   idxless = list(np.where(c < anow)[0])
   
   idpop = idxless[0]
   
   if (anow != anext):
    d = c[idpop]
    c = np.delete(c,idpop)
    c = np.insert(c,[0],anow)
    c = np.sort(c)
   else:
    d = c[idpop]
   
   anow = d
  

  anew.append(anow)
 
 anew = np.array(anew)[ida_sort_sort]
 
 return(anew)