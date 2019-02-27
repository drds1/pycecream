# function to find the nearest element in an array to some value

import numpy as np

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def find_nearest_idx(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

np.vectorize(find_nearest)    
  
#xin2 and tin2 have fewest points
def find_nearest_fvg(tin1,xin1,tin2,xin2):
    
    xop = find_nearest_vector(xin1,xin2)
    
    # this is the condition for the array having no repeated elements
    if np.array_equal( np.sort(xop), np.unique(xop) ):
        top = tin2
        return(top,xop)
    
    else:
        xdup = xop[ xop[1:] == xop[:-1] ]   #the values of xop which are duplicated
        
        for val in xdup:
            idx = np.where(xop == val)[0]   # identify the indicies of the duplicated values
            idxkeep=[]
            for i in range(t1in):
                idxkeep.append( idx[ np.amin(t1in[i] - t2in[idx]) ] )
                
    # finish the function
    
    return()
                
    
    
    