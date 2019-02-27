# convert time,data,error from mjd,mjd-5000,tohjd
import numpy as np
import os


########################definitions
def hjd2mjd(x):
    x=x-2400000.5
    return(x)
     
def mjd2hjd(x):
    x=x+2400000.5
    return(x)
 
########################

#inps the directory of light curves dir//diragn, fname input file name, sring_func what to do to the times

def mydateconv(dir,diragn,fname,string_func):
 #dir='./5548data/gisella_18nov/'#'../fort/fortcode/mcmcmulti3/'
 #diragn=''#'5548ltonly_july28/'
 pwd=os.getcwd()
 
 #a=raw_input('Enter dir path here if you want to change from :- '+str(dir))
 #if (a != ''):
 #    dir=a
 
 #a=raw_input('Enter agn path here if you want to change from :- '+str(diragn))
 #if (a != ''):
 #    diragn=a
 
 #fname=raw_input('Enter file name here :-')
 fnameop='converted_'+fname
 #a=raw_input('Enter new file name here if you dont like :-'+fnameop)
 #if (a != ''):
 #    fnameop=a
 
 os.chdir(str(dir)+str(diragn))
 
 #print 'write 1 of the following:'
 #print '    hm_num : convert HJD to MJD +num'
 #print '    mh_num : convert MJD to HJD +num'
 #print ':-'
 #a=raw_input()
 a = string_func
 numadd=int(a[a.index('_')+1:])
 a=a[:2]
 

 
 
 dat=np.loadtxt(fname)
 
 ## check the data is in ascending order by date
 idxsort=np.argsort(dat[:,0])
 dat[:,:]=dat[idxsort,:]
 ##
 
 ## Act on input instructions
 if (a=='hm'):
     dat[:,0]=hjd2mjd(dat[:,0])+numadd
 elif (a=='mh'):
     dat[:,0]=mjd2hjd(dat[:,0])+numadd
 else:
     print 'I know not of this jibber jabber you speak, try again.'
 
 
 np.savetxt(fnameop,dat)
 
 os.chdir(pwd)