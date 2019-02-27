# python script to convert swift light curves into mJy
# from erg/s/cm2/A
# NOTE if unknown units, anything ~ 10^-14 is most likely in erg/s/cm^2/A (at least for 5548)

import numpy as np


def erg2mjy(fname,convwav):
 dat=np.loadtxt(fname)
 x=dat[:,1]
 sigx=dat[:,2]
 
 
 conv=2.99292458e-5
 
 #see astro.wku.edu/strolger/UNITS.txt  the 1000 goes from jansky to mJy
 a=convwav*convwav/conv * 1000
 
 x=x*a
 sigx=sigx*a
 # save data
 
 dat[:,1]=x
 dat[:,2]=sigx
 
 # save the converted data
 np.savetxt('converted_'+str(fname),dat)
 



def mjy2erg_now(fmjy,wav):

 conv=wav*wav/2.99292458e-5
 
 op = fmjy/conv/1000.
 
 return(op)
 
#input ferg = ferg / 1e-15
def erg2mjy_now(ferg,wav):

 conv=wav*wav/2.99292458e-5
 
 op = ferg*1000*conv/1.e15
 
 return(op)