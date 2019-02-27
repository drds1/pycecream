# python code to take an input light curve add signal to noise and save resave
import myrandom
import numpy as np

#inputs:  Dat 2d array of x,y,sig,
#         SNR desired signal to noise ratio
#         File output filename
# output: the modified data with altered with gaussian noise with the input snr

def myaddnoise(dat,snr):
    x   = dat[:,0]
    y   = dat[:,1]
    sig = dat[:,2]
    
    
    signew=np.abs(y/snr)
    ynew=myrandom.normdis(1,y,signew)
    datnew=np.ones((dat[:,0].shape[0], 3))
    datnew[:,0]=x
    datnew[:,1]=ynew
    datnew[:,2]=signew
    
    
    #for i in range(datnew[:,0].shape()[0]):
    # print datnew[i,0],datnew[i,1],datnew[i,2]
     
     
    return(datnew)    
    
    
    
    
