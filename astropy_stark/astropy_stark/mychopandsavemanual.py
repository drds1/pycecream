import os
import numpy as np
# this code takes an input light curve
# truncates it into new light curves of length specified by the user.
# resaves the files

# hardwired options
equaldat = 1 #if equaldat = 1, the code will try to pick roughly even chunks but force them all to have the same number of data points

# filter twist th long kepler light curve periodically rotates, if filtertwist = 1, the light curve will be divided into segments which do not have a break in the filter 
filtertwist = 1

dir_default = os.getcwd()
dir = raw_input('Enter directory of data: ')
if not dir:
    dir = dir_default


os.chdir(dir)
fname = raw_input('Enter file name: ')




dat = np.loadtxt(fname)
x   = dat[:,0]
y   = dat[:,1]
sig = dat[:,2]




# manually select the indicies of the x values to chop
idxchophi = [123,1234,2256,34556]
idxchophi = np.array(idxchoplo)
idxchophi  = np.sort(idxchoplo)

nidx = idxchoplo.shape[0]

os.mkdir('mychopandsave_manual')
os.chdir('mychopandsave_manual')

for i in range(nidx):
    idxlo = idxchophi[i]
    idxhi = idxchophi[i+1]
    
    xnew = x[idxlo:idxhi]
    ynew = y[idxlo:idxhi]
    signew = sig[idxlo:idxhi]
    
    nnew = xnew.shape[0]
    
    datnew = np.zeros((nnew,3))
    
    datnew[:,0] = xnew
    datnew[:,1] = ynew
    datnew[:,2] = signew
    
    np.savetxt('chop_'+str(i)+fname, datnew)
    
os.chdir('../')