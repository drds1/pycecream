import os
import numpy as np
from pylab import *
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

ndat = x.shape[0]

xdiff = np.abs(y[1:] - y[:-1])
xdiffmean = np.mean(xdiff)
xdiffsort = np.sort(xdiff)

nsections = 12
nchops = nsections - 1
cutoff = xdiffsort[-nchops:]
cutoff = np.insert(cutoff,0,0)

#idxchop = np.where( xdiff > cutoff)[0]#xdiffmean*10)[0]
#idxchop = np.insert(idxchop,0,0)

#nsections = idxchop.shape[0]
dirsave = 'chopped_mean'
direxist = os.path.isdir(dirsave)

if (direxist == True):
    os.system('rm -rf '+dirsave)

os.mkdir('chopped_mean')
os.chdir('chopped_mean')
for i in range(nsections-1):

    if (i == 0):
        idxlo = 0
    else:
        idxlo = np.where(xdiff == cutoff[i])[0]
    
    idxhi = np.where(xdiff == cutoff[i+1])[0]

    xnew = x[idxlo:idxhi]
    ynew = y[idxlo:idxhi]
    signew = sig[idxlo:idxhi]
    nnew = xnew.shape[0]
    
    datnew = np.zeros((nnew,3))
    datnew[:nnew,0] = xnew
    datnew[:nnew,1] = ynew
    datnew[:nnew,2] = signew
    np.savetxt('chopped_'+str(i)+'.dat',datnew)
    
    errorbar(xnew,ynew,signew,linestyle='None')
    savefig('chopped_'+str(i)+'.png')
    clf()
    
os.chdir('../')
