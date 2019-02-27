import os
import numpy as np
# this code takes an input light curve
# truncates it into new light curves of length specified by the user.
# resaves the files

# hardwired options
equaldat = 1 #if equaldat = 1, the code will try to pick roughly even chunks but force them all to have the same number of data points

dir_default = os.getcwd()
dir = raw_input('Enter directory of data: ')
if not dir:
    dir = dir_default


os.chdir(dir)
fname = raw_input('Enter file name: ')


print 'light cruve will be chopped at multiples of the truncation time...'
trunctime = raw_input('Enter truncation time: ')
trunctime = np.double(trunctime)


dat = np.loadtxt(fname)
x   = dat[:,0]
y   = dat[:,1]
sig = dat[:,2]
ndat = x.shape[0]


xlen = x[-1] - x[0]

nchops = np.int( np.floor(xlen/trunctime) )

if (equaldat == 1):
    nequal = ndat/nchops

idxlo = 0
os.mkdir('chopped_'+fname[:fname.find('.')])
os.chdir('chopped_'+fname[:fname.find('.')])
for i in range(nchops):
    os.mkdir(str(i)+'_of_'+str(nchops))
    os.chdir(str(i)+'_of_'+str(nchops))
    
    if (equaldat ==1):
        idx = list(np.arange(idxlo,idxlo+nequal,1))
        idxlo = idxlo + nequal
        if (idxlo + nequal > ndat):
            break
    
    else: 
    	xchoplo  = i*trunctime
    	xchophi  = (i+1)*trunctime
    	idx      = np.where( (x >= xchoplo) and (x < xchophi) )[0]
    
    nidx     = len(idx)
    
    xnew   		= x[idx]
    ynew   		= y[idx]
    signew 		= sig[idx]
    datnew 		= np.zeros((nidx,3))
    datnew[:,0] = xnew
    datnew[:,1] = ynew
    datnew[:,2] = signew
    
    fnamenew = str(i)+'_of_'+str(nchops)+'_'+fname
    
    np.savetxt(fnamenew, datnew)
    os.chdir('../')

os.chdir(dir_default)
    
    
