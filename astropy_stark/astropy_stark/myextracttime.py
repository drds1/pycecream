import numpy as np
import os



dir = '/Users/ds207/Documents/standrews/sta/fort/fortcode/mcmcmulti3/aas_rebinned/'
file_list = 'mcmcmultinames.dat'
tlo = 6690.0
thi = 6810.0



pwd = os.getcwd()
os.chdir(dir)
f = open(file_list)
files = f.read().splitlines()
f.close()
nf = len(files)


for fnow in files:
 print 'reading ',fnow
 dat = np.loadtxt(fnow)
 idxtrue = ((dat[:,0] >= tlo) & (dat[:,0] <= thi))
 datnew = dat[idxtrue,:]
 
 np.savetxt('extracted_'+fnow,datnew)
 

os.chdir(pwd) 
 