# redistribute noise on gaussian light curve
import numpy as np
import os

pwd=os.getcwd()

dir='../fort/fortcode/fakelc/'
os.chdir(dir)

		
	
fname=['resamp_fake_3556.51_nonoise.dat','resamp_fake_4702.48_nonoise.dat','resamp_fake_6175.57_nonoise.dat','rebinned_fake_6300.00_nonoise.dat','resamp_fake_7489.97_nonoise.dat','resamp_fake_8946.70_nonoise.dat']
nfile=len(fname)

for i in range(nfile):
    dat=np.loadtxt(fname[i])
    n=dat[:,0].shape[0]
    dat[:,1]=dat[:,1]+np.random.randn(n)*dat[:,2]

    np.savetxt('refuxed_'+fname[i],dat)


os.chdir(pwd)