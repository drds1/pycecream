#python  code to give me the range of a set of data
import numpy as np
import os

pwd=os.getcwd()

dir='../fort/fortcode/mcmcmulti3/5548ltonly_july28/'
os.chdir(dir)

		
	
			
	
fname=['converted_u_mjd.dat','converted_g_mjd.dat','converted_r_mjd.dat','converted_i_mjd.dat','converted_z_mjd.dat']#['resamp_fake_3556.51_nonoise.dat','resamp_fake_4702.48_nonoise.dat','resamp_fake_6175.57_nonoise.dat','rebinned_fake_6300.00_nonoise.dat','resamp_fake_7489.97_nonoise.dat','resamp_fake_8946.70_nonoise.dat']
nfile=len(fname)
fracvar=np.zeros(nfile)
mean=np.zeros(nfile)
sd=np.zeros(nfile)


for i in range (nfile):
    dat=np.loadtxt(fname[i])
    sd[i]=np.std(dat[:,1])
    mean[i]=np.mean(dat[:,1])
    fracvar[i]=sd[i]**2/mean[i]**2
    
print 'the mean s...'
print mean
print 'the sd s'
print sd
print 'the fracvar s...'
print fracvar
print 'the end'

os.chdir(pwd)