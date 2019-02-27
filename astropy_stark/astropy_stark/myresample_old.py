#### code to randomly resample a set of input light curves

##avesamp is the average length of time between the random samples

import numpy as np
import os
import sys

def myresample(dir,fname,dtave,sampmin=0.8,sampcode=2):
 dirpy=os.getcwd()
 #dir   = sys.argv[1] #the directory storing the light curves e.g '../fort/fortcode/fakelc/kep_18mar'
 #fname = sys.argv[2] # a list of the files .e.g ['file1.dat','file2.dat','...'] etc
 #dtave = sys.argv[3] #e.g 0.5 will resample with mean half day cadence
 #sampmin   = 0.8
 #sampcode = 2
 
 #!!### user arguments above. Don't change sampcode or sampmin unless you know what they do (I don't and I wrote the code).###
 
 os.chdir(dir)
 
 Nfile=len(fname)
 
 for ifile in range(Nfile):
     dat=np.loadtxt(fname[ifile])
     t=dat[:,0]
     x=dat[:,1]
     sig=dat[:,2]
     Ndat=t.shape[0]
     dt = (t[-1] - t[0])/(Ndat-1)
     
 # below are two versions of the code (the 2nd should be more sophisticated and consider the approximate spacing between each point when making its idxsamp selection    
     if sampcode == 1:
     	nidx=(1.-sampmin)*np.random.ranom_sample(1)[0]+sampminn
     	idxsamp=np.random.rand(low=0,high=Ndat,size=nidx)
     	datsamp=np.zeros((nidx,3))
     	datsamp[:,0]=t[idxsamp]
     	datsamp[:,1]=x[idxsamp]
     	datsamp[:,2]=sig[idxsamp]
     	
     elif sampcode == 2:
         idxcount=0
         tthen=t[0]
         idxsamp=[]
         while (idxcount < Ndat) & (tthen < t[-1]):
             
 
             a = np.random.randn(1)*dt
             
             tnow = tthen + dtave + a
             idxsamp.append(np.abs(t-tnow).argmin()) ## index of closest time to tnow
             
             tthen=tnow
             idxcount=idxcount+1
             
         idxsamp=np.array(idxsamp)
         datsamp=np.zeros((idxsamp.shape[0],3))
         datsamp[:,0]=t[idxsamp]
         datsamp[:,1]=x[idxsamp]
         datsamp[:,2]=sig[idxsamp]
 
 
 
     np.savetxt('resamp_'+fname[ifile],datsamp)   
 
 os.chdir(dirpy) # change back to python directory
 return()