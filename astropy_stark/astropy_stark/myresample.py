#### code to randomly resample a set of input light curves
#10/9/2017 sampcode 3 update includes sampmin the minum spacing between data points
#sampcode 4, dtave indicates the minimum space between data points, data points will be selected (with no interpolation) from the parent sample, skipping points until the minimum spacing dtave is achieved
##avesamp is the average length of time between the random samples
#set dir = '' and fname=[''] to have this function work on datin[nt,3] and output array rather than save to file

#new 10/9/2017 added option dtsdin need mean and standard deviation of spacing between points e.g setting
#dtsdin very small will give very regularly spaced points.
#if negative then the absolute value is the fraction relative to the mean e.g -0.2 will set the 
#standard deviation as a fifth of the mean spacing between points


import numpy as np
import os
import sys
import astropy_stark.myrandom as mr


def myresample(dir,fname,dtave,dtsdin = -0.2, sampmin=0.8,sampcode=3,datin=[]):
 
 if (dtsdin < 0):
  dtsd = np.abs(dtsdin)*dtave
 else:
  dtsd = dtsdin
 
 dirpy=os.getcwd()
 #dir   = sys.argv[1] #the directory storing the light curves e.g '../fort/fortcode/fakelc/kep_18mar'
 #fname = sys.argv[2] # a list of the files .e.g ['file1.dat','file2.dat','...'] etc
 #dtave = sys.argv[3] #e.g 0.5 will resample with mean half day cadence
 #sampmin   = 0.8
 #sampcode = 2
 
 #!!### user arguments above. Don't change sampcode or sampmin unless you know what they do (I don't and I wrote the code).###
 if (dir != ''):
  os.chdir(dir)
 
 Nfile=len(fname)
 
 for ifile in range(Nfile):
     if (fname[ifile] == ''):
      dat = datin
     else:
      dat=np.loadtxt(fname[ifile])
     t=dat[:,0]
     x=dat[:,1]
     sig=dat[:,2]
     Ndat=t.shape[0]
     dt = (t[-1] - t[0])/(Ndat-1)
     
 # below are two versions of the code (the 2nd should be more sophisticated and consider the approximate spacing between each point when making its idxsamp selection    
     if sampcode == 1:
     	nidx=(1.-sampmin)*np.random.ranom_sample(1)[0]+sampmin
     	idxsamp=np.random.rand(low=0,high=Ndat,size=nidx)
     	datsamp=np.zeros((nidx,3))
     	datsamp[:,0]=t[idxsamp]
     	datsamp[:,1]=x[idxsamp]
     	datsamp[:,2]=sig[idxsamp]
     	
     elif sampcode == 2:
         idxcount=0
         tthen=t[0]
         idxsamp=[]
         xn   = []
         sign = []
         tn   = []
         while (idxcount < Ndat) & (tthen < t[-1]):
             
 
             a = np.random.randn(1)*dt*2
             
             tnow = tthen + dtave + a
             tn.append(tnow)
             xn.append(np.interp([tnow],t,x)[0])
             sign.append(np.interp([tnow],t,sig)[0])
             #idxsamp.append(np.abs(t-tnow).argmin()) ## index of closest time to tnow
             
             tthen=tnow
             idxcount=idxcount+1
             
         #idxsamp=np.array(idxsamp)
         tn = np.array(tn)
         xn = np.array(xn)
         sign = np.array(sign)
         nn = xn.shape[0]
         datsamp=np.zeros((nn,3))
         datsamp[:,0]=tn[:,0]
         datsamp[:,1]=xn[:,0]
         datsamp[:,2]=sign[:,0]
 
 
 
     elif sampcode == 3:
         #print 'ararar'
         idxcount=0
         tthen=t[0]
         idxsamp=[]
         tlast = t[-1]
         while (idxcount < Ndat-1) & (tthen < tlast - 4*sampmin):
             
             
             #a = np.random.randn(1)*dtave
             a = np.random.normal(dtave,dtsd,1)[0]
             
             tnow = tthen +  np.abs(a)
             idxtemp = np.abs(t-tnow).argmin()
             
             #print tnow,tthen,'before mrs'
             #print tnow, tthen, idxcount, sampmin
             if ((idxtemp not in idxsamp) and ((tnow - tthen > sampmin) or (tnow > tlast - sampmin))):
              idxsamp.append(idxtemp) ## index of closest time to tnow
              idxcount=idxcount+1
              a= 1.*tthen-tnow
              tthen=tnow
             
             #print idxcount, Ndat, tthen, t[-1],'mrs'
         idxsamp=np.array(idxsamp)
         datsamp=np.zeros((idxsamp.shape[0],3))
         datsamp[:,0]=t[idxsamp]
         datsamp[:,1]=x[idxsamp]
         datsamp[:,2]=sig[idxsamp]
         
         ttemp_space = datsamp[1:,0]-datsamp[:-1,0]
         #print('min,max,and ave spacing between elements',ttemp_space.min(), ttemp_space.max(), np.mean(ttemp_space))
 
     
 
 
     elif sampcode == 4:
         idxcount=0
         tthen=t[0]
         idxsamp=[]
         while (idxcount < Ndat) & (tthen < t[-1]):
             
             
             tnow = tthen +  dtave
             
             b = t>tnow
             idxtemp = [i for i, elem in enumerate(b, 1) if elem]
             if (len(idxtemp) ==0):
              break
             
             idxtemp = idxtemp[0]
             if (idxtemp >= t.shape[0]):
              break
              
             if (idxtemp not in idxsamp):
              idxsamp.append(idxtemp) ## index of closest time to tnow
              idxcount=idxcount+1
             
             a= tnow - tthen
             tthen=t[idxtemp]
             
             #if (t[idxtemp] - t[idxtemp-1] < dtave):
              #print 'PROBLEM!!',
              #print t[idxtemp] - t[idxtemp-1], dtave#print tnow, idxcount, dtave, dt, tthen, a
              #raw_input()
         #print idxsamp
         idxsamp=np.array(idxsamp)
         datsamp=np.zeros((idxsamp.shape[0],3))
         datsamp[:,0]=t[idxsamp]
         datsamp[:,1]=x[idxsamp]
         datsamp[:,2]=sig[idxsamp]
         
         ttemp_space = datsamp[1:,0]-datsamp[:-1,0]
         ns = len(idxsamp)
         #for i in range(ns):
         # print i,datsamp[i,:]
         #print('min,max,and ave spacing between elements',ttemp_space.min(), ttemp_space.max(), np.mean(ttemp_space))
         #print('locations...',ttemp_space.argmin(),ttemp_space.argmax())
 
 
     
     np.savetxt('resamp_'+fname[ifile],datsamp)   
 
 os.chdir(dirpy) # change back to python directory
 return(datsamp)
 
 
 
 
 

