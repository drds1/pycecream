#NOTE when rebinning if you have a SNR in mind, remember what cadence the SNR applies to. i.e ~500 for 0.5 day sampled kepler light curves. If you have finely sampled data and you bin it up. SNR will improve. 
#program to take an input time series and rebin it according to the parameters in the input
# rebin optimally averages all data within the new cadence range and caluclates errors according to keith ada notes lecture 4
import numpy as np
import os
import sys


def myrebin(dir,files,cadence):
 #dir=sys.argv[1]       #the directory storing the light curves e.g '../fort/fortcode/fakelc/kep_18mar'
 #files = sys.argv[2]   # a list of the files .e.g ['file1.dat','file2.dat','...'] etc
 #cadence = sys.argv[3] #e.g 0.5 will rebin into half day cadence
 
 
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 3 user inputs defined above!!!!!!!!!!!!!!!!!
 
 
 dirpy = os.getcwd()
 os.chdir(dir)
 
 dirreal = os.getcwd()
 
 nfiles=len(files)
 
 
 #fname=raw_input('Please enter file name:-')
 for ifile in range(nfiles):    
     fname=files[ifile]
     dat=np.loadtxt(fname) #load the data
     
     t=dat[:,0]
     x=dat[:,1]
     sig=dat[:,2]
     ndat=t.shape[0]
     trange=t[-1]-t[0]
     dtbefore=trange/(ndat-1)
     nbin=int(np.ceil(trange/cadence))
     
     op=np.zeros((nbin,3))
     
      
     for i in range(nbin):
         binlo=t[0]+i*cadence
         binhi=binlo+cadence
         binmid=binlo+cadence/2
         idx=np.where((t>=binlo) & (t<binhi))[0] # array of indices of all elements to go in this bin
         bintime = np.median(t[idx])
         
         op[i,0]=bintime
         top=(x[idx]/sig[idx]**2).sum()
         bot=(1./sig[idx]**2).sum()
         op[i,1]=top/bot
         sigbin=1/np.sqrt(bot)
         op[i,2]=sigbin
         
     print 'consider your data re-binned!' 
     
 #    os.chdir(path)
     op=op[np.where(op[:,1] == op[:,1])[0],:]
     
     #os.chdir(dirnew)
     np.savetxt('rebinned_'+str(fname),op)
     #os.chdir(dirreal)
    
 os.chdir(dirpy)
 return()