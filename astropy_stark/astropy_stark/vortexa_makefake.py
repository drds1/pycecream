#test with fake data

import numpy as np
import matplotlib.pylab as plt
import glob
import pandas as pd
import datetime
import statsmodels.api as sm
import mylcgen as mlc
import astropy.convolution as apc


def makefake(tlo,thi,dt,shift = [2,3],noise=0.1,iseed=-1):

 
 dat=mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=tlo,thi=thi,dt=dt,ploton=0,iseed=iseed,meannorm = 0, sdnorm = 1.0)
  
 nlag = len(shift)
 ntimes = np.int((thi - tlo)/dt + 1)
 output = []
 
 lg = np.arange(ntimes)*dt
 today = pd.Timestamp.today()
 for i in range(nlag):
  lagcent = shift[i]
  lagwide = dt
  
  
  datelo = today - pd.Timedelta(days=thi-tlo)
  dates = [datelo + pd.Timedelta(days=i) for i in range(ntimes)]
  times = np.arange(tlo,thi+dt,dt)
  shape_k = np.int(2*(ntimes-1) + 1)
  kernel = np.zeros(shape_k)
  k0 = np.int(np.floor(shape_k/2))
  laggrid = np.zeros(shape_k)
  print(k0,shape_k)
  laggrid[k0:]=lg
  laggrid[:k0]=lg[-1:0:-1]
  
  
  kernel= np.exp(-0.5*((laggrid-lagcent)/lagwide)**2)/(2*np.pi*lagwide)**2
  response = apc.convolve(dat[:,1],kernel)
  
  response_std = np.std(response)
  response_mean = np.mean(response)
  response = (response - response_mean)/response_std
  response = response + np.random.randn(ntimes)*noise
  
  
  if (i == 0):
   dsave_mean = np.mean(dat[:,1])
   dsave_std  = np.std(dat[:,1])
   dsave = (dat[:,1] - dsave_mean)/dsave_std
   dsave = dsave + np.random.randn(ntimes)*noise
   output.append(pd.DataFrame(list(np.array([dates,dsave]).T)))
   

  output.append(pd.DataFrame(list(np.array([dates,response]).T) )) 
  
  
 return(output)
 
 
#test
#a = makefake(0,100,1.0)



