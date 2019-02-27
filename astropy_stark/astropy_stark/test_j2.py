#code to test a more javelin styke way of running cream
import numpy as np
import mylcgen as mg
from myrandom import *
import matplotlib.pylab as plt
import myredcisq as mrc


#generate test light curve and add noise
dat      = mg.mylcgen(tlo=0,thi=100,dt=2.0)
ndat     = np.shape(dat[:,0])[0]
datmean  = np.std(dat[:,1])
sig      = np.ones(ndat)/10*datmean
for i in range(ndat):
 dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
y        = dat[:,1]
time     = dat[:,0]

tlo = np.min(time)
thi = np.max(time)
dt  = np.mean(time[1:]-time[:-1])

dtres = dt/2



#fakelc making done now do actual routine
tgrid = np.arange(tlo,thi+dtres,dtres)
ngrid = np.shape(tgrid)[0]
nits  = 1000



#make a fake light curve for each teration and compute a fit statistic for each
ymsum  = np.zeros(ngrid)
ym2sum = np.zeros(ngrid)
wsum   = 0
cisq0  = 0
cisqsave = []
for i in range(nits):
 #print i
 dat      = mg.mylcgen(tlo=tlo,thi=thi,dt=dtres)
 cisq     = mrc.myredcisq(time,y,sig,dat[:,0],dat[:,1],npar = 0)[0]
 if (i == 0):
  cisq0 = cisq
 
 cisqsave.append(cisq)
 print i, cisq, cisq0
 w        = cisq#np.exp(-(cisq-cisq0)/2)
 wsum     = wsum   + w
 ymsum    = ymsum  + w*dat[:,1]
 ym2sum   = ym2sum + w*dat[:,1]**2

yplot   = ymsum / wsum / nits
ysdplot = (ym2sum - ymsum) / wsum / nits










#plot the results
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(tgrid,yplot,label=None,color='b')
ax1.errorbar(time,y,sig,color='k',label=None,ls='')
plt.savefig('test_j2.pdf')




