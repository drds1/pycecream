import numpy as np
from mylcgen import *
from myrandom import *
import matplotlib.pylab as plt
import mysinecostrans as msc

#generate test light curve and add noise
dat      = mylcgen(tlo=0,thi=100,dt=2.0,iseed=42342332)
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

flo = 0.5/(thi-tlo)
fhi = 5./dt

msc.sct(time,y,sig,freq=np.arange(flo,fhi+flo,flo),costrans = 0)


#test the routine above and do it manually to compare

tlo = np.min(time)
thi = np.max(time)
dt  = np.mean(time[1:]-time[:-1])




#if cosine tranform, use only coses, else if sine transform use only sines
costrans = 0




if (costrans == 1):
 phase = 0
elif (costrans == 0):
 phase = -np.pi/2


#input some minimum and maximum frequency range for the Fourier terms
flo = 0.5/(thi-tlo)
fhi = 5./dt

#angular frequency
w = np.arange(flo,fhi+flo,flo)*2*np.pi

#now compute the hessian matrix
nw = np.shape(w)[0]
np2 = 2*nw

cvec = np.ones(nw)
hes = np.ones((nw,nw))
sig2 = sig*sig
#save the sine and cosine value in  2d array
cwt = np.ones((ndat,nw))
swt = np.ones((ndat,nw))
for iw in range(nw):
 cwt[:,iw] = np.cos(w[iw]*time +phase)

for ikp in range(nw):
 wp = w[ikp]
 
 wp2 = wp*wp
 cvec[ikp]   = np.sum(y*cwt[:,ikp]/sig2)
 #print ikp,nw
 for ik in range(nw):
  #cosine column
  hnow_c = np.sum(cwt[:,ik]*cwt[:,ikp]/sig2)
  if (ik == ikp):
   hnow_c = hnow_c + wp2
  hes[ikp,ik]   = hnow_c

   

      
#invert the matrix to get the covariance matrix
cov = np.linalg.inv(hes)


#find the parameters
parm  = cov.dot(cvec)
ckout = parm









#plot the input and output light curves to see if they match
freq = w/2/np.pi
dtres = dt/10
tplot = np.arange(tlo,thi+dtres,dtres)
xplot = []
for tnow in tplot:
 xplot.append( np.sum(ckout*np.cos(w*tnow + phase) ) )

fig = plt.figure()
ax1 = fig.add_subplot(211) 
ax1.errorbar(time,y,sig,ls='',color='k')
ax1.plot(tplot,xplot,color='b')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('flux')


ax2 = fig.add_subplot(212) 
ax2.plot(freq,ckout**2,ls='',marker='o')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Frequency (cyc/day)')
ax2.set_ylabel(r'power ($S_k^2 + C_k^2$)')
plt.savefig('test_r_hesfit_cosine.pdf')


#
#
#
#
#
#
#


    