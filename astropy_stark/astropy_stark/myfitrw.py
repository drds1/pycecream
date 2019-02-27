import numpy as np
from mylcgen import *
from myrandom import *
import matplotlib.pylab as plt
import myresample as mrs

#generate test light curve and add noise
datpre      = mylcgen(tlo=0,thi=100,dt=0.03,iseed=34245)

 

ndat     = np.shape(datpre[:,0])[0]
datmean  = np.std(datpre[:,1])
sig      = np.ones(ndat)/10*datmean

dat = mrs.myresample(dir='',fname=[''],dtave=2.0,sampmin=0.8,sampcode=3,datin=np.array((datpre[:,0],datpre[:,1],sig)).T)
ndat = np.shape(dat[:,0])[0]
sig = dat[:,2]


for i in range(ndat):
 dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
y        = dat[:,1] +50
time     = dat[:,0]

tlo = np.min(time)
thi = np.max(time)
dt  = np.mean(time[1:]-time[:-1])

#input some minimum and maximum frequency range for the Fourier terms
flo = 0.5/(thi-tlo)
fhi = 2./dt

#angular frequency
w = np.arange(flo,fhi+flo,flo)*2*np.pi
w = np.arange(0,fhi+flo,flo)*2*np.pi
w[0] = flo/50


#now compute the hessian matrix
nw = np.shape(w)[0]
np2 = 2*nw

cvec = np.ones(np2)
hes = np.ones((np2,np2))
sig2 = sig*sig
#save the sine and cosine value in  2d array
cwt = np.ones((ndat,nw))
swt = np.ones((ndat,nw))
for iw in range(nw):
 cwt[:,iw] = np.cos(w[iw]*time)
 swt[:,iw] = np.sin(w[iw]*time)

iwp = 0
for ikp in range(0,np2,2):
 wp = w[iwp]
 
 wp2 = wp
 cvec[ikp]   = np.sum(y*cwt[:,iwp]/sig2)
 cvec[ikp+1] = np.sum(y*swt[:,iwp]/sig2)
 print iwp,nw
 iw = 0
 for ik in range(0,np2,2):
  #cosine column
  hnow_c = np.sum(cwt[:,iw]*cwt[:,iwp]/sig2)
  hnow_s = np.sum(swt[:,iw]*cwt[:,iwp]/sig2)
  if (ik == ikp):
   hnow_c = hnow_c + wp2
  hes[ikp,ik]   = hnow_c
  hes[ikp,ik+1] = hnow_s 
   
  #sine column
  hnow_c = np.sum(cwt[:,iw]*swt[:,iwp]/sig2)
  hnow_s = np.sum(swt[:,iw]*swt[:,iwp]/sig2)
  if (ik == ikp):
   hnow_s = hnow_s + wp2
  hes[ikp+1,ik]   = hnow_c
  hes[ikp+1,ik+1] = hnow_s 
  iw = iw + 1
  
 iwp = iwp + 1
      
#invert the matrix to get the covariance matrix
cov = np.linalg.inv(hes)


#find the parameters
parm  = cov.dot(cvec)
skout = parm[1::2]
ckout = parm[::2]










#plot the input and output light curves to see if they match
freq = w/2/np.pi
dtres = dt/10
tplot = np.arange(tlo,thi+dtres,dtres)
xplot = []
for tnow in tplot:
 xplot.append( np.sum(ckout*np.cos(w*tnow) + skout*np.sin(w*tnow)) )

fig = plt.figure()
ax1 = fig.add_subplot(211) 
ax1.errorbar(time,y,sig,ls='',color='k')
ax1.plot(tplot,xplot,color='b')
ax1.plot(tplot,ckout[0]*np.cos(w[0]*tplot) + skout[0]*np.sin(w[0]*tplot))
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('flux')


ax2 = fig.add_subplot(212) 
ax2.plot(freq,ckout**2 + skout**2,ls='',marker='o')
ax2.plot(freq,skout**2,ls='',color='g',marker='o')
ax2.plot(freq,ckout**2,ls='',color='b',marker='o')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Frequency (cyc/day)')
ax2.set_ylabel(r'power ($S_k^2 + C_k^2$)')
ax2.grid(True)
plt.savefig('test_r_hesfit.pdf')











    