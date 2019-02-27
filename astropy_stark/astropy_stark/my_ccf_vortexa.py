import numpy as np


import numpy as np
import my_kateccf as mccf
import mylcgen as mlc
import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import scipy.signal as ss

ilag = 5
sigma = 0.3

#generate synthetic data
drive = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=100,dt=1.0,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
 
t = drive[:,0]
x = drive[:,1]
xsd = np.std(x)
nx = np.shape(x)[0] 
x = (x - np.mean(x))/xsd
xsd = np.std(x)
sigx = np.ones(nx)*xsd*sigma
x = x + np.random.randn(nx)*sigx

lc1 = np.array([t,x,sigx]).T

#define convolution kernel
conv = np.zeros(nx)
conv[ilag-1:ilag+1] = 1



echo  = np.convolve(x, conv, mode='same')
echo2 = np.convolve(x, conv, mode='full')[:nx]
 
echo2 = (echo2 - np.mean(echo2))/np.std(echo2) 
esd = np.std(echo2)
ne2=np.shape(echo2)[0]
sige = np.zeros(ne2)+esd*sigma
echo2 = echo2 + np.random.randn(ne2)*esd*sigma
lc2 = np.array([t,echo2,sige]).T




frac = 0.8

dt = 0.1
tmin = min(np.min(lc1[:,0]),np.min(lc2[:,0]))
tmax = max(np.max(lc1[:,0]),np.max(lc2[:,0]))

n1 = np.shape(lc1[:,0])[0]
n2 = np.shape(lc2[:,0])[0]
#interpolate time series onto regular sampled grid
tgrid = np.arange(tmin.tmax,dt)
ngrid = np.shape(tgrid)[0]
id1 = np.arange(n1)
id2 = np.arange(n2)


nsub_1 = np.int(frac*n1)
nsub_2 = np.int(frac*n2)



for iteration in range(nits):
 
 idsub1 = np.random.choice(id1,size = nsub_1)#replace true by default
 idsub2 = np.random.choice(id2,size = nsub_2)
 
 y1_itp = np.interp(tgrid,lc1[idsub1,0],lc1[idsub1,1])
 y2_itp = np.interp(tgrid,lc2[idsub2,0],lc2[idsub2,1])

 ccf = ss.correlate(x,echo2)
 nccf = np.shape(ccf)[0]
 tccf = (np.arange(nccf) - nccf/2)*dt











