import numpy as np
from pylab import *



G=6.673e-11
h=6.626e-34
k=1.38e-23
c=2.9979e8
pi=np.pi
msun=2e30
year=3600*24.*365
day=year/365
def taumeanface(lam,mmdot):
    
    x=lam*1e-10
    y=mmdot*msun*msun/year
    tdel=((45*G*y*x**4)/(16*pi**6 * c**5 *h))**(1./3)
    return(tdel/day)



def dlflat(mmdot,wav,cosi,fvb,fvf):

    epsilon=fvf/fvb
    b=(1.-epsilon)/(1.-epsilon**1.5)
    
    taucent=taumeanface(mmdot,wav)
    print mmdot,taucent,wav
    
    dlf = 6.3 * taucent * (wav/1.e4)**(-1.5) * 1/sqrt(fvb/cosi) * b*b

    return(dlf)




dat=np.loadtxt('outputpars.dat')
nwav=6
ncol=dat[:,0].shape[0]
wav=np.zeros(nwav)
fvb=np.zeros(nwav)
fvf=np.zeros(nwav)
mmdot=dat[:,2]

wav[:]=1367.20000, 3556.50928, 4702.48389, 6175.57031, 7489.96973, 8946.70312
fvb[:]=3.4,3.6,1.14,1.18,1.05,1.20
fvf[:]= 2.2,2.4,0.9,1.0,0.90,0.95

dl=np.zeros((ncol,nwav))
cosi=1.0
for i in range(ncol):
    for i2 in range(nwav):
         dl[i,i2]=dlflat(mmdot[i],wav[i2],cosi,fvb[i2],fvf[i2])
         
             
    
