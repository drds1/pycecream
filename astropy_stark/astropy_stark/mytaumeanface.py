## Python script to calculate the mean time delay for a face on disk
## this is not the mean delay. Just the delay where x=hc/(lamkT) =1
import numpy as np

G=6.673e-11
h=6.626e-34
k=1.38e-23
c=2.9979e8
pi=np.pi
msun=2e30
year=3600*24.*365

def taumeanface(lam,mmdot):
    
    x=lam*1e-10
    y=mmdot*msun*msun/year
    tdel=((45*G*y*x**4)/(16*pi**6 * c**5 *h))**(1./3)
    return(tdel)

