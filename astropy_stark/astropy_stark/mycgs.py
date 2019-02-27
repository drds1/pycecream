import numpy as np


##  function of a number of physical constant i

pi=np.pi
c=2.997925*10**10
G=6.673*10**-8
k=1.3806*10**-16
#R=8.314
#mu_0=4*pi*10**-7
#ep_0=8.85*10**-12
#w_c=2.8977686e-3
plank=6.626e-27


### astrophysics constant
Msun=2*10**33
ly=9.4605284*10**17
ps=3.258*ly
sb=5.670373*10**-5
H0=1.70*100000/10**6/ps
#au = 149597871000

##black hole scwarzschild radius
def r_sch(M):
	rs=2*G*M/c**2
	return(rs)
