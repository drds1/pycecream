import numpy as np
from pylab import *

# script to make a skewed gaussian.

n=2000
lo=-10.0
hi=10.0
a=4


def phi(x):
    phi=1/np.sqrt(2*np.pi) *e**(-x*x/2)
    return(phi)
    


x=np.arange(lo,hi,(hi-lo)/n)

    
y=2*phi(x)*phi(a*x)


plot(x,y)

show()



def gaus(x,x0,w,b,dx):
    y=e**(-np.log(2*(np.log(1.+2.*b*(x-x0)/dx)*1/b)**2))
    return(y)
    
    
x0=1.0
w=0.5
b=1.0
dx=0.1

y=gaus(x,x0,w,b,dx)