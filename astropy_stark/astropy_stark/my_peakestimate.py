import numpy as np
import matplotlib.pylab as plt
#code to test weather I have the correct formula to predict the peak position of 
#a sinusoid with power-law-variable frequency and amplitude

n = np.arange(10)
f_slope = 1.0

f0 = 1.0
r0 = 1.0
A0 = 0.1
amp_slope = 1.0
wav0 = 1./f0

bg0 = 1.0
bg_slope = 2.0

rp = (n*wav0*r0**f_slope)**(1./(f_slope+1.))


rpmax = np.max(rp)
x = np.linspace(0,rpmax,1000)

amp = A0*(x/r0)**amp_slope
f = f0*(x/r0)**f_slope
bg = bg0*(x/r0)**bg_slope

y = amp * np.cos(2*np.pi * f  * x) + bg




fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x,y)


yp = np.interp(rp,x,y)


ax1.scatter(rp,yp,color='r')
#ax2 = fig.add_subplot(212)
#ax2.scatter(n,rp)


plt.show()













