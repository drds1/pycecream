import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as mcfit






#generate a test polynomial function and add errors
ain = 2.
bin = 5.3
cin = 10.
x    = np.arange(100)
y    = ain*x**2 + bin*x + cin
sigy = np.ones(100)
y    = y + np.random.randn(100)*sigy #error of sigy about mean 0+ 0.



#define polynoimal function and fit
def func(x, a, b ,c):
 return(a*x**2 + b*x + c)
 

#fit function
popt, pcov = mcfit.curve_fit(func, x, y, sigma=sigy)
 


#output
xres = np.arange(0,100,0.1) 
yres = popt[0]*xres**2 + popt[1]*xres + popt[2] 

#sample probability space to generate errors






#plot the result
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(xres,yres)
ax1.errorbar(x,y,sigy,ls='')
plt.show()




