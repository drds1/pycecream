import astropy_stark.myfake as mf
import matplotlib.pylab as plt
lc1 = mf.myfake(
    [-1.0,-1.0],
    [50.0,100.],
    [1.0,2.0],
    sdforcein=[1.0,1.0],
    meanforcein = [0.0,0.0],
    thcent = 20.0,
    iseed = 12345
)['echo light curves']

lc2= mf.myfake(
    [-1.0],
    [50.0],
    [1.0],
    thcent = 5.0,
    sdforcein=[2.0],
    meanforcein = [5.0],
    iseed = 12345
)['echo light curves']
dat = lc1+lc2
#plot the light curves here to demonstrate
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in range(len(dat)):
    ax1.errorbar(dat[i][:,0],dat[i][:,1],dat[i][:,2],ls=None)
plt.show()
