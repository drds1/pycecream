import kdestats as kde
import numpy as np
import matplotlib.pylab as plt
from scipy import stats
from scipy import interpolate


#test script how to identify confidence regions from a kde useful for outlier code!

covmat = [[1., 11.5], [1.5, 4.]]

#fop_pca = np.random.multivariate_normal([-3, 5], covmat, 10000)

a = np.random.multivariate_normal([-3, 5], covmat, 5000)
b = np.random.multivariate_normal([0, 0], covmat, 5000)
fop_pca = np.vstack((a,b))



xmin,xmax = np.min(fop_pca[:,0]),np.max(fop_pca[:,0])
ymin,ymax = np.min(fop_pca[:,1]),np.max(fop_pca[:,1])
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = fop_pca.T
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)



#kdehist = kde.kdehist2(xy[:,0], xy[:,1], [30, 30])
clevels = kde.confmap(Z.T, [.9973,.9545,.6827])



fig = plt.figure()
ax1 = fig.add_subplot(111)

xmod = X[:,0]
ymod = Y[0,:]
nxmod = np.shape(xmod)[0]
#ax1.contour(xmod,ymod,Z,clevels)
ax1.scatter(fop_pca[:,0],fop_pca[:,1])

#f = interpolate.interp2d(xmod, ymod, Z, kind='cubic')
#zdat = f(values[0,:],values[1,:])
#find Z values of actual scatter points
idx_x = np.array(np.interp(values[0,:],xmod,np.arange(nxmod)),dtype='int')
idx_y = np.array(np.interp(values[1,:],ymod,np.arange(nxmod)),dtype='int')
zdat = Z[idx_x,idx_y]


for i in range(len(clevels)):
 idlike = np.where(zdat > clevels[i])
 ax1.scatter(fop_pca[idlike,0],fop_pca[idlike,1])
plt.show()



