from sklearn import decomposition
import numpy as np





#generate n fake light curves with different levels of correlation to the 'flows'
r = [0.8,0.5,0.2,0,0,0,0,0]

import mylcgen as mlc
flows = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=1100,dt=30.0,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
flowmean = np.mean(flows[:,1])
flowrms  = np.std(flows[:,1])


nflow = np.shape(flows[:,0])[0]
x_op = [flows[:,0]]

#the constituent signals will be a combination of gaussian noise and the original signal
#The proportion of the original signal represented by the time series is represented by
#r. r = 1 means the signal has no noise and entirely represents the original time series
#r = 0 means the signal is pure noise.
for rnow in r:
 noise = np.random.randn(nflow)*flowrms*np.sqrt(1.-rnow)
 x     = noise + (flows[:,1]-flowmean)*np.sqrt(rnow)/flowrms
 x_op.append(x)





#can we retrieve the explained variance information
var = np.array(pca.explained_variance_)
varnorm = var/var[0]


#!!!!!#!!!!!!!
x_op = np.array(x_op).T
x_op = (x_op - np.mean(x_op,axis=0))/np.std(x_op,axis=0)
nt,ndim = np.shape(x_op)
pca = decomposition.PCA(n_components=ndim)
pca.fit(x_op)
X = pca.transform(x_op)

