import numpy as np
import mydiptest_2 as md2
import myrandom as mr
import matplotlib.pylab as plt

nboot = 10
cent = np.arange(0,10,1)

pplot = []
dplot = []

xpdf1 = mr.normdis(5000,0,1)

#try with gaussian
ncen = np.shape(cent)[0]
for ib in range(ncen):
 xpdf2 = mr.normdis(5000,cent[ib],1)
 xpdf  = np.concatenate((xpdf1,xpdf2))
 print 'before dtest',ib
 dtest = md2.DipTest(xpdf)#md2.DipTestSig(xpdf,1)
 pplot.append(dtest[0])
 dplot.append(dtest[0])


fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(cent,pplot)
ax1.set_ylabel('p value')

ax1 = fig.add_subplot(212)
ax1.plot(cent,dplot)
ax1.set_ylabel('dip test')


plt.savefig('fig_diptesttest.pdf')

