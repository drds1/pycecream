import numpy as np
import matplotlib.pylab as plt
import myvarlc as mv

file = '/Users/ds207/Google Drive/cont_222_thlag/katelc_group_1febmerge_8/rm369_gimerge/rm369_g_bok_G10_ccd4.dat'
#load light curve
dat = np.loadtxt(file)




#run script
time = dat[:,0]
flux = dat[:,1]
err  = dat[:,2]


a = mv.linvarfit(time,flux,err,diag=0)
print 'linvarfitop',a




#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.errorbar(time,flux,err,ls='')
#plt.show()