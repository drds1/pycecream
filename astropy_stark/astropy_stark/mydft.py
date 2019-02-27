# code to load and input light curve, and calculate the fourier transform and power spectrum, plotting the output on the plot.

import numpy as np
import os
import pylab as plt
from myoptscal import *
import scipy
import myfft

pwd = os.getcwd()
print 'Current directory :- ', pwd
dir = raw_input('Please input directory :- ')
os.chdir(dir)

nfiles = np.int(raw_input('Please input number of LC to analyse :- '))
fname=[]
for i in range(nfiles):
    print 'Please enter file',i+1,'of',nfiles,' :- '
    file = raw_input()
    fname.append(file)



#### hardwired parms

twopi=2*np.pi
nalong=1
include_phase = 0 # 0 or 1 depending if you want to include the phase spectrum in the plot
if (include_phase == 0):
    nvert = 2
else:
    nvert = 3

# set up the plot window
properdtmin = 0.2
fontsize=12

autoylim = 0
ylim=[-6,3]

col=['r','b','g','c','o','p','k','y']
fig=plt.figure()
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.3)
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fontsize}

plt.rc('font', **font)
loge10 = np.log(10)

##### end of hardwired parms

for i in range(nfiles):
    dat =np.loadtxt(fname[i])
    t   =dat[:,0]
    x   =dat[:,1]
#    x=x-x.mean()
    sig =dat[:,2]
    

    nt  = t.shape[0]
    
    tfft=1.*np.arange(nt)
    ntmod = nt*10
    tlo=t[0]
    thi=t[-1]
    tlen = thi-tlo
    dt = tlen/(nt-1)
    dtmodplot = tlen/(ntmod-1)
    tmodplot = np.arange(t[0],t[-1],dtmodplot)  
    ntmod=tmodplot.shape[0]
    xmodplot = np.zeros(ntmod)
    
    #nt=tlen1/dt1 = tlen2/dt2

    ntfftmod = 1.*ntmod
    dtfftmod = (nt-1)/tlen * dtmodplot
    tfftmod=np.arange(tfft[0],tfft[-1]+dtfftmod, dtfftmod)

    
#    ntmodplot = (tmod[-1]-tmod[0])/dtmod + 1
    
    
#    tmodplot = 
#    ntmod=tmod.shape[0]

#calculate the fourier transform of the data must be evenly spaced, if not evenly spaced, interpolate.
    
    
    ft = np.fft.fft(x) / nt
    
    
    ftreal = np.real(ft)
    ftimag = np.imag(ft)
    freq = np.fft.fftfreq(nt,d=dt)
    w = twopi*freq
    ftabs = np.abs(ft)
    
    for it in range(ntmod):
        xmodplot[it] = xmodplot[it]  + np.sum( ftreal*np.cos(w*tfftmod[it]) - ftimag*np.sin(w*tfftmod[it]) )
    
#    xmodplot = xmodplot
 
 
 
    ft = myfft.myft(t,x,sig)
    ftabs = ft[4]
    sigftabs=ft[5]
    freq = ft[6]
    
    ift = myfft.myift(ftreal,ftimag,tlo,thi)   
    tmodplot = ift[0]
    xmodplot = ift[1]

#this is for the plotting
    ax1=fig.add_subplot(nvert,nalong,0)             # plot the light curve and model fit
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel(r'f$_{\nu}$(t)')
    ax1.errorbar(t,x,sig,linestyle='None',color=col[i])
    ax1.plot(tmodplot,xmodplot,color=col[i])
   
    ax1=fig.add_subplot(nvert,nalong,1)             # plot the amplitude spectrum
    ax1.set_xlabel('frequency (cyc/day)')
    ax1.set_ylabel(r'$Amp^2$')
#    ax1.set_yscale('log', basey=10)
    ax1.set_xscale('log', basex=10)
    if (autoylim == 0):
        ax1.set_ylim(ylim[0],ylim[1])
    ax1.errorbar(freq,np.log10(ftabs),sigftabs,linestyle='None',marker='x',color=col[i])

ax1.legend(fname)

os.chdir(pwd) #change back to initial directory    

fig.show()
#fig.savefig('mydftfig.png')
       
    

#    ax1=fig.add_subplot(nalong,nvert,idxinc)
#    ax1.set_ylabel(r'$\psi$($\tau$ | $\lambda$)')
#    ax1.set_xlabel(r' $\tau$ (days)')
#    ax1.errorbar(t,x,sig,linestyle=None,color=col[i])
