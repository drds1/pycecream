### subroutine to call a light curve and construct a power spectrum from the data
# NOTE when specifying the w(iw) = wlo + dw*(iw-1), the -1 is very important, if you dont include it, the model will not fit the data very well, (see you 1st year report for an example of what happens if you leave off the - 1.

# snippet from fortran code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NW=ceiling(0.5*NT)
#whi=pi/dt
#wlo=twopi*1/tlen!*0.01
#dw=(whi-wlo)/(NW-1)
#wlo=dw
#NW=ceiling((whi-wlo)/dw)-2


#allocate(w(NW),dtau(NW),sigdtau(NW))
#do i=1,NW
#w(i)=wlo+dw*(i-1)
#enddo
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import numpy as np
import os
import pylab as plt
from myoptscal import *

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
    op = mypowerspecdef.mypowerspec(fname[i])
#    dat=hesfit(x,y,sig,p):
#    dat[0][0]

#this is for the plotting
    ax1=fig.add_subplot(nvert,nalong,0)             # plot the light curve and model fit
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel(r'f$_{\nu}$(t)')
    ax1.errorbar(t,x,sig,linestyle='None',color=col[i])
    ax1.plot(tmod,xop,color=col[i])
    
    ax1=fig.add_subplot(nvert,nalong,1)             # plot the amplitude spectrum
    ax1.set_xlabel('frequency (cyc/day)')
    ax1.set_ylabel(r'$Amp^2$')
#    ax1.set_yscale('log', basey=10)
    ax1.set_xscale('log', basex=10)
    ax1.errorbar(f[:nf],ampspeclog,sigampspeclog,linestyle='None',color=col[i])
    if (autoylim == 0):
        ax1.set_ylim(ylim[0],ylim[1])
    
    
    if (include_phase == 1):
        ax1=fig.add_subplot(nvert,nalong,2)             # plot the phase spectrum
        ax1.set_xlabel('frequency (cyc/day)')
        ax1.set_ylabel('phase')
        ax1.errorbar(f[:nf],phasespec,sigphasespec,linestyle='None',color=col[i])

ax1.legend(fname)

	# save the output
opdat = np.zeros((nf,3))
opdat[:,0]=np.log10(f[:nf])
opdat[:,1]=ampspeclog
opdat[:,2]=sigampspeclog
np.savetxt('mypowerspecop_'+fname[i],opdat)

os.chdir(pwd) #change back to initial directory    




fig.show()

#fig.savefig('mypowerspecfig.png')
       
    

#    ax1=fig.add_subplot(nalong,nvert,idxinc)
#    ax1.set_ylabel(r'$\psi$($\tau$ | $\lambda$)')
#    ax1.set_xlabel(r' $\tau$ (days)')
#    ax1.errorbar(t,x,sig,linestyle=None,color=col[i])


