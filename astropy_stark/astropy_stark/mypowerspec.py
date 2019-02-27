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
    dat =np.loadtxt(fname[i])
    t   =dat[:,0]
    thi =t.max()
    tlo = t.min()
    tlen = thi - tlo
    x   =dat[:,1]
    sig =dat[:,2]
#    sig = sig*0 + 0.1
    nx  = x.shape[0]
    nxmod = nx*10
    dtmod = (t[-1]-t[0])/(nxmod-1)
    tmod=np.arange(t[0],t[-1]+dtmod, dtmod)
    nxmod=tmod.shape[0]
    
    xop = np.zeros(nxmod)
    sop=[]
    sigsop=[]
    cop=[]
    sigcop=[]
    ampspec=[]
    sigampspec=[]
    phasespec=[]
    sigphasespec=[]
    
#specify the lower and upper frequency limits
    dt=t[1:] - t[:-1]
    dt=dt[np.nonzero(dt)] #problem fixed if 2 datapoints taken at same time
    dtmin = dt.min() #this is the interval between points

# this fixes data points that are super close together making the nyquist frequency higher than can be computed in a reasonable time.    
    if (dtmin < properdtmin):
        print 'Time resolution too high, using average time resolution...'
        dtmin = (tlen/(nx-1))
        if (dtmin < properdtmin):
            print 'The time resolution on this lightcurve is better than sense! I dont like it and will use',properdtmin,'instead.'
            dtmin = properdtmin
    


    nf  = np.int( np.ceil(0.5*nx) )
    fhi = 0.5/dtmin              
    flo = 1/(t[-1]-t[0])
    df = 0.5*(fhi - flo)/(nf-1)
    f = np.arange(flo,fhi+df,df)
    nf = f.shape[0]
    
    w = 2*np.pi*f
# subtract the mean from the data (to avoid worry about the constant offset
    xmean = np.median(x)
    xinp = x - xmean
    xop = np.ones(nxmod)*xmean
    
    
    for i_f in range(nf):
        optscalop = myoptscal(t,xinp,sig,f[i_f])
        s    = optscalop[0]
        sigs = optscalop[1]
        c    = optscalop[2]
        sigc = optscalop[3]
        sop.append(s)
        sigsop.append(sigs)
        cop.append(c)
        sigcop.append(sigc)
        

             
#these bits are updated every i_f
        fit = s*np.sin(w[i_f]*t) + c*np.cos(w[i_f]*t)
        xinp = xinp - fit
        fitmod = s*np.sin(w[i_f]*tmod) + c*np.cos(w[i_f]*tmod)
        
#        print fitmod.mean(),i_f,'test',np.std(fitmod)
# if you fit to too high frequencies, the model goes to crazy high values, introduce a check for this below.
        if (np.abs( np.std(fitmod) ) > 100*xmean):
            print 'The fit has gone crazy! you have tried to fit to too high a frequency... stopping here.'
            nf = len(ampspec)
            break

# amplitude spectrum        
        sigssq = 2*s*sigs
        sigcsq = 2*c*sigc
        a =s*s+c*c
        siga = np.sqrt(sigssq*sigssq + sigcsq*sigcsq)
        ampspec.append(np.sqrt(a))
        sigampspec.append( 0.5/np.sqrt(a)* siga )

# phase spectrum
        soverc    = s/c
        sigsoverc = soverc * np.sqrt( (sigs/s)**2 + (sigc/c)**2 )
        phasespec.append( np.arctan(-1.*soverc) )
        sigphasespec.append( np.abs( sigsoverc/(1. + soverc*soverc) ) ) 

# update the plot model
        xop = xop + fitmod
        


# fit a broken power law to the amplitude spectrum
    ampspeclog = np.log10(np.array(ampspec))
    sigampspeclog = np.abs( np.array(sigampspec)/(np.array(ampspec)*loge10 ) )
    
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


