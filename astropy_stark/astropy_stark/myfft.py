#manually written code to obtain the fourier transform of a real input signal
#it also does the inverse fft too. Will hopefully code this to give me the uncertainties
# of the spectrum given error bars on the signal. numpy and scipy's code does not do this.
#Thie FT gives both plus and - frequencies. The amplitude shown on the (e.g +ve frequency end) plots are half 
#the true amplitude as both +ve and -ve terms contribute half to the true amplitude

from pylab import *
import numpy as np

twopi=2*np.pi



# subroutine to calculate the fourier transform and uncertainties of an input data set
# this will go into the routine
# INPUT: x[:N],y[:N] signal
# OUTPUT: ftreal,ftimag,ftabs, sigftreal,sigftimag,sigftabs the sin and cosines amplitudes,and absolute amplitude spectrum and uncertaintes of the signal
# freq[:N] the frequencies of the signal


def myft(x,y,sigy=0):
    N=y.shape[0]
    dt = (x[-1] - x[0])/(N -1)
    if (isinstance(sigy,list)):
     sigyin = sigy
    else:
     sigyin = np.zeros(N)
     
    ftreal=np.zeros(N)
    ftimag=np.zeros(N)
    ftabs=np.zeros(N)
    sigftreal=np.zeros(N)
    sigftimag=np.zeros(N)
    sigftabs=np.zeros(N)

# calculate the real, imaginary, absolue, and uncertainties in the amplitude spectrum    
    for ik in range(N):
    
        w = twopi*(ik)*(np.arange(N)) /N
        cw = np.cos(w)
        sw = np.sin(w)
        
        realik = np.sum( y[:] * cw ) /N
        realik2 = realik*realik
        sigrealik = np.sqrt( np.sum( (cw*sigyin)**2) )/N
        sigrealik2 = sigrealik*sigrealik
        
        imagik = -np.sum( y[:] * sw ) /N
        sigimagik = np.sqrt( np.sum( (sw*sigyin)**2) )/N
        imagik2 = imagik*imagik
        sigimagik2 = sigimagik*sigimagik
        
        ftreal[ik] = realik
        sigftreal[ik]=sigrealik
        
        ftimag[ik] = imagik
        sigftimag[ik]= sigimagik
        
        ftabs[ik] = np.sqrt(realik2 + imagik2)
        sigftabs[ik] = 2*np.sqrt( realik2*sigrealik2 + imagik2*sigimagik2 )

    freq = np.fft.fftfreq(N,d=(N-1.)/N*dt)
    
    return(ftreal,sigftreal,ftimag,sigftimag,ftabs,sigftabs,freq)

#dftreal=1.*ftreal
#dftimag=1.*ftimag
#NT = N    
#t=1.*x
#tlo=t[0]
#thi=t[-1]  


#! subroutine to calculate the inverse fourier transform of an input set of real,
#! and imaginary fourier coefficients and a desired time grid to evaluate the ift at.
# inputs tlo,thi: first and last data point on the time series, dt: desired precision of output signal  
# outputs: tmod (the modified time array), oprealplot (the function after inverse ft)
  
def myift(ftreal,ftimag,tlo,thi,dtmod=-1):    
    tlen = thi - tlo 
    if (dtmod == -1):
        dtmod = tlen/(ftreal.shape[0]-1)/10
    
    
    
    NT = ftreal.shape[0]
    dt = tlen/(NT-1) # the separation between
    #dtmod = dt/10
    #NTmod = np.int(tlen/dtmod +1 )
    tmod = np.arange(tlo,thi+dtmod,dtmod)
    NTmod = np.int(tmod.shape[0])
#    print NTmod
    opreal=[]
    Ntfft = NTmod
    dtfft = (NT-1.)/(NTmod-1)
    tfft = np.arange(0,NT-1+dtfft,dtfft)
#    print tfft.shape[0]
#tfft = np.arange(NT)
#t_op  = np.arange(tlo,thi+dtmod,dtmod)

    #the NT/tfft[-1] is introduced to fix a small offset with the test lightcurve after the ift. NOT sure why this happened
    freq = np.fft.fftfreq(NT,d=NT/tfft[-1]) # give me the frequencies of the fourier transofmr
    wfft = 2*twopi*freq

#    oprealplot=[]
#    for it in range(NT):
#        w = twopi*np.arange(NT)*it /NT
#        sumreal = np.sum( ftreal[:NT]*np.cos(w) ) - np.sum( ftimag[:NT]*np.sin(w) )
#        opreal.append(sumreal)
    
    oprealplot=np.zeros(NTmod)
    #print NTmod,tmod.shape[0]
    for it in range(NTmod):
        oprealplot[it] = oprealplot[it]+np.sum( ftreal*np.cos(wfft*tfft[it]) - ftimag*np.sin(wfft*tfft[it]) )

    return(tmod,oprealplot)

    

#tlo=10.0
#thi=20.0
#dt=1.0
#x=np.arange(tlo,thi+dt,dt)
#y=1*np.cos(twopi*x) + 1*np.sin(twopi*2*x) + 4
#nx=x.size
#sigy=np.ones(nx)
#
#dat = myfft(x,y,sigy)
#realamp = dat[0]
#sigrealamp = dat[1]
#imagamp = dat[2]
#sigimagamp = dat[3]
#
#absamp = dat[4]
#sigabsamp=dat[5]
#freq=dat[6]
#
#errorbar(freq,absamp,sigabsamp)
#show()
#
##plot(y)
##plot(opreal)
##plot(tfft,oprealplot)
#
##show()
#    
#dat= myift(realamp,imagamp,tlo,thi)    
#plot(dat[0],dat[1])
#plot(x,dat[2])
#plot(x,y)
#show()
#

#    