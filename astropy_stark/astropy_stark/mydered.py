#insert fluxes and wavelengths
#deredden these to nu^1/3 and report back E(B-V) gaskell 2004

import numpy as np
import myextmag as em
import matplotlib.pylab as plt
import scipy.optimize as so
import myrandom as mr

wav = np.arange(2000,9000,1000.)
nw = np.shape(wav)[0]
wav0 = 4000.0
fnu0 = 10.0

ebmvtrue = 0.16
A0 = em.extmag_agn(wav0,1.0)[0]
Awav = em.extmag_agn(wav,1.0)
Awav = np.array(Awav)
fnu = fnu0*(wav/wav0)**(-1./3)*10**(-0.4*ebmvtrue*(Awav)) 


#fnu[0:3] = fnu[0:3] - 5.0

sigfnu = np.ones(nw)
fnu = fnu + mr.normdis(1,0.,sigfnu)

logfnu = np.log10(fnu)
logwav = np.log10(wav)

rln10 = np.log(10.)
siglog = 1./rln10/fnu * sigfnu





nw = wav.shape[0]
#wav0 = wav[nw/2]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! define function to fit!!!!!!!!!!!!!!!!!!!!!!!
#wav_wav0_log = log(wav/wav0)
def fnu_red_log(wav_wav0_log,fnu0_log,ebmv): 
 return(fnu0_log - 0.33333*wav_wav0_log - 0.4*ebmv*(Awav))
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


wav_wav0_log = np.log10(wav/wav0)

#now fit the reddening function (Get back ebmv).
fit_coef = so.curve_fit(fnu_red_log,wav_wav0_log,logfnu,sigma=siglog)
sig_coef=np.sqrt([fit_coef[1][0,0],fit_coef[1][1,1]])
Nres = 100

wavhi = wav[-1]*1.3
wavlo = wav[0]*0.7
wavmod = wavlo + 1.*np.arange(Nres)/Nres*(wavhi-wavlo)
logwavmod=np.log10(wavmod)

fnu0_log = fit_coef[0][0]
sigfnu0_log = sig_coef[0]
sigebmv     = sig_coef[1]
ebmv     = fit_coef[0][1]
Awavmod = np.array(em.extmag_agn(wavmod,1.0))
logfnumod = fnu0_log - 0.33333*(logwavmod - np.log10(wav0)) #- 0.4*ebmv*(Awavmod )
siglogfnu = np.sqrt(sigfnu0_log**2 + (0.4*(Awavmod-A0)*sigebmv)**2) 

rln10 = np.log(10.)
lfnumod10 = 10**logfnumod

fnumod = lfnumod10*(1. + 0.5*rln10*rln10*siglogfnu*siglogfnu)
sigfnumod = rln10*lfnumod10*siglogfnu
#


#plot the result

fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.errorbar(wav,fnu,sigfnu,linestyle='',color='r',linewidth=4)
fnumodlo = fnumod-sigfnumod
fnumodhi = fnumod+sigfnumod
ax1.plot(10**logwavmod,fnumodhi,color='r')
ax1.plot(10**logwavmod,fnumodlo,color='r')
ax1.fill_between(10**logwavmod,fnumodlo,fnumodhi,color='r',alpha=0.4)

l=[r'$ E(B-V)='+str(round(ebmv,2))+' \pm '+str(round(sigebmv,2))+'$']
ax1.plot(10**logwavmod,fnumod,color='r',label=l[0])
ax1.legend(l)

ax1.set_title('Mean and Variable Component Spectrum')   
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xlim(wavlo,wavhi)
ax1.set_xlabel(r'$ \mathrm{Wavelength} / \AA$')
ax1.set_ylabel(r'$f_{\nu} \left( \lambda \right)$')
xtl = np.arange(wavlo+(1000.-np.mod(wavlo,1000.)),wavhi-np.mod(wavhi,1000),2000.)#np.logspace(np.log10(wavlo),np.log10(wavhi),5)#[100,200,500,1000]
ax1.set_xticks(xtl)
xtls = ["%d" % f for f in xtl]
ax1.set_xticklabels(xtls)


plt.show()
 
 












