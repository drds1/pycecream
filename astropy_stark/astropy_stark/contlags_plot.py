#rsinc ds207@lucky2.st-andrews.ac.uk:/microlens/ds207/fortcode/cream_export/katelc_group_1febmerge_* ./sd_hilag_111117/
import numpy as np
import glob
import matplotlib.pylab as plt
import os
import mymedquart as mm
import myrandom as mr
import myedlum as me
from astropy import constants as con
from astropy import units as u
import scipy.optimize

recalc = 1
plotmoddat = 1
plotedrat = 1
#tplotlim = 
logscal = [1,0]
xpl = [[2,99],[0.1,99.9]]
#xpl = [[1.e-5,10000.],[0.,1.]]
#ypl=[[0.,10.],[0.,1.]]

loglbol=0

labsave = ['mdot','cosinc']
axlab   = [r'$M\dot{M}$ $(10^7 M_\odot^2\mathrm{yr}^{-1})$',r'$\cos i$']



#make plots
#load data (if just redooing plot can skip to here)

parsave_plot = np.loadtxt('parsave_plot.txt')
parsigsave_plot = np.loadtxt('parsigsave_plot.txt')
lagsave_plot = np.loadtxt('lagsave_plot.txt')
lagsigsave_plot = np.loadtxt('lagsigsave_plot.txt')

lagsaveth_plot = np.loadtxt('lagsaveth_plot.txt')
lagsigsaveth_plot = np.loadtxt('lagsavesigth_plot.txt')




cisqred_plot = np.loadtxt('cisqred_plot.txt')
#fvar_plot = np.loadtxt('fvar_plot.txt')
snrstatsave_plot = np.loadtxt('snrstatsave_plot.txt')
fvm2_plot =  np.loadtxt('xsnorm_plot.txt')
fvm_plot = np.sqrt(fvm2_plot)
sigxs_plot = np.loadtxt('xs_plot.txt')

fname = 'snr_kate.txt'
with open(fname) as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content] 
f.close()
snr_kate = [np.float(con.split()[1]) for con in content]


nsave,npar = np.shape(parsigsave_plot)


#if plotedrat on then load the black hole mass estimates for the sdss quasars, match the ids to the rm ids analysed by cream and calculate the eddington ratios
if (plotedrat == 1):
 os.system('python sdss_bhmasslist.py')
 datbh = np.loadtxt('sdss_bhmasslist.op')
 embhlog    = datbh[:,1]
 embhsiglog = datbh[:,2]
 rmid = datbh[:,0].astype(int)
 ummdot = parsave_plot[:,0]*1.e7
 ummdotsig = parsigsave_plot[:,0]*1.e7
 ummdotlog    = np.log10(ummdot)
 ummdotsiglog = 1./np.log(10)/ummdot * ummdotsig #np.log10(parsave_plot[:,0]*1.e7)

 #match rmid of cream runs to black hole masses in fits file
 targsave = np.loadtxt('targsave_op.txt')
 targload = datbh[:,0]
 idnow = np.array( [np.where(targload == tsnow)[0][0] for tsnow in targsave] )
 embhlog = embhlog[idnow]
 embhsiglog = embhsiglog[idnow]
 
 #calculate eddington ratios
 eta = 0.1
 a = np.array(me.edd2(embhlog, embhsiglog, ummdotlog, ummdotsiglog, eta)).T
 eddrat = a[:,2]
 sigeddrat = a[:,3]   
 
 
 #add a parameter to the parsave_plot column to instruct code to make extra plot
 a = np.ones((nsave,1))
 a[:,0] = a[:,0]*eddrat 
 parsave_plot = np.hstack((parsave_plot,a))
 a = np.ones((nsave,1))
 a[:,0] = a[:,0]*sigeddrat 
 parsigsave_plot = np.hstack((parsigsave_plot,a))
 labsave.append('eddrat')
 axlab.append(r'$L/L_\mathrm{edd}$')
 logscal.append(1)
 #xmed = np.percentile(eddrat,50)
 #xlo  = np.percentile(eddrat,1)
 #xhi  = np.percentile(eddrat,99)
 xpl.append([1,99])
 nsave,npar = np.shape(parsigsave_plot)


for ip in range(npar):
 fig = plt.figure()
 ax1 = fig.add_subplot(111)
 
 ynow = []
 xnow = []
 znow = []
 
 
 alphanow = 1/parsigsave_plot[:,ip]
 alphanow = alphanow/np.max(alphanow)
 for it in range(nsave):
  #xnow.append(cisqred_plot[it])
  ynow.append(snrstatsave_plot[it])
  xnow.append(parsave_plot[it,ip])
  znow.append(parsigsave_plot[it,ip])
  #znow.append(alphanow[it])
 
 xnow = np.array(xnow)
 ynow = np.array(ynow)
 znow = np.array(znow)
 idinc = np.arange(nsave)#np.where((xnow < xpl[ip][1]) &  (ynow < ypl[ip][1]))[0]
 xnow = xnow[idinc]
 ynow = ynow[idinc]
 znow = znow[idinc]
 nsnew = np.shape(xnow)[0] 
 #for it in range(nsnew): 
 #ax1.scatter(xnow,ynow,znow,color='b')
 ax1.scatter(xnow,ynow,color='b',s=3.0)
 #ax1.errorbar(xnow,ynow,xerr=znow,ls='',capsize = 10,color='k',alpha=0.3)

 #ax1.set_xlabel(r'$\chi^2/n$')
 ax1.set_xlabel(axlab[ip])
 ax1.set_ylabel('SNR = (Max flux - Min flux)/Median Error) \n averaged over all telscopes per filter')

 if (logscal[ip] == 1):
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  #if (xpl[ip][0] != -1.0):
  # print 'xlim log plot',labsave[ip]+' ',xpl[ip]
  # ax1.set_xlim(xpl[ip])
 #elif (xpl[ip][0] == 0 & xpl[ip][1] == 0):
 print 'fsfsfd',xpl[ip]
 xlo = np.percentile(parsave_plot[:,ip],xpl[ip][0])
 xhi = np.percentile(parsave_plot[:,ip],xpl[ip][1])
 
 ax1.set_xlim([xlo,xhi])

 #else:
 # ax1.set_xlim(xpl[ip])
 # ax1.set_ylim(ypl[ip])
 plt.savefig('fig_'+labsave[ip]+'.pdf')






#make plot of snrstat (max(flux) - min(flux) / median errorbar) vs i - g continuum lag
try:
 idlaghi = 1
 idlaglo = 0
 fig = plt.figure()
 ax1 = fig.add_subplot(111)
 for it in range(nsave):
  x = lagsave_plot[it,idlaghi] - lagsave_plot[it,idlaglo]
  y = snrstatsave_plot[it]
  sigx = np.sqrt( lagsigsave_plot[it,idlaghi]**2 + lagsigsave_plot[it,idlaglo]**2 )
  ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
  ax1.scatter(x,y,color='k',s=0.5)
 
 ylim = list(ax1.get_ylim())
 ax1.plot([0.,0.],ylim,color='k',linewidth=3)
 ax1.set_ylim(ylim)
 ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
 ax1.set_ylabel('SNR = (Max flux - Min flux)/Median Error) \n averaged over all telscopes per filter')
 plt.savefig('fig_snrplot.pdf')
except:
 print 'problem making lag plot fig_snrplot.pdf. Maybe these lcs were fitted with top hat contiuum lags'


#make plot of top hat continuum lags
idlaghi = 1
idlaglo = 0
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
xsave = []
ysave = []
for it in range(nsave):
 x = lagsaveth_plot[it]
 if (x[0] == 0):
  x = x[1] - x[0]
 elif (x[1] == 0):
  x = x[0] - x[1]

 
 y = snrstatsave_plot[it]
 xsave.append(x)
 ysave.append(y)
 sigx = np.sqrt( np.sum(lagsigsaveth_plot[it]**2)  )
 ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
 ax1.scatter(x,y,color='k',s=0.5)

xl = ax1.get_xlim()
ax2.hist(xsave,bins=np.linspace(xl[0],xl[1],60))
ylim = list(ax1.get_ylim())
ax1.plot([0.,0.],ylim,color='k',linewidth=3)
ax1.set_ylim(ylim)
ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
ax1.set_ylabel('SNR = (Max flux - Min flux)/Median Error) \n averaged over all telscopes per filter')
ax2.set_xlabel(r'$\langle \tau_i \rangle$ (days)')
fig.tight_layout()
plt.savefig('fig_lagthsnrplot.pdf')



#make plot of top hat continuum lags
idlaghi = 1
idlaglo = 0
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
xsave = []
ysave = []
for it in range(nsave):
 x = lagsaveth_plot[it]
 if (x[0] == 0):
  x = x[1] - x[0]
 elif (x[1] == 0):
  x = x[0] - x[1]
 
 y = snr_kate[it]
 xsave.append(x)
 ysave.append(y)
 sigx = np.sqrt( np.sum(lagsigsaveth_plot[it])**2  )
 ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
 ax1.scatter(x,y,color='k',s=0.5)

ylim = list(ax1.get_ylim())
ax1.plot([0.,0.],ylim,color='k',linewidth=3)
ax1.set_ylim(ylim)
ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
ax1.set_ylabel(r'SNR = RMS$/\sigma_\mathrm{RMS}$')
xl = ax1.get_xlim()
ax2.hist(xsave,bins=np.linspace(xl[0],xl[1],60))
ax2.set_xlabel(r'$\langle \tau_i \rangle$ (days)')
fig.tight_layout()
plt.savefig('fig_lagthrmsplot.pdf')

#make plot of snrstat (max(flux) - min(flux) / median errorbar) vs i - g continuum lag
idlaghi = 1
idlaglo = 0
fig = plt.figure()
ax1 = fig.add_subplot(111)
for it in range(nsave):
 x = lagsave_plot[it,idlaghi] - lagsave_plot[it,idlaglo]
 y = fvm_plot[it]
 sigx = np.sqrt( lagsigsave_plot[it,idlaghi]**2 + lagsigsave_plot[it,idlaglo]**2 )
 ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
 ax1.scatter(x,y,color='k',s=0.5)

ylim = list(ax1.get_ylim())
ax1.plot([0.,0.],ylim,color='k',linewidth=3)
ax1.set_ylim(ylim)
ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
ax1.set_ylabel(r'$F_\mathrm{var}$ averaged over all telscopes per filter')
plt.savefig('fig_fvarplot.pdf')



#make plot of snrstat (max(flux) - min(flux) / median errorbar) vs i - g continuum lag
idlaghi = 1
idlaglo = 0
fig = plt.figure()
ax1 = fig.add_subplot(111)
for it in range(nsave):
 x = lagsave_plot[it,idlaghi] - lagsave_plot[it,idlaglo]
 y = sigxs_plot[it]
 sigx = np.sqrt( lagsigsave_plot[it,idlaghi]**2 + lagsigsave_plot[it,idlaglo]**2 )
 ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
 ax1.scatter(x,y,color='k',s=0.5)

ylim = list(ax1.get_ylim())
ax1.plot([0.,0.],ylim,color='k',linewidth=3)
ax1.set_ylim(ylim)
ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
ax1.set_ylabel(r'$\sigma_\mathrm{XS}$ \n averaged over all telscopes per filter')
plt.savefig('fig_sigxs.pdf')


   
   
   
   
#make a plot of eddington ratio obtained from fits file bh mass bolometric correction
# vs cream value
eddrat_fits = datbh[idnow,10]
sigeddrat_fits = datbh[idnow,11]
ider = -1
eddrat_cream = parsave_plot[:,ider]
sigeddrat_cream = parsigsave_plot[:,ider]

fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.errorbar(eddrat_fits,eddrat_cream,xerr=sigeddrat_fits,yerr=sigeddrat_cream,ls='',marker='',color='k',alpha=0.2)
ax1.scatter(eddrat_fits,eddrat_cream)
#ax1.set_xlim(np.percentile(eddrat_fits,[1,99]))
#ax1.set_ylim(np.percentile(eddrat_cream,[1,90]))
ax1.set_xlim([0,1])
ax1.set_ylim([0,1])

xl  = ax1.get_xlim()
yl = ax1.get_ylim()
ax1.plot(xl,yl,color='k',ls='--')
ax1.set_ylabel(r'$L/L_\mathrm{Edd}$ (CREAM)')
ax1.set_xlabel(r'$L/L_\mathrm{Edd}$ (Bolometric Correction)')
plt.savefig('fig_eddrat_comp.pdf')





fig = plt.figure()
ax1=fig.add_subplot(111)

y = eddrat_cream/eddrat_fits
idinc = np.where(y < np.inf)[0]
x = snrstatsave_plot
sigy = y * np.sqrt( (sigeddrat_cream/eddrat_cream)**2 + (sigeddrat_fits/eddrat_fits)**2 )
ax1.errorbar(x[idinc],y[idinc],sigy[idinc],color='k',alpha=0.3,ls='')
ax1.scatter(x[idinc],y[idinc],color='k',s=0.5)
#ax1.set_ylim([0,2])
ax1.set_yscale('log')
#ax1.set_ylim(np.percentile(eddrat_cream,[10,90]))
ax1.set_ylabel(r'$L/L_\mathrm{Edd}$ (CREAM) / $L/L_\mathrm{Edd}$ (fits)')
ax1.set_xlabel('SNR = (Max flux - Min flux)/Median Error) \n averaged over all telscopes per filter')
plt.savefig('fig_eddratrat_snr.pdf')




#ummdot
#ummdotsig 
em = 10**embhlog
emsig = np.log(10)*em*embhsiglog
eta = 0.1
mdot = ummdot/em
sigmdot = mdot * np.sqrt((emsig/em)**2 + (ummdotsig/ummdot)**2 )
syr = 31557600.0
conv = eta * con.c.value**2 * con.M_sun.value / syr * (1* u.joule).to(1* u.erg).value
lbol = mdot*conv
siglbol = sigmdot*conv 
x = datbh[idnow,12]
sigx = datbh[idnow,13]
#now plot Lbol vs nufnu





idgood = np.where((np.abs(np.log10(x)) < np.inf) & (np.abs(np.log10(x)) > 0))[0]
x = x[idgood]
sigx = sigx[idgood]



if (loglbol == 0):
 xlog = x
 sigxlog = sigx#1./np.log(10)/x * sigx
 xlogorth = np.mean(x)#np.mean(xlog)
 ylog = lbol[idgood]#np.log10(lbol)[idgood]
 sigylog = siglbol[idgood]#1./np.log(10)/lbol[idgood] * siglbol[idgood]
else:
 xlog = np.log10(x)
 sigxlog = 1./np.log(10)/x * sigx
 xlogorth = np.mean(xlog)
 ylog = np.log10(lbol)[idgood]
 sigylog = 1./np.log(10)/lbol[idgood] * siglbol[idgood]

 
def func(x,a,b): 
 return(a+b*x)
fit_coef=scipy.optimize.curve_fit(func,xlog-xlogorth,ylog,sigma=sigylog,bounds=[[-np.inf,-np.inf],[np.inf,np.inf]])
npar = np.shape(fit_coef[1])[0]
sig_coef = np.sqrt([[fit_coef[1][i,i]] for i in range(npar)])[:,0]
fit_coef = fit_coef[0]
xres = np.linspace(0.0,xlog.max(),500)#np.linspace(xlog.min(),xlog.max(),500)
yres = fit_coef[0] + (xres - xlogorth)*fit_coef[1]
yrlo = fit_coef[0] - sig_coef[0] + (xres - xlogorth)*(fit_coef[1] - sig_coef[1])
yrhi = fit_coef[0] + sig_coef[0] + (xres - xlogorth)*(fit_coef[1] + sig_coef[1])


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(xlog,ylog,xerr=sigxlog,yerr=sigylog,color='k',alpha=0.3,ls='')
ax1.scatter(xlog,ylog,color='k',s=0.5)
ax1.plot(xres,yres)
ax1.fill_between(xres,yrlo,yrhi,alpha=0.1)
if (loglbol == 1):
 ax1.set_ylabel(r'log $\eta \frac{M\dot{M}_\mathrm{CREAM}} {M_\mathrm{BH}} c^2$ (erg $\mathrm{s^{-1}}$)')
 ax1.set_xlabel(r'log $\nu F_\nu$ (erg cm$^{-2}$ $\mathrm{s^{-1}}$)')
else:
 ax1.set_ylabel(r'$\eta \frac{M\dot{M}_\mathrm{CREAM}} {M_\mathrm{BH}} c^2$ (erg $\mathrm{s^{-1}}$)')
 ax1.set_xlabel(r'$\nu F_\nu$ (erg cm$^{-2}$ $\mathrm{s^{-1}}$)')

#ax1.set_yscale('log')
#ax1.set_xscale('log')
ax1.set_ylim([0,1.e46])
ax1.set_xlim([0,1.4e46])
tit = r'$L_\mathrm{BOL} = ('+np.str(np.round(fit_coef[0],2))+'\pm '+np.str(np.round(sig_coef[0],2))+')$ $'+' \\nu L_\\nu$ (erg cm$^{-2}$ $\mathrm{s^{-1}}$)'
ax1.set_title(tit)
plt.savefig('fig_lbovsnufnu.pdf')

