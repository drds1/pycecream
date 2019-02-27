#rsinc ds207@lucky2.st-andrews.ac.uk:/microlens/ds207/fortcode/cream_export/katelc_group_1febmerge_* ./sd_hilag_111117/
import numpy as np
import glob
import matplotlib.pylab as plt
import os
import mymedquart as mm
import myrandom as mr
import myvarlc as mvar

recalc = 1
plotmoddat = 0
contth = 1
#if contth = 1 and detect ncontfail lags (registers as failed) just in case acidentally set too many lags varying and not kept a reference lag
ncontfail = 3
#tplotlim = 
#if idburnin = 0.66 then use only final third of points to calculate uncertainties and parameters
idburnin = 0.66
labsave = ['mdot','cosinc']
axlab   = [r'$\dot{M}$ $(M_\odot^2\mathrm{yr}^{-1})$',r'$\cos i$']


pwd = os.getcwd()

if (recalc == 1):
 dmain = '/microlens/ds207/fortcode/cream_export/cream_python/katelc_group_1febmerge_*'
 #'/Users/ds207/Google Drive/cont_222_thlag/katelc_group_1febmerge_1'#'/Users/ds207/Google Drive/sd_hilag_111117/katelc_group_1febmerge_*'#'/Users/ds207/Google Drive/cont_lag_tophat_241117/katelc_group_1febmerge_*'#'/Users/ds207/Google Drive/sd_hilag_111117/katelc_group_1febmerge_*'#'/Users/ds207/Documents/standrews/sta/lightcurves/sd_hilag_111117/katelc_group_1febmerge_*'
 
 nits = 1000
 
 
 dsubs   = glob.glob(dmain)
 targ    = []
 targ_good = []
 cisqred = []
 snrstatsave = []
 fvm = []
 sigxs = []
 #cisave      = []
 #cisigsave   = []
 parsave     = []
 parsigsave  = []
 
 failsave = []
 #mdotsave    = []
 #mdotsigsave = []
 lagsave     = []
 siglagsave = []
 targsave = []
 cadmed = []
 frac2tot_plot = []
 dtcombmedsave_plot = []
 dtcombmed_sig_plot = []
 lagsaveth = []
 siglagsaveth = []
 snr_kate = []
 
 
 for dirnow in dsubs:
  os.chdir(dirnow)
  
  targnow = glob.glob('rm*')
 
  idrej = []
  ntn = len(targnow)
  for it in range(ntn):
   os.chdir(dirnow)
   tnow = targnow[it]
   
   os.chdir(tnow)
   
   try:
    #figure out how many light curves are at the same wavelength
    wav = []
    with open('creamnames.dat') as f:
     parline = f.readlines()
     parline = [x.strip() for x in parline]
    f.close()
    
 
    nlctemp = len(parline) 
    for i in range(nlctemp):
     pl = parline[i].split()
     wav.append(np.float(pl[-1]))
    
    wavu = np.unique(wav)
    nwavu = np.shape(wavu)[0]
    
    try:
     opdir = sorted(glob.glob('output_2*'))[-1]
    except:
     idrej.append(it)
     continue
    
    
    
    dnow = np.loadtxt(opdir+'/outputpars.dat')
    nitsnow = np.shape(dnow[:,0])[0]
    idit = int(nitsnow*idburnin)
    
    mdotnow = dnow[idit:,2]
    cincnow = dnow[idit:,3]
    
    a = mm.mymedquart(mdotnow,conlev=0.68)
    mdotmed = a[0]
    mdotlo  = a[1]
    mdothi  = a[2]
    
    
    a = mm.mymedquart(cincnow,conlev=0.68)
    cincmed = a[0]
    cinclo  = a[1]
    cinchi  = a[2]
    
    
    
    #if we have line lags load those too
    try:
     tfline = np.loadtxt(opdir+'/outputpars_th.dat')[idit:,:]
     linelag = 1
    except:
     linelag = 0 
    
    
    
    tfmod = np.loadtxt(opdir+'/plots/modeltf.dat')
    tfmod_sig = np.loadtxt(opdir+'/plots/modeltf_sig.dat')
    taumod = tfmod[:,0]
    tfmod = tfmod[:,1:]
    tfmod_sig = tfmod_sig[:,1:]
    ntau = np.shape(taumod)[0]
    
    wavnow = -22.0
    nlc= np.shape(tfmod)[1]
    #save the lags and use monte carlo techniques to retrieve error on lags
    
    lsnow = []
    lssignow = []
    print ''
    print 'target',tnow
    idnon0 = 0
 
 
    print 'hererre contth',contth     
    if (contth == 1):
     print 'a cont'
     thinfo = np.percentile(tfline[:,:nlc],[50,16,84],axis=0)
     #list the indicees of different mean lags (as same filter tied together), and non zero lags 
     print 'a2 cont'
     idu = np.unique(thinfo[0,:],return_index=True)[1]
     idus = np.sort(idu)
     #print 'b cont'
     if (np.shape(idu)[0] == ncontfail):
      print 'failsave',failsave
      failsave.append(tnow+' lag lc label prob')
      continue
     
     #print 'checking lagsave'
     #print thinfo[0,:]
     #print wav
     #print tnow
     #print idu
     #print idus
     #raw_input()
       
     lagsaveth.append(thinfo[0,idus])
     siglagsaveth.append((thinfo[2,idus] - thinfo[1,idus])/2)
     print 'lags',thinfo[0,idu]
     
     #failed if ncontfail lags
  
 
    if (contth ==0):
     for ilc in range(nlc):
      wn = wav[ilc]
      print 'wav wavnow',wn,wavnow
      #skip if shared wavelength
      if (wn == wavnow and ilc != idnon0):
       continue 
      
      #line
      elif (wn == -1):
       if (ilc == 0):
        tfeg = tfline[idit,:]
        idnon0 = np.where(tfeg != 0)[0][0]
       tmline = tfline[idit:,ilc]
       a = np.percentile(tmline,[50,16,84])
       print 'lag',a
      #continuum 
      else:    
       taumean = []
       print 'monte carloing lags...'
       for iteration in range(nits):
        psinew = mr.normdis(ntau,tfmod[:,ilc],tfmod_sig[:,ilc])
        taumean.append( np.sum(taumod*psinew)/np.sum(psinew) )
       tm = np.array(taumean)
       a = mm.mymedquart(tm,conlev=0.68)
  
      
      
      lsnow.append(a[0])
      lssignow.append( (a[2] - a[1])/2 )
      wavnow = wn
     
    
    #if we are stepping continuum together using top hat lags do this now
    
    
    
     if (len(lsnow) > nwavu):
      lsnow = lsnow[:nwavu]
     if (len(lssignow) > nwavu):
      lssignow = lssignow[:nwavu]
  
     lagsave.append(lsnow)
     siglagsave.append(lssignow)  
    
    
    
      
    print 'Done with monte carloing lags'
    
    
    
     
    
    #load model and data and plt
    model = np.loadtxt(opdir+'/plots/modellc.dat')
    ntmod = np.shape(model[:,0])[0]
    
    tmod = model[:,0]
    xmod = model[:,1:]
    
    fdat = sorted(glob.glob(opdir+'/plots/data_echo_dat_*'))
    cisqsum = 0
    nsum    = 0
    snrstatsum = 0
    
    if (plotmoddat == 1):
     fig = plt.figure()
     idpl = np.where(xmod[:,0] > 0.0)[0]
     tplotlim = [tmod[idpl[0]],tmod[idpl[-1]]]
    
    fvar = []
    sigfvar = []
    xsnorm=[]
    frac2tot = []
    #ntot = 0
    #top = 0.d0
    #bot = 0.d0
    timeu = [[] for i in range(nwavu)]
    for ilc in range(nlc):
     wavnow = wav[ilc]
     idu = np.where(wavu == wavnow)[0][0]
     
     dat = np.loadtxt(fdat[ilc])
     ndat = np.shape(dat[:,0])[0]
     tdat   = dat[:,0]
     timeu[idu].append(tdat)
     
     xdat   = dat[:,1]
     sigdat = dat[:,2]
     xm = np.mean(xdat)
     #ntot = ntot + ndat
     
     #calculate fvar
     print 'calculate fvar'
     variance = np.sum((xdat - xm)**2)/(ndat-1)
     sigmean2 = np.sum(sigdat**2)/ndat
     sigxsnow = variance - sigmean2
     fvarnow = sigxsnow/xm**2 #np.sqrt( np.sum((xdat - xm)**2 - sigdat**2 )/ndat ) / xm
     #print 'fvarnow',fvarnow
     frac2totnow = sigxsnow/variance
     frac2tot.append(frac2totnow)     
     print 'here -2'
     fvar.append(fvarnow)
     xsnorm.append(sigxsnow)
     #now calculate kate snr 
     b = fdat[ilc]+' '+ np.str(mvar.linvarfit(tdat,xdat,sigdat)[1])
     snr_kate.append(b+'\n')
     
     #sdm2=np.mean(sigdat**2)
     #print 'sdm2',sdm2,ndat,fvarnow,np.sum((xdat - xm)**2 - sigdat**2)
     #print (np.sqrt(sdm2/ndat) / xm )**2,'asdfad'
     #print (np.sqrt(1/2./ndat)*sdm2/xm**2/fvarnow)**2,'asdasd'
     fvarsignow = (np.sqrt(1/2./ndat)*sigmean2/xm**2/np.sqrt(fvarnow))**2 - (np.sqrt(sigmean2/ndat) / xm )**2
     #sigfvar.append(fvarsignow)
     print 'variance,mean square error',variance,sigmean2
     print 'normxs,xs out',fvarnow,sigxsnow
     print 'Fvar,sig Fvar out',np.sqrt(fvarnow),np.sqrt(fvarsignow)
     print 'frac2tot', frac2totnow
     #print 'i am here',ilc,nlc
     if (plotmoddat == 1):
      ax1 = fig.add_subplot(nlc,1,ilc+1)
      ax1.errorbar(tdat,xdat,sigdat,ls='',color='k',label=None)
      #ax1.plot(tmod,xmod[:,ilc],label=np.str(wav[ilc]))
      ax1.plot(tmod,xmod[:,ilc],label=None)
      xdm = np.mean(xdat)
      xsd = np.std(xdat)
      ax1.set_ylim([xdm-3.5*xsd,xdm+3.5*xsd])
      ax1.set_xlim(tplotlim)
      ax1.set_ylabel(r'flux $'+np.str(int(wav[ilc]))+' \AA$')
      #ax1.legend()
      if (ilc == nlc-1):
       fig.tight_layout()
       plt.savefig('creamplot.pdf')
     
     #print 'i am here',ilc 
     snrstatnow = (np.max(xdat) - np.min(xdat))/np.median(sigdat)
     snrstatsum = snrstatnow + snrstatsum
     xmitp = np.interp(tdat,tmod,xmod[:,ilc])
     cisqsum = cisqsum + np.sum(((xmitp - xdat)/sigdat)**2)
     nsum = nsum + ndat
    
    #print 'here -1'
    frac2tot_plot.append(np.mean(frac2tot))
    #calculate the cadence of each light curve
    dtcombmed_save    = []
    dtcombmedsig_save = []
    #print 'here 0 ',nwavu
    
    for i in range(nwavu):
     tcp = []
     tcp = tcp + [list(timeu[i][it]) for it in range(len(timeu[i]))]
     tcp = np.concatenate(tcp)
     tcomb = np.sort(tcp)
     #print 'here 1'
     dtcomb = tcomb[1:] - tcomb[:-1]
     #print 'here 2'
     tcombmed = np.percentile(dtcomb,[50,16,84])
     ntot = np.shape(tcp)
     tcpmax = np.max(tcp)
     tcpmin = np.min(tcp)
     dtcombmed_sig = (tcombmed[2] - tcombmed[1])/2
     dtcombmed = (tcpmax - tcpmin)/ntot#tcombmed[0]
     dtcombmed_save.append(dtcombmed)
     dtcombmedsig_save.append(dtcombmed_sig)
     #print i,nwavu,'here 3', dtcombmed_save
     #raw_input()
     
    #print 'i am here'
    dtcombmedsave_plot.append(dtcombmed_save)
    dtcombmed_sig_plot.append(dtcombmedsig_save)
    snrstatsave.append(snrstatsum/nlc) 
    parsave.append([mdotmed,cincmed])
    parsigsave.append([(mdothi - mdotlo)/2,(cinchi - cinclo)/2])  
    #mdotsave.append(mdotmed)
    #mdotsigsave.append((mdothi - mdotlo)/2)
    #cisave.append(cincmed)
    #cisigsave.append((cinchi - cinclo)/2)
    cisqred.append(cisqsum/nsum)
    fvarmean = np.mean(fvar)
    sigxsmean = np.mean(xsnorm)
    #print 'fvarmean info start'
    fvm.append(fvarmean)
    sigxs.append(sigxsmean)
    #print 'normxs info',fvm
    #print 'xs info',sigxs
    targ = targ + [targnow[it] for it in range(ntn) if it not in idrej ] 
    targsave.append(tnow)
   except:
    failsave.append(tnow)
 
   print ''
   print ''
 
 os.chdir(pwd)
 
 
 
 #save data 
 f = open('contlags_fail.txt','w')
 for fnow in failsave:
  f.write(fnow+'\n')
 f.close()
 
 parsave_plot = np.array(parsave)
 np.savetxt('parsave_plot.txt',parsave_plot)
 
 parsigsave_plot = np.array(parsigsave)
 np.savetxt('parsigsave_plot.txt',parsigsave_plot)
 
 lagsave_plot = np.array(lagsave)
 np.savetxt('lagsave_plot.txt',lagsave_plot)
 
 lagsigsave_plot = np.array(siglagsave)
 np.savetxt('lagsigsave_plot.txt',lagsigsave_plot)
 
 cisqred_plot = np.array(cisqred)
 np.savetxt('cisqred_plot.txt',cisqred_plot)
 
 snrstatsave_plot = np.array(snrstatsave)
 np.savetxt('snrstatsave_plot.txt',snrstatsave_plot)
 
 frac2tot_plot = np.array(frac2tot_plot)
 np.savetxt('frac2tot_plot.txt',frac2tot_plot)
  
 
 f = open('snr_kate.txt','w')
 for fn in snr_kate:
  f.write(fn)
 f.close()
 
 fvm_plot = np.array(fvm)
 np.savetxt('xsnorm_plot.txt',fvm_plot)

 
 sigxs_plot = np.array(sigxs)
 np.savetxt('xs_plot.txt',sigxs_plot)
 
 
 dtcombmedsave_plot = np.array(dtcombmedsave_plot)
 np.savetxt('cadence_plot.txt',dtcombmedsave_plot)
 
 dtcombmed_sig_plot = np.array(dtcombmed_sig_plot)
 np.savetxt('dtcombmed_sig_plot.txt',dtcombmed_sig_plot)

 
 lagsaveth_plot = np.array(lagsaveth)
 np.savetxt('lagsaveth_plot.txt',lagsaveth_plot)
 
 lagsavesigth_plot = np.array(siglagsaveth)
 np.savetxt('lagsavesigth_plot.txt',lagsavesigth_plot)
 

 
 
#extract numbers from targsave array
ts = [int(t[2]+t[3]+t[4]) for t in targsave]
ts = np.array(ts)

np.savetxt('targsave_op.txt',ts)

#nsave = len(parsave)
#npar = len(parsave[0])
#f = open('contlags_op.dat','w')
#f.write('target cisq/n '+np.str(npar)+' prameters and errors '+np.str(nlc)+' lags and errors (one for each light curve)\n')
#
#for i in range(nsave):
# a = ''
# for ip in range(npar):
#  a = a+np.str(parsave[i][ip])+' '+np.str(parsigsave[i][ip])+' '
# b = ''
# for ip in range(nwavu):
#  b = b+np.str(lagsave[i][ip])+' '+np.str(siglagsave[i][ip])+' '
# 
# f.write(targ[i]+ ' '+np.str(cisqred[i])+' '+ a + b +'\n')
# 

# targw = targw + ' '+targ[i]
# lnow = ''
# lagw = ['']*nlc
# siglagw=['']*nlc
# for il in range(nlc):
#  lagw[il]= lagw[il]+' '+np.str(lagsave[i][il])
#  siglagw[i] = siglagw[i]+ ' '+np.str(siglagsave[i][il])
#
# parw = ['']*npar
# sigparw = ['']*npar
# for ip in range(npar):
#  parw[ip] = parw[ip] + 

#f.close()






#finish this later to be able to load pre calculated lags from earlier up to remake plots rather than recalculate everythong
#with open('contlags_op.dat') as f:
#    content = f.readlines()
#f.close()
#
#ntarg = len(content)
#targ=[]
#cisqred= [] 
#for i in range(ntarg):
# cnow = content[i].split()
# targ.append(cnow[0])
# cisqred.append(
#


#
##make plots
#os.chdir(pwd)
##load data (if just redooing plot can skip to here)
#
#parsave_plot = np.loadtxt('parsave_plot.txt')
#parsigsave_plot = np.loadtxt('parsigsave_plot.txt')
#lagsave_plot = np.loadtxt('lagsave_plot.txt')
#lagsigsave_plot = np.loadtxt('lagsigsave_plot.txt')
#cisqred_plot = np.loadtxt('cisqred_plot.txt')
#snrstatsave_plot = np.loadtxt('snrstatsave_plot.txt')
#
#
#
#
#
#
#nsave,npar = np.shape(parsigsave_plot)
#for ip in range(npar):
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# 
# ynow = []
# xnow = []
# znow = []
# 
# 
# alphanow = 1/parsigsave_plot[:,ip]
# alphanow = alphanow/np.max(alphanow)
# for it in range(nsave):
#  xnow.append(cisqred_plot[it])
#  ynow.append(parsave_plot[it,ip])
#  znow.append(alphanow[it])
#  ax1.scatter(xnow[it],ynow[it],znow[it],color='b')
# 
# 
# ax1.set_xlabel(r'$\chi^2/n$')
# ax1.set_ylabel(axlab[ip])
# 
# plt.savefig('fig_'+labsave[ip]+'.pdf')
#
#
#
#
#
#
##make plot of snrstat (max(flux) - min(flux) / median errorbar) vs i - g continuum lag
#idlaghi = 1
#idlaglo = 0
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#for it in range(nsave):
# x = lagsave_plot[it,idlaghi] - lagsave_plot[it,idlaglo]
# y = snrstatsave_plot[it]
# sigx = np.sqrt( lagsigsave_plot[it,idlaghi]**2 + lagsigsave_plot[it,idlaglo]**2 )
# ax1.errorbar(x,y,xerr=sigx,ls='',capsize = 10,color='k',alpha=0.3)
# ax1.scatter(lagsave_plot[it,idlaghi] - lagsave_plot[it,idlaglo],snrstatsave_plot[it],color='k')
#
#ylim = list(ax1.get_ylim())
#ax1.plot([0.,0.],ylim,color='k',linewidth=3)
#ax1.set_ylim(ylim)
#ax1.set_xlabel(r'$\langle \tau_i \rangle - \langle \tau_g \rangle$ (days)')
#ax1.set_ylabel('SNR = (Max flux - Min flux)/Median Error) \n averaged over all telscopes per filter')
#plt.savefig('fig_snrplot.pdf')
#
#   