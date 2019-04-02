#code to compare the model and lag estiamtes for cream fits with different
#error bar treatments for a few example targets
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import glob
import astropy_stark.cream_posterior as cpos
idnow = 0

class plot_library:

 def __init__(self,directory):
  self.dnow = directory
  self.title = ''
  self.idburnin = 2. / 3
  self.justth = 0
  self.justcont = 0
  self.plotinfo = 1
  self.plottrace = 0
  self.plots_per_page = 5
  self.xlclab = 'Time (HJD - 50,000)'
  self.xtflab = 'lag (days)'
  self.forcelab = []
  self.forcelag = []
  self.sameplotdrive = 1
  self.extents = []
  self.justnewsig = 0
  self.taumeanplot = 1
  self.tau90plot = 0
  self.postplot = 1
  self.header = ''
  self.tauplot0 = 0
  self.gplot = 1
  self.true = ['', '', np.log10(0.75)]

  self.load_info()


 def load_info(self):
  idburnin = self.idburnin
  title = self.title
  header = self.header
  taumean = 1
  idcount = 1
  dirres = []
  tit = []
  head = []
  dnow = self.dnow
  if 'output_20' in dnow:
   dirres.append(dnow)
  else:
   a = sorted(glob.glob(dnow + '/output_20*'))[-1]
   dirres.append(a)
  if (title == ''):
   tit.append('lcplot_' + np.str(idcount))
  elif isinstance(title, str):
   tit.append(title + '_' + np.str(idcount))
  else:
   tit = list(title)
  head.append(header)
  self.head = head

  '''
  load data from relevant files
  '''
  print('cream_lcplot plotting results from...', dnow)
  fileth = '/outputpars_th.dat'
  filetf = '/plots/modeltf.dat'
  filetfsig = '/plots/modeltf_sig.dat'
  filemod = '/plots/modellc.dat'
  filemodsig = '/plots/modellc_sig.dat'
  filevarexp = '/outputpars_varexpand.dat'
  filesigexp = '/outputpars2.dat'
  filedrive = '/plots/modeldrive.dat'
  filebof = '/testbofs.dat'
  filepspec = '/cream_furparms.dat'

  ndres = len(dirres)
  for idnow in range(ndres):

   # put this in a loop to make plots for all tested targets
   dnow = dirres[idnow]

   # read cream names file
   f = open(dnow + '/../creamnames.dat')
   filedat = f.readlines()
   filedat = [i.strip() for i in filedat]  # remove \n
   wav = [i.split(' ')[-1] for i in filedat]
   wav = [float(x) for x in wav]
   filedat = [i.split(' ')[0] for i in filedat]
   filedat = [i[1:-1] for i in filedat]
   f.close()
   nwav = len(wav)

   # !!!!!!!!!!!!!!!!!!!!!! load the model (all wavelengths in one file) !!!!!!!!!
   try:
    mod = np.loadtxt(dnow + '/' + filemod)
   except:
    mod = np.loadtxt(dnow + filemod)

   idxinc = np.where(mod[:, 1] != 0)[0]
   tlo, thi = mod[idxinc[0], 0], mod[idxinc[-1], 0]
   mod = mod[idxinc, :]
   sigmod = np.loadtxt(dnow + '/' + filemodsig)
   tmod = mod[:, 0]
   mod = mod[:, 1:]
   sigmod = sigmod[idxinc, 1:]
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   # !!!!!!!!!!!!!!!!!!!!! load the response function !!!!!!!!!!!!!!!!!!!!!!!!!!!
   modtf = np.loadtxt(dnow + '/' + filetf)
   sigtf = np.loadtxt(dnow + '/' + filetfsig)
   tau = modtf[:, 0]
   modtf = modtf[:, 1:]
   sigtf = sigtf[:, 1:]

   # try to load top hat parameters, replace the disk response if non zero
   modth = np.loadtxt(dnow + '/' + fileth)[:, :nwav]
   nits = np.shape(modth[:, 0])[0]
   modth = modth[int(idburnin * nits):, :]
   thstats = np.percentile(modth, [15.865, 50, 84.135], axis=0)
   thmed = thstats[1, :]
   idxth = np.where(thmed != 0)[0]
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   # !!!!!!!!! load the error bar modification parameters calculate the median parameter and apply it
   try:
    sigexp = np.loadtxt(dnow + '/' + filesigexp)[:, 2 * nwav:]
    nits = np.shape(sigexp[:, 0])[0]
    sigexp = sigexp[int(idburnin * nits):, :]
    semed = np.median(sigexp, axis=0)
   except:
    print('no error bar f factors')
    semed = np.ones(nwav)

   try:
    varexp = np.loadtxt(dnow + '/' + filevarexp)
    nits = np.shape(varexp[:, 0])[0]
    varexp = varexp[int(idburnin * nits):, :]
    vamed = np.median(varexp, axis=0)
   except:
    print('no var expand parameters')
    vamed = np.zeros(nwav)

   try:
    # data required for outlier rejection
    fileoutrej = [dnow + '/plots/data_echo_dat_' + filedat[ilc] for ilc in
                  range(len(wav))]  # sorted(glob.glob(dnow+'/plots/data_echo_dat_*'))
    drej = [np.loadtxt(fileoutrej[idx]) for idx in range(nwav)]
   except:
    print('no outlier rejection found in...', dnow + '/plots/data_echo_dat_*')
    drej = [np.zeros((1, 4)) for idx in range(nwav)]
  self.idxinc = idxinc
  self.semed = semed
  self.vamed = vamed
  self.tau = tau
  self.modtf = modtf
  self.drej = drej
  self.idxth = idxth
  self.sigtf = sigtf
  self.tmod = tmod
  self.mod = mod
  self.tlo = tlo
  self.thi = thi
  self.sigmod = sigmod
  self.taumean = taumean
  self.modth = modth
  self.nwav = nwav
  self.tit = tit




 def plot_trace(self):
  plt.close()
  head = self.head
  nwav = self.nwav
  '''
  :return:
  '''
  dnow = self.dnow
  filetrace = ['/outputpars.dat','/outputpars2.dat','/outputpars_th.dat','/outputpars_varexpand.dat']
  idtrace = [[2,3,4],list(np.arange(3*nwav)),list(np.arange(2*nwav)),list(np.arange(nwav))]
  tit_trace = [[r'$\log \dot{M}$',r'$\cos i$',r'$\alpha$'],
  ['stretch '+np.str(idx) for idx in range(nwav)] + ['offset '+np.str(idx) for idx in range(nwav)] + ['sig f '+np.str(idx) for idx in range(nwav)],
  [r'$\tau_\mathrm{cent}$ '+np.str(idx) for idx in range(nwav)] +  [r'$\tau_\mathrm{width}$ '+np.str(idx) for idx in range(nwav)],
  ['sig_V '+np.str(idx) for idx in range(nwav)]
  ]

  nft = len(filetrace)
  idtot  = 0
  ptraceop = []
  ttraceop = []
  for idx in range(nft):
   npar = len(idtrace[idx])
   try:
    dtrace = np.loadtxt(dnow+filetrace[idx])[:,idtrace[idx]]
    vartrace = np.std(dtrace,axis=0)
    idinc = np.where(vartrace > 0)[0]
    d_trace = [dtrace[:,idx2] for idx2 in idinc]
    t_trace = [tit_trace[idx][id2] for id2 in range(npar) if id2 in idinc]
    ptraceop.append(d_trace)


    ttraceop.append(t_trace)
   except:
    print('cannot find file for trace plots in...')
    print(dnow+filetrace[idx])

  ptrace_plot = [j for i in ptraceop for j in i]
  ttrace_plot = [j for i in ttraceop for j in i]

  ntrace = len(ptrace_plot)
  nalong = min(4,ntrace)


  ndown = np.int(np.ceil(1.*ntrace/nalong))
  gs1 = gridspec.GridSpec(ndown, nalong)
  gs1.update(left=0.1, right=0.9, wspace=0.05,hspace = 0.0,bottom=0.1,top=0.99)

  iynow = -1
  for idx in range(ntrace):
   ixnow = np.mod(idx,nalong)
   if (ixnow == 0):
    iynow = iynow + 1
   ax1 = plt.subplot(gs1[iynow, ixnow])
   ax1.plot(ptrace_plot[idx])
   ax1.set_ylabel(ttrace_plot[idx])
  ax1.set_title(head[idnow])

  return(ax1)






 def plot_driver(self):
  plt.close()
  filedrive = '/plots/modeldrive.dat'
  filebof = '/testbofs.dat'
  filepspec = '/cream_furparms.dat'
  dnow = self.dnow
  boffile = np.loadtxt(dnow + filebof)
  cisq = boffile[:,0]
  idxinc = self.idxinc
  head = self.head

  if (self.plotinfo == 1):
   boffur = boffile[:,1]
   boftot = cisq + boffur

   pspecfile = np.loadtxt(dnow + filepspec,skiprows=1)
   freq = pspecfile[1:,0]/2/np.pi
   pspec = pspecfile[1:,[1,2]]

   datdrive = np.loadtxt(dnow + filedrive)



   ndown = 2
   nalong = 2
   gs1 = gridspec.GridSpec(ndown, nalong)
   gs1.update(left=0.1, right=0.9, wspace=0.3,hspace = 0.3,bottom=0.1,top=0.99)

   ax1 = plt.subplot(gs1[0, :])
   ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1],color='k',linewidth=0.5)
   ax1.fill_between(datdrive[idxinc,0],datdrive[idxinc,1]-datdrive[idxinc,2],datdrive[idxinc,1]+datdrive[idxinc,2],alpha = 0.3,color='k')
   ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]-datdrive[idxinc,2],color='k')
   ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]+datdrive[idxinc,2],color='k')
   ax1.set_xlabel('time (days)')
   ax1.set_ylabel('X(t) (arbitrary units)')

   ax1 = plt.subplot(gs1[1,0])
   ax1.plot(boftot,label=r'$\chi^2+\mathrm{BOF}_\mathrm{fur}$' )
   ax1.plot(cisq,label=r'$\chi^2$')
   ax1.plot(boffur,label=r'$\mathrm{BOF}_\mathrm{fur}$')
   ax1.set_xlabel('Iteration')
   ymin = np.min(boffur)
   ymax = np.mean(boftot) + np.std(boftot)*2#ymin + np.max(
   ylim = [ymin,ymax]
   ax1.set_ylim(ylim)
   ax1.set_ylabel('BOF')
   plt.legend(fontsize='small')
   ax1 = plt.subplot(gs1[1,1])
   ps = pspec[:,0]**2 + pspec[:,1]**2
   ax1.scatter(freq,ps)
   try:
    sigpspec = pspecfile[1:,[3,4]]
    ps_sd = 2* np.sqrt( (pspec[:,0]*sigpspec[:,0])**2 + (sigpspec[:,1]*sigpspec[:,1])**2 )
    ax1.errorbar(freq,ps,ps_sd)
    #ax1.fill_between(freq,ps,ps-ps_sd,ps+ps_sd)
   except:
    print('cannot find uncertanties for power spectrum plot')

   xlim = [freq[0],freq[-1]]
   ax1.set_xlim(xlim)
   ax1.set_xscale('log')
   ax1.set_yscale('log')
   ax1.set_xlabel('frequency (cycles/day)')
   ax1.set_ylabel('P(f)')
   ax1.set_title(head[idnow])

   return(ax1)
   #plt.savefig('fitinfo_'+np.str(tit[idnow])+'.pdf')





 def plot_lightcurves(self):

  plt.close()
  axtf = None
  dnow = self.dnow
  justcont = self.justcont
  justth = self.justth
  plots_per_page = self.plots_per_page
  sameplotdrive = self.sameplotdrive
  semed = self.semed
  vamed = self.vamed
  justnewsig = self.justnewsig
  taumeanplot = self.taumeanplot
  tau90plot = self.tau90plot
  tauplot0 = self.tauplot0
  forcelab = self.forcelab
  forcelag = self.forcelag
  idxinc = self.idxinc
  idxth = self.idxth
  tlo, thi = self.tlo, self.thi
  mod = self.mod
  sigmod = self.sigmod
  tit = self.tit

  fig_objects = []



  # read cream names file
  f = open(dnow + '/../creamnames.dat')
  filedat = f.readlines()
  filedat = [i.strip() for i in filedat]  # remove \n
  wav = [i.split(' ')[-1] for i in filedat]
  wav = [float(x) for x in wav]
  filedat = [i.split(' ')[0] for i in filedat]
  filedat = [i[1:-1] for i in filedat]
  f.close()
  nwav = len(wav)


  nwavplot = nwav
  if (justcont == 1):
   idth = np.where(np.array(wav) == -1)[0]
   nth = np.shape(idth)[0]
   nwavplot = nwav - nth

  if (justth == 1):
   idcont = np.where(np.array(wav) != -1)[0]
   ncont = np.shape(idcont)[0]
   nwavplot = nwav - ncont

  ####!!!!       initialise the figure
  #grid spec gives better customisation of plot layout. Copied and pasted from
  #pyplot_cream_easy_2.py in mcmcmulti3 directory other codes below customise for this case






  iwavinc = [i for i in range(nwav)]
  if (justth == 1):
   iwavinc = [i for i in range(nwav) if wav[i] == -1]
  if (justcont == 1):
   iwavinc = [i for i in range(nwav) if wav[i] != -1]
  nwavinc = len(iwavinc)

  if (plots_per_page == -1):
   wavu = np.sort(np.unique(wav))
   npages = len(wavu)
  else:
   npages   = np.int(np.ceil(1.0*nwavinc/plots_per_page))


  #organise how to distribute the plots
  if (plots_per_page == -1):
   iwavinc_page = []
   for iw2 in range(npages):
    wavunow = wavu[iw2]
    iwnow = [iwavinc[i2] for i2 in range(nwavinc) if wav[iwavinc[i2]] == wavunow ]
    iwavinc_page.append(iwnow)

  elif (nwavinc > plots_per_page):
   iwavinc_page = [iwavinc[i:i+npages] for i  in range(0, len(iwavinc), npages)]

  else:
   iwavinc_page = [iwavinc]



  for ipage in range(npages):
   nwavplot = len(iwavinc_page[ipage])
   fig = plt.figure()
   if (sameplotdrive == 1 and ipage == 0):
    ipnow = 1
    gs1 = gridspec.GridSpec(nwavplot+1, 4)
    nshapegs = nwavplot+1
   else:
    ipnow = 0
    gs1 = gridspec.GridSpec(nwavplot, 4)
    nshapegs = nwavplot


   gs1.update(left=0.1, right=0.9, wspace=0.05,hspace = 0.0,bottom=0.15,top=0.99)
   ###!!!!!


   #if ipage == 0 then plot the driver up top
   if (sameplotdrive == 1 and ipage == 0):
    filedrive = '/plots/modeldrive.dat'
    datdrive = np.loadtxt(dnow + filedrive)
    ax1 = fig.add_subplot(gs1[0, 1:])
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1],color='k',linewidth=0.5)
    ax1.fill_between(datdrive[idxinc,0],datdrive[idxinc,1]-datdrive[idxinc,2],datdrive[idxinc,1]+datdrive[idxinc,2],alpha = 0.3,color='k')
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]-datdrive[idxinc,2],color='k')
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]+datdrive[idxinc,2],color='k')
    ax1.set_xticklabels([])
    ax1.tick_params(direction='in')
    ax1.set_yticks([])
    ax1.set_xlim([tlo,thi])


   for ilc in iwavinc_page[ipage]:
    ax1 = fig.add_subplot(gs1[ipnow, 1:])
    axtf = fig.add_subplot(gs1[ipnow,:1])
    #ax1 = plt.subplot(gs1[ipnow, 1:])
    #axtf = plt.subplot(gs1[ipnow,:1])

    dat = np.loadtxt(dnow+'/../'+filedat[ilc])
    tdat = dat[:,0]
    xdat = dat[:,1]
    sigdat = dat[:,2]
    sigdat_new = np.sqrt((semed[ilc]*sigdat)**2 + vamed[ilc])
    xm = mod[:,ilc]
    sm = sigmod[:,ilc]

    #decide order of old/new error bar plotting
    for it in range(np.shape(xdat)[0]):
     signew = sigdat_new[it]
     sigold = sigdat[it]
     if (signew > sigold):
      ax1.errorbar([tdat[it]],[xdat[it]],[signew],ls='',color='b')
      if (justnewsig == 0):
       ax1.errorbar([tdat[it]],[xdat[it]],[sigold],ls='',color='r')
     else:
      if (justnewsig == 0):
       ax1.errorbar([tdat[it]],[xdat[it]],[sigold],ls='',color='r')
      ax1.errorbar([tdat[it]],[xdat[it]],[signew],ls='',color='b')



    #plot outlier rejected points if any
    drnow = self.drej[ilc]
    idrej = np.where(drnow[:,3] > 0)[0]
    ax1.scatter(drnow[idrej,0],drnow[idrej,1],color='r',marker='o')

    print(self.tmod[10],xm[10:15],ilc)
    ax1.plot(self.tmod,xm)
    ax1.fill_between(self.tmod,xm-sm,xm+sm,alpha=0.2)

    #set ylim by mod or data
    ymodlo = np.min(xdat[:]-sigdat[:])
    ymodhi = np.max(xdat[:]+sigdat[:])
    yrange = ymodhi - ymodlo
    ylim = [ymodlo-yrange/10,ymodhi+yrange/10]
    ax1.set_ylim(ylim)


    #check of we have top hat or disk response function
    if (ilc in idxth):
     axtf.hist(self.modth[:,ilc],bins=self.tau)
     thlags = np.percentile(self.modth[:,ilc],[15.865,50,84.135])
     lo = np.str(np.round(thlags[1]-thlags[0],2))
     med = np.str(np.round(thlags[1],2))
     hi = np.str(np.round(thlags[2] - thlags[1],2))
     lagtxt=r'$\tau_\mathrm{cent}='+med+'^{+'+hi+'}_{-'+lo+'}$'

     tfyl = list(axtf.get_ylim())
     axtf.set_ylim(tfyl)

     if (self.taumean == 1):
      axtf.text(0.9,0.8,lagtxt,ha='right',transform=axtf.transAxes,fontsize=8)
      axtf.plot([thlags[0]]*2,tfyl,ls=':',label=None,color='k')
      axtf.plot([thlags[1]]*2,tfyl,ls='-',label=None,color='k')
      axtf.plot([thlags[2]]*2,tfyl,ls=':',label=None,color='k')
    else:
     #disk response function
     #response function only saved once in wavelength file
     #if have multiple lc's at same wavelength, need to plot the response function
     #from the first light curve at that wavelength
     wavnow = wav[ilc]
     idsamewav = np.where(np.array(wav) == wavnow)[0]
     if (np.shape(idsamewav)[0] > 0):
      itf = idsamewav[0]
     else:
      itf = ilc
     axtf.plot(self.tau,self.modtf[:,itf])
     axtf.fill_between(self.tau,self.modtf[:,itf]-self.sigtf[:,itf],self.modtf[:,itf]+self.sigtf[:,itf],alpha=0.2)
     yltf = list(axtf.get_ylim())
     if (taumeanplot == 1):
      taumean = np.sum(self.tau*self.modtf[:,itf])/np.sum(self.modtf[:,itf])
      axtf.plot([taumean]*2,yltf,color='k')
     if (tau90plot > 0):
      Ntau = np.shape(self.tau)[0]
      psi90 = tau90plot
      psimax = np.max(self.modtf[:,itf])
      for itau in range(Ntau-1,0,-1):
       psinow = self.modtf[itau,itf]/psimax
       if (psinow > psi90):
        tau90 = self.tau[itau]
        id90 = itau
        break
      axtf.plot([self.tau[0],self.tau[-1]],[psi90]*2,color='r',ls='--')
      #axtf.plot([self.tau90]*2,yltf,color='r',ls='--')
      if (tauplot0 == 1):
       axtf.plot([self.tau[0],self.tau[-1]],[0]*2,color='k')
    #miscellaneous label positions/tick mark info etc
    if (ipnow == nshapegs - 1):
     axtf.set_xlabel(self.xtflab)
     ax1.set_xlabel(self.xlclab)
    else:
     ax1.set_xticklabels([])
     axtf.set_xticklabels([])
    if (ipnow == nshapegs/2):
     axtf.set_ylabel('response (relative units)')
     ax1.set_ylabel('flux (relative units)')





    ax1.tick_params(direction='in')
    axtf.set_yticklabels([])
    axtf.set_yticks([])
    ax1.set_yticks([])
    ax1.yaxis.set_label_position('right')
    ax1.set_yticklabels([])

    if (forcelag != []):
     axtf.set_xlim(forcelag)

    if (forcelab == []):
     ax1.text(0.9,0.8,filedat[ilc],ha='right',transform=ax1.transAxes,fontsize=8)
    elif (forcelab == 0):
     pass
    else:
     ax1.text(0.9,0.8,forcelab[ilc],ha='right',transform=ax1.transAxes,fontsize=8)

    #title
    if (ipnow == 0):
     ax1.text(0.5,1.2,tit[idnow],ha='center',transform=ax1.transAxes,fontsize=14)

    ax1.set_xlim([tlo,thi])
    ipnow = ipnow + 1

   ax1.set_title(self.head[idnow])
   fig_objects.append(fig)
  return(fig_objects)
   #plt.savefig('page_'+np.str(np.int(ipage))+'_'+'lcplot_'+np.str(tit[idnow])+'.pdf')



 def plot_posterior(self):
  '''
  plot the disk posterior covariances
  :return:
  '''
  dnow = self.dnow
  head = self.header
  extents = self.extents
  tit = self.title
  true = self.true
  #fsave = 'posterior_'+tit[idnow]+'.pdf'
  fsave = ''

  try:
   fig, axes = cpos.cream_posterior(dnow,true=true,header=head[idnow],extents_in=extents,fsave=fsave)
  except:
   print('unable to make covariance plot for disc posteriors. Please check at least some of these are set to vary'
         'in the fit.')
   fig = None
   axes = None

  return(fig,axes)
