#code to compare the model and lag estiamtes for cream fits with different
#error bar treatments for a few example targets
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import glob
import astropy_stark.cream_posterior as cpos





def lcplot(dnow,title='',idburnin=2./3,justth=0,justcont=0,plotinfo=1,
           plottrace=0,plots_per_page=5,xlclab = 'Time (HJD - 50,000)',xtflab ='lag (days)',forcelab=[],forcelag=[],sameplotdrive=1,extents=[],justnewsig=0,taumeanplot=1,tau90plot=0,postplot=1,header='',tauplot0=0,gplot=1,true=['','',np.log10(0.75)]):
 '''
 #input dirres (list of target directories from which to make plot (1plot per list element))
 #... if no ../output_20xxxx at end of each string, cream will look for these automatically and use the most
 #... recent simulation in te directory
 #tit... list of plot titles defaults to lcplot_x.pdf where x is the integer of list element
 plot the results of a cream simulation
 #if justth = 1 then only plot lines
 #if justcont = 1 then only plot continuum
 #if fitinfo = 1 then plot the parameter trace plots
 #if plottrace = 1 then plot the trace plots of the parameters
 #if plotsperpage == -1 then have one page per wavelength
 :param dnow:
 :param title:
 :param idburnin:
 :param justth:
 :param justcont:
 :param plotinfo:
 :param plottrace:
 :param plots_per_page:
 :param xlclab:
 :param xtflab:
 :param forcelab:
 :param forcelag:
 :param sameplotdrive:
 :param extents:
 :param justnewsig:
 :param taumeanplot:
 :param tau90plot:
 :param postplot:
 :param header:
 :param tauplot0:
 :param gplot:
 :param true:
 :return:
 '''
 taumean = 1
 idcount = 1
 dirres = []
 tit = []
 head = []
 if 'output_20' in dnow:
  dirres.append(dnow)
 else:
  a = sorted(glob.glob(dnow+'/output_20*'))[-1]
  dirres.append(a)
 if (title == ''):
  tit.append('lcplot_'+np.str(idcount))
 elif isinstance(title,str):
  tit.append(title+'_'+np.str(idcount))
 else:
  tit = list(title)
 head.append(header)

 '''
 load data from relevant files
 '''
 print('cream_lcplot plotting results from...',dnow)
 fileth = '/outputpars_th.dat'
 filetf = '/plots/modeltf.dat'
 filetfsig = '/plots/modeltf_sig.dat'
 filemod = '/plots/modellc.dat'
 filemodsig = '/plots/modellc_sig.dat'
 filevarexp = '/outputpars_varexpand.dat'
 filesigexp = '/outputpars2.dat'
 filedrive = '/plots/modeldrive.dat'
 filebof= '/testbofs.dat'
 filepspec = '/cream_furparms.dat'










 ndres = len(dirres)
 for idnow in range(ndres):


  #put this in a loop to make plots for all tested targets
  dnow = dirres[idnow]

  #read cream names file
  f = open(dnow+'/../creamnames.dat')
  filedat = f.readlines()
  filedat = [i.strip() for i in filedat]#remove \n
  wav     = [i.split(' ')[-1] for i in filedat]
  wav = [float(x) for x in wav]
  filedat = [i.split(' ')[0] for i in filedat]
  filedat = [i[1:-1] for i in filedat]
  f.close()
  nwav = len(wav)





  #!!!!!!!!!!!!!!!!!!!!!! load the model (all wavelengths in one file) !!!!!!!!!
  try:
   mod = np.loadtxt(dnow+'/'+filemod)
  except:
   mod = np.loadtxt(dnow+filemod)

  idxinc = np.where(mod[:,1] != 0)[0]
  tlo, thi = mod[idxinc[0],0],mod[idxinc[-1],0]
  mod = mod[idxinc,:]
  sigmod = np.loadtxt(dnow+'/'+filemodsig)
  tmod = mod[:,0]
  mod  = mod[:,1:]
  sigmod = sigmod[idxinc,1:]
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  #!!!!!!!!!!!!!!!!!!!!! load the response function !!!!!!!!!!!!!!!!!!!!!!!!!!!
  modtf = np.loadtxt(dnow+'/'+filetf)
  sigtf = np.loadtxt(dnow+'/'+filetfsig)
  tau   = modtf[:,0]
  modtf = modtf[:,1:]
  sigtf = sigtf[:,1:]

  #try to load top hat parameters, replace the disk response if non zero
  modth = np.loadtxt(dnow+'/'+fileth)[:,:nwav]
  nits = np.shape(modth[:,0])[0]
  modth = modth[int(idburnin*nits):,:]
  thstats = np.percentile(modth,[15.865,50,84.135],axis=0)
  thmed = thstats[1,:]
  idxth  = np.where(thmed != 0)[0]
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  #!!!!!!!!! load the error bar modification parameters calculate the median parameter and apply it
  try:
   sigexp = np.loadtxt(dnow+'/'+filesigexp)[:,2*nwav:]
   nits = np.shape(sigexp[:,0])[0]
   sigexp = sigexp[int(idburnin*nits):,:]
   semed = np.median(sigexp,axis=0)
  except:
   print('no error bar f factors')
   semed = np.ones(nwav)

  try:
   varexp = np.loadtxt(dnow+'/'+filevarexp)
   nits = np.shape(varexp[:,0])[0]
   varexp = varexp[int(idburnin*nits):,:]
   vamed = np.median(varexp,axis=0)
  except:
   print('no var expand parameters')
   vamed = np.zeros(nwav)

  try:
   #data required for outlier rejection
   fileoutrej = [dnow+'/plots/data_echo_dat_'+filedat[ilc] for ilc in range(len(wav))]#sorted(glob.glob(dnow+'/plots/data_echo_dat_*'))
   drej = [np.loadtxt(fileoutrej[idx]) for idx in range(nwav)]
  except:
   print('no outlier rejection found in...',dnow+'/plots/data_echo_dat_*')
   drej = [np.zeros((1,4)) for idx in range(nwav)]
















  #!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#
  #!!!!!!!!!!!!!!!# trace plots !!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!
  #!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#
  if (plottrace == 1):
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
    print('traceplot...',iynow,ixnow,np.mean(ptrace_plot[idx]), np.std(ptrace_plot[idx]))
   #plt.tight_layout()
   ax1.set_title(head[idnow])
   plt.savefig('traceplot_'+np.str(tit[idnow])+'.pdf')





  #!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#
  #!!!!!!!!!!!!!!!# driver, power spectrum and BOF plots !!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!
  #!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!#
  boffile = np.loadtxt(dnow + filebof)
  cisq = boffile[:,0]

  if (plotinfo == 1):
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

   plt.savefig('fitinfo_'+np.str(tit[idnow])+'.pdf')















  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





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
    ax1 = plt.subplot(gs1[0, 1:])
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1],color='k',linewidth=0.5)
    ax1.fill_between(datdrive[idxinc,0],datdrive[idxinc,1]-datdrive[idxinc,2],datdrive[idxinc,1]+datdrive[idxinc,2],alpha = 0.3,color='k')
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]-datdrive[idxinc,2],color='k')
    ax1.plot(datdrive[idxinc,0], datdrive[idxinc,1]+datdrive[idxinc,2],color='k')
    ax1.set_xticklabels([])
    ax1.tick_params(direction='in')
    ax1.set_yticks([])
    ax1.set_xlim([tlo,thi])


   for ilc in iwavinc_page[ipage]:

    ax1 = plt.subplot(gs1[ipnow, 1:])
    axtf = plt.subplot(gs1[ipnow,:1])

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
    drnow = drej[ilc]
    idrej = np.where(drnow[:,3] > 0)[0]
    ax1.scatter(drnow[idrej,0],drnow[idrej,1],color='r',marker='o')

    print(tmod[10],xm[10:15],ilc)
    ax1.plot(tmod,xm)
    ax1.fill_between(tmod,xm-sm,xm+sm,alpha=0.2)

    #set ylim by mod or data
    ymodlo = np.min(xdat[:]-sigdat[:])
    ymodhi = np.max(xdat[:]+sigdat[:])
    yrange = ymodhi - ymodlo
    ylim = [ymodlo-yrange/10,ymodhi+yrange/10]
    ax1.set_ylim(ylim)


    #check of we have top hat or disk response function
    if (ilc in idxth):
     axtf.hist(modth[:,ilc],bins=tau)
     thlags = np.percentile(modth[:,ilc],[15.865,50,84.135])
     lo = np.str(np.round(thlags[1]-thlags[0],2))
     med = np.str(np.round(thlags[1],2))
     hi = np.str(np.round(thlags[2] - thlags[1],2))
     lagtxt=r'$\tau_\mathrm{cent}='+med+'^{+'+hi+'}_{-'+lo+'}$'

     tfyl = list(axtf.get_ylim())
     axtf.set_ylim(tfyl)

     if (taumean == 1):
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
     axtf.plot(tau,modtf[:,itf])
     axtf.fill_between(tau,modtf[:,itf]-sigtf[:,itf],modtf[:,itf]+sigtf[:,itf],alpha=0.2)
     yltf = list(axtf.get_ylim())
     if (taumeanplot == 1):
      taumean = np.sum(tau*modtf[:,itf])/np.sum(modtf[:,itf])
      axtf.plot([taumean]*2,yltf,color='k')
     if (tau90plot > 0):
      Ntau = np.shape(tau)[0]
      psi90 = tau90plot
      psimax = np.max(modtf[:,itf])
      for itau in range(Ntau-1,0,-1):
       psinow = modtf[itau,itf]/psimax
       if (psinow > psi90):
        tau90 = tau[itau]
        id90 = itau
        break
      axtf.plot([tau[0],tau[-1]],[psi90]*2,color='r',ls='--')
      axtf.plot([tau90]*2,yltf,color='r',ls='--')
      if (tauplot0 == 1):
       axtf.plot([tau[0],tau[-1]],[0]*2,color='k')
    #miscellaneous label positions/tick mark info etc
    if (ipnow == nshapegs - 1):
     axtf.set_xlabel(xtflab)
     ax1.set_xlabel(xlclab)
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

   ax1.set_title(head[idnow])
   plt.savefig('page_'+np.str(np.int(ipage))+'_'+'lcplot_'+np.str(tit[idnow])+'.pdf')
   plt.close()


  if (postplot == 1):
   print('making posterior plot....','posterior_'+tit[idnow]+'.pdf')
   try:
    cpos.cream_posterior(dnow,true=true,header=head[idnow],extents_in=extents,fsave='posterior_'+tit[idnow]+'.pdf')
   except:
    print('unable to make covariance plot for disc posteriors. Please check at least some of these are set to vary'
          'in the fit.')




 #make the g plot
 if (gplot == 1):
  with open(dnow+'/cream_furparms.dat') as f:
   content = f.readlines()
   content = [x.strip() for x in content]
   f.close()
  x = [np.float(cnow) for cnow in content[0].split()]
  p0,dw,w0 = x
  dfur = np.loadtxt(dnow+'/cream_furparms.dat',skiprows=1)
  fnow,sfur,cfur,sd_sfur,sd_cfur = dfur.T
  sig0_2 = 0.5*p0*dw*((w0/2/np.pi)/fnow)**2
  sig_c2 = sd_cfur**2
  sig_s2 = sd_sfur**2
  Gsk = sig0_2 / ( sig0_2 + sig_s2 )
  Gck = sig0_2 / (sig0_2 + sig_c2)
  nfurtot = np.shape(Gck)[0]
  fig = plt.figure()
  ax2 = fig.add_subplot(111)
  ax2.plot(fnow/2/np.pi,Gsk,label='Sk',color='r')
  ax2.plot(fnow/2/np.pi,Gck,label='Ck',color='b')
  sum_gs = np.sum(Gsk)
  sum_gc = np.sum(Gck)
  sum_ave = (sum_gc + sum_gs)/2
  ax2.set_xlabel('frequency cycles/day')
  ax2.set_ylabel(r'$G=\sigma_0^2 / \left( \sigma_0^2 + \sigma^2 \right)$')
  ax2.set_title('G factor')
  ax2.text(0.98,0.55,r'$N_{\mathrm{eff}}=\sum_k \frac{G_{Sk} + G_{Ck}}{2} = '+np.str(np.int(sum_ave))+'$ of '+np.str(np.int(nfurtot)),ha='right',transform=ax2.transAxes,fontsize=14)
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  plt.legend()
  print(dnow+'/G_plot.pdf')
  plt.tight_layout()
  plt.savefig(dnow+'/G_plot.pdf')


  #save the effective chi squared and calculate the effective number of Fourier (and non Fourier parameters)
  cisqmed = np.median(cisq[int(idburnin*nits):])

  Nth = np.shape(idxth)[0]
  fdisk = np.loadtxt(dnow+'/outputpars.dat')[:,[2,3,4]]
  stddisk = np.std(fdisk,axis=0)
  Ndisk = np.shape(np.where(stddisk != 0)[0])[0]
  print('Nth ',Nth,' Ndisk',Ndisk)
  Np_notfur = 2*nwav + Ndisk + Nth
  if (np.mean(semed) == 0):
   semed_on = 0
  else:
   semed_on = 1
   Np_notfur = Np_notfur + nwav
  if (np.mean(vamed) == 0):
   vamed_on = 0
  else:
   vamed_on = 1
   Np_notfur = Np_notfur + nwav

  Np_fur = np.int(sum_ave)
  cisqmed_reduce = cisqmed/(Np_fur+Np_notfur)

  f = open(dnow+'/G_info.txt','w')
  f.write('N_freq_tot	N_freq_eff	N_non_freq	cisq(med)	cisq(med)/(N_non_freq + N_freq_eff)')
  op = np.str(nfurtot)+' '+np.str(Np_fur)+' '+np.str(np.int(Np_notfur))+' '+np.str(cisqmed_reduce)
  f.write(op+'\n')
  f.close()

 return()
