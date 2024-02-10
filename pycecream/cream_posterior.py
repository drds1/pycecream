import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
from astropy_stark.mycorner_jul6_16 import *
from astropy_stark.myedlum import *
from astropy_stark.mytemp0 import *
import matplotlib

#font = {'family' : 'normal',
#        'size'   : 20}
#
#matplotlib.rc('font', **font)


def cream_posterior(dirin,ibin = 2./3,
                    npoints = -1,
                    idT1=2,filepar='outputpars.dat',
                    extents_in=[],
                    header='',
                    idxres=[2,3,4],
                    special=[1,0,0],
                    true=['','',np.log10(0.75)],
                    fsave=''):
 #dat = np.loadtxt('/Users/ds207/Documents/standrews/sta/fort/fortcode/mcmcmultiresults/16feb_4151_alls/output_20170228_003/outputpars.dat')
 truei = list(true)
 #stack posterior probability distributions together
 icheck = isinstance(dirin,(list)) 
 ic = 0
 if (icheck == 1):
  for dir in dirin:
   try:
    d = np.loadtxt(dir+'/'+filepar,ndmin=2)
    print(dir+'/'+filepar,np.shape(d))
    ndown,nalong = np.shape(d)
    idlo = np.int(ibin*ndown) 
    if (ic == 0):
     dat = d[idlo:,:]
    else:
     dat = np.vstack((dat,d))
    ic = ic + 1 
   except:
    print('problem with output parameters in', dir+'/'+filepar)
  dir = dirin[0]
 else:
  dat = np.loadtxt(dirin+'/'+filepar,ndmin=2) 
  dir = dirin

 
 
 if (filepar == 'outputpars.dat'):
  datnorm = 1
 #first off do the accretion disc parameters then try the top hat centroids 
 #try:
 # 
 # datnorm = 1
 # idxres  = [2,3,4]
 # special = [1,0,0]#is special = 2 take arc cos of data before potting
  
 #except:
 # print 'no outputpars.dat file must be in tfx mode',dir+'/outputpars.dat'
 # dat = np.loadtxt(dir+'/outputpars_tfxdisk.dat')
 # datnorm = 0
 # idxres = [0,1,2,3,4]
 # special = [0,0,1,1,0]
 
 
 try:
  ed   = np.loadtxt(dir+'/plots/cream_miscpar.dat')
  embh = ed[0]
  eta  = ed[1]
 except:
  embh = 1.e7
  eta  = 0.1
 
 
 
 
 
 
 #find if mmdot is varied
 if (datnorm == 1):
  try:
   imdot = np.where(np.array(idxres) == 2)[0][0]
   if (special[imdot] == 1):
    logmmdot = 1
   else:
    logmmdot = 0
  except:
   logmmdot = 0
   imdot = -1
  
  
  
  
  
  #if T1 = 1 express temperatures as T1 here
  if (idT1 > 0):
   idslope = 4
   emlog = np.log10(embh)
   dotmmlog = emlog+np.log10(dat[:,idT1])
   alphain  = dat[:,idslope]
   temp = temp0(dotmmlog, emlog, sigdotmmlog =0, sigemlog = 0.0, alpha_in = 0.75, sig_alphain=0.0,eta=0.1, alb=0, hxrs = 3,r0ld = 1)
   temp_info = temp[0]
   if (truei[0] != ''):
    truei[0] = np.log10(temp0(emlog+np.log10(truei[0]), emlog, sigdotmmlog =0, sigemlog = 0.0, alpha_in = 0.75, sig_alphain=0.0,eta=0.1, alb=0, hxrs = 3,r0ld = 1)[0])
 
 
 
 nver = np.shape(dat)[0]
 
 
 idxlo = int(np.floor(ibin*nver))
 dat = dat[:,idxres]
 nver  = np.shape(dat[:,0])[0]
 nhor = np.shape(dat)[1]
 
 
 #calculate edingon ratios
 #embh = 3.e7
 #eta = 0.1
 
 
 
 
 
 

 
 
 ivar = [1]*nhor
 for i in range(nhor):
  idspec = special[i]
  if (idspec == 2):
   dat[:,i] = np.arccos(dat[:,i])*180/np.pi
  if (idspec == 1):
   dat[:,i] = np.log10(dat[:,i])
  
   
   
  sd = np.std(dat[:50,i])
  if (sd < 1.e-10):
   ivar[i]=0
 
 ivar = np.array(ivar)
 idxvar = list(np.where(ivar > 0)[0])
 #idxvar = idxvar[:-1]

 
 

 dat = dat[:,idxvar]

 
 #function to cast number in standard form
 def sci_notation(number, sig_fig=2):
   ret_string = "{0:.{1:d}e}".format(number, sig_fig)
   a,b = ret_string.split("e")
   b = int(b) #removed leading "+" and strips leading zeros too.
   return "$"+a + " \\times 10^" + str(b) +"$"
 
 
 
 if (((embh > 0) or (imdot != -1)) and (datnorm == 1)):
  if (logmmdot == 1):
   dotm = 10**dat[:,0]
  else:
   dotm = dat[:,0]
  edd_rat = edd(embh,dotm,eta)[2]
  idxsort     = np.argsort(edd_rat)
  ns          = np.shape(idxsort)[0]
  ers         = edd_rat[idxsort]
  idlo        = int(max(0,np.floor(0.17*ns)))
  idhi        = int(min(np.ceil((1.-0.17)*ns),ns-1))
  ermed       = ers[ns/2]
  erlo        = ers[idlo]
  erhi        = ers[idhi]
  siglo       = ermed - erlo
  sighi       = erhi  - ermed
  aner        = ', $L/L_{edd}='+np.str(np.round(ermed,2))+'_{-'+np.str(np.round(siglo,2))+'}^{+'+np.str(np.round(sighi,2))+'}$ for $M_\\mathrm{BH}=$'+sci_notation(embh)+'$M_\odot$'
 else:
  aner = ''
 
 if (npoints == -1):
  iskip = 1
 else:
  iskip = nver/npoints
 
 
 
 #truths = [p0,w0,alpha,beta]
 #annotate = [r'$\dot{M}/M_\odot$yr$^-1$',r'$\cos \left( i \right)$',r'$\alpha$']
 
 
 
 costick = np.cos(np.pi/180*np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90])) 
 #coslab = [np.str(np.arccos(costick[i])*180/np.pi) for i in range(len(costick))]
 coslab = ['0','','','','','','','','40','','','','60','','','','80','','']
 #coslabblank = ['' for i in range(len(costick))]
 
 if (datnorm == 1):
  if (extents_in == []):
   extents=[0.9999,(0.,1.),0.9999]
  else:
   extents=extents_in
  
  annotate = [r'$\log ( T_1 $/K$  )$',r'$i$',r'$\log \alpha$']
  ann      = ['','','']#[aner,'','']
  smpi     = [r'$\log ( T_1 $/K$  )$',r'cos $i$',r'$\log \alpha$']
  specax = [[],costick,[]]
  specaxlab = [[],coslab,[]]
  if (idT1 > 0):
   annotate[0] = r'$\log T_1$ (K)'
   smpi[0] = r'$\log T_1$ (K)'
 
 else:
  if (extents_in == []):
   extents=[0.9999,(0.,1.),0.9999]
  else:
   extents=extents_in

  annotate = [r'$\log \alpha_v$',r'$\log \alpha_i$',r'$T_{1v}$',r'$T_{1i}$',r'$i$']
  ann      = ['']*5
  smpi     = [r'$\log \alpha_v$',r'$\log \alpha_i$',r'$T_{1v}$',r'$T_{1i}$',r'$i$']
  specax = [[],[],[],[],costick]
  specaxlab = [[],[],[],[],coslab]
  
 
 #print 'idx var...',idxvar
 a = [annotate[i] for i in idxvar]
 b = [smpi[i] for i in idxvar]
 c = [extents[i] for i in idxvar]
 d = [ann[i] for i in idxvar]
 specax = [specax[i] for i in idxvar]
 specaxlab = [specaxlab[i] for i in idxvar]
 truth_in = [truei[i] for i in idxvar]
 
 #extents=[(0.,1.5),(0.,90.)]
 ann_in   = [annotate[i] for i in idxvar]
 extra_an = [' ' for i in idxvar]#[ann[i] for i in idxvar]
 #a = mc.corner_1(dat,plot_contours=0,sigconts=[100.-68,100.-95,100-99.7],skip = iskip,plot_datapoints=True,
 #sigmedann_pre_in=b,labels=a,extents=[0.99,(0.,90.)],figname = 'parmparmplot.pdf')
 
 
 
 datcorner = np.array(dat)
 if ((idT1> 0) and (datnorm == 1)):
  datcorner[idxlo:,0] = np.log10(temp_info[idxlo:])
 
 #tr slope should be plotted in log
 for idnow in range(len(b)):
  if ('alpha' in b[idnow]):
   datcorner[idxlo:,idnow] = np.log10(dat[idxlo:,idnow])
  
  
 if (fsave==''):
  figname = 'posterior.pdf'
 else:
  figname = fsave
 
 print('making posterio plots...cream_posterior', truth_in)
 a = corner_1(datcorner[idxlo:,:],plot_contours=0,sigconts=[100.-68],skip = iskip,plot_datapoints=True,sigmedann_pre_in=b,sigmedann_post_in = d, 
 labels=a,extents=c,annotate = extra_an, sigmedann=1, 
 xtxtcoord = 1.25,ytxtcoord=1.05,lbm_in = 0.65,ltr_in = 0.4,xlab_coord=-0.3,ylab_coord=-0.3,
 specax=specax,specaxlab=specaxlab,truths=truth_in,title=header, 
 truth_color='r', figname = fsave)
 #except:
 # print 'cannot make continuum mdot vs inc posterior plot (maybe no continuum light curves)'
 # pass
 
 
 #try to make line plots of all the varied parameters
 try:
  nvar = len(idxvar)
  fig = plt.figure()
  for i in range(nvar):
   ax1 = fig.add_subplot(nvar,1,i+1)
   ax1.plot(dat[:,i])
   ax1.set_xlabel('Iteration')
   ax1.set_ylabel(annotate[i])
  plt.savefig('fig_cream_lines.pdf')
 except:
  pass
 
 #now do for the top hat parameters
 try:
  npoints = 200
  dat = np.loadtxt('../outputpars_th.dat') 
  nlc = np.shape(dat[0,:])[0]/2
  anthcen = ['lc (cent) '+str(i) for i in range(nlc)]
  anthwid = ['lc (fwhm) '+str(i) for i in range(nlc)]
  
  ann = anthcen + anthwid
  
  
  
  #check which points to plot
  npar = 2*nlc
  idxin = [i for i in range(npar) if np.shape(np.unique(dat[:,i]))[0] > 1]
  dat = dat[:,idxin]
  ann = [ann[i] for i in range(npar) if i in idxin ]
  
  nver = np.shape(dat[:,0])[0]
  iskip = max(1,nver/npoints)
   
  a = corner_1(dat[idxlo:,:],plot_contours=0,sigconts=[100.-68],skip = iskip,plot_datapoints=True,
  labels=ann, title=title, xtxtcoord=1.45,ytxtcoord=0.95,lbm_in = 0.65,ltr_in = 0.4,
  xlab_coord=-1.9,ylab_coord=-1.9,figname = 'posterior_th.pdf')
 except:
  print('cannot make top hat posterior plot (maybe no BLR top hat light curves)')
  pass
  
  
 return()
