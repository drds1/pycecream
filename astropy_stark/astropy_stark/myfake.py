#code to make fake light curves
#makes driver
#convolves with response function
#adds noise
#resamples appropriately

#INPUT wav[...nwav], snr[...nwav] lists of wavelengths and desired snr for the light curves
#..... dtmeanin[...nwav] list of the mean cadence for each light curve
#..... embh, degi, thi the balck hole mass, inclination and length of the light curve (starts at t = 0)

#OPTIONAL INPUT dl, z, er, emdot the luminosity distance in Mpc, redshift, eddington ratio and accretion rate emdot
#if emdot is -ve, use er (eddington ratio) to calculate mdot appropriate for this er.
#if emdot +ve ignore eddington ratio (er)

#iseed set to -1 to automatically generate driving light curve. Can manually set to a large positive
#...integer to get the same driving light curve each call


#OUTPUT tlcout[...nt], xlcout[...nt], echosave[...nt,...nwav],sigsave[...nt,...nwav] 
# 1 and 2-d lists of the timegrid of the driver and echolight curfves respectively
# no errors calculated for driving light curve include later if needed

import numpy as np
import pylab as plt
import os
import astropy_stark.mytfb_quick as tfb
import astropy_stark.myedlum as me
import astropy_stark.myconvolve as mc
import astropy_stark.mydisksim as mds
import astropy_stark.myfake_amp as mfa
import astropy_stark.myrandom as mr
import astropy_stark.mylcgen as mlg
import astropy_stark.myresample as mrs


def myfake(wavin, snr, dtmeanin, embh = 1.e7, degi = 0.0,
           thi = 100.0,gap = [], dlMpc = 70., z = 0.0,
           er=0.1, emdot=-1, eta = 0.1, dtres = 0.05,
           tauhi = 50., diag_plot = 0, iseed = -1,sampmin=0.8,
           meanforcein=[],sdforcein=[],sigforcein = [],
           mean_x = 0., sd_x = 1.,thfwhm = 0.2,
           thcent = 1.0,fakeplotnorm = 1,affine=0,
           paraffine=[-1,-2,-3,-4],
           dirsave = None,
           naffine=1000,noise_gap=0,
           tfx=0,T0v=1.e4,T0x=1.e4,sv=0.75,sx=0.75,
           verbose = False):
 '''
 generate fake light curve using a random walk driver and covolving with the accretion disc response fucntion
 described in Starkey et al 2016 'https://academic.oup.com/mnras/article-abstract/456/2/1960/1066664?redirectedFrom=PDF'
 :param wavin: list of wavelengths
 :param snr: list of signal to noise of each disk light curve as a relative to light curve standard deviation
 :param dtmeanin: average cadence for each response light curve
 OPTIONAL ARGUMENTS
 :param dtmeanin: mean cadence default = 1 day
 :param embh: black hole mass default = 1.e7
 :param degi: disk inlincation (default 0)
 :param thi: light curve time span (default = 100 days)
 :param gap:
 :param dlMpc:
 :param z:
 :param er:
 :param emdot:
 :param eta:
 :param dtres:
 :param tauhi:
 :param diag_plot:
 :param iseed: set to be a positive integer e.g 3423452 to always get the same results from the function call,
 else will get a new random light curve each time
 :param sampmin:
 :param meanforcein:
 :param sdforcein:
 :param sigforcein:
 :param mean_x:
 :param sd_x:
 :param thfwhm:
 :param thcent:
 :param fakeplotnorm:
 :param pycall:
 :param affine:
 :param paraffine:
 :param dirsave: if not None then save all output light curves to this directory in the
 directory structure required by the
 fortran version of cream. Will also general the required creamnames.dat file within the directory
 :param naffine:
 :param noise_gap:
 :param tfx:
 :param T0v:
 :param T0x:
 :param sv:
 :param sx:
 :return: output (dictionary of driver, echo light curve and response fn's)
 '''

 if dirsave is not None:
  os.system('rm -rf '+dirsave)
  os.system('mkdir '+dirsave)


 if (wavin == -1):
  wav = [-4000.0]
 else:
  wav = wavin
  
 nwav = len(wav)
 if (len(meanforcein) == 0): 
  meanforce = [-1]*nwav
 else:
  meanforce = list(meanforcein)
  
 if (len(sdforcein) == 0):
  sdforce = [-1]*nwav
 else:
  sdforce = list(sdforcein) 

 if (len(sigforcein) == 0):
  sigforce = [-1]*nwav
 else:
  sigforce = list(sigforcein)
  
 echosave = []
 sigsave  = []
 psi_out  = []
 if (emdot < 0):
  emdot = me.ermin_mdotout(embh,er,eta=eta)
 #print 'what is coming out here i dont understand',emdot
 
 #obtain disk spectrum
 fnu = []
 idw = 0
 for wavnow in wav:
  if (meanforce[idw] == -1): 
   fnunow = mds.mds([np.abs(wavnow)],embh,emdot,degi, dl = dlMpc, radlosch=3.0,radhisch=10000., ngrid = 1000, mjy = 1)
  else:
   fnunow = meanforce[idw]
  fnu.append(fnunow)
 fnu = np.array(fnu)  
 
 #generate fake light curve
 dat = mlg.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=thi+tauhi+2*dtres,dt=dtres,ploton=0,iseed=iseed,sdnorm = 1., meannorm = 0.)
 tlc = dat[:,0]
 xlc = dat[:,1]
 
 taunow = np.arange(0,tauhi,dtres)
 ntau   = np.shape(taunow)[0]
 idxlo  = np.where(tlc > tauhi)[0][0]#chop off lower limit where lookback is below the start time of the light curve
 
 tlclo  = tlc[idxlo]
 tlcout = tlc[idxlo:] #- tlclo
 xlcout = xlc[idxlo:]
 
 
 
 
 
 
 #define a list t_out to store the light curve times resampled to simulate uneven - random sampling
 t_out = []
 
 
 ntlc = np.shape(tlc[idxlo:])[0]
 for i in range(nwav):
  wavnow = wav[i]
  
  #what to do if we want the driver out
  if (wavnow == 0):
   echolc = 1.*xlc[idxlo:]
   sdwav = sd_x
   fdiskmean = mean_x
   psinow = np.zeros(ntau)
  else:
   
   #calculate expected mean
   if (meanforce[i] == -1):
    fdiskmean = fnu[i]
   else:
    fdiskmean = meanforce[i]
   #calculate expected variability amplitude
   if (sdforce[i] == -1):
    if fdiskmean == 0:
     fdfm = 1.0
    else:
     fdfm = fdiskmean
    sdwav = mfa.mfamp(embh,np.abs(wavnow),fdfm,thi, dMpc = dlMpc, z= z)
   else:
    sdwav = sdforce[i]
   
   #generate response function
   #print 'before response',i,wavnow,fdiskmean,emdot,embh,wavnow,'response function info',pycall
   psinow = tfb.pytfb_sub(taunow,embh,emdot,wavnow, degi, t0vin=-1, t0iin = -1, alpha_visc = -0.75, hxsch = 3.0, alpha_irad = -0.75, eta = 0.1, rlosch = 3.0, norm = 1,thcent = thcent, thfwhm = thfwhm)


   #convolve with response function for each wavelength
   echolc = mc.mc3(tlc,xlc,taunow,psinow)[idxlo:]
    
  psi_out.append(psinow)  
  fnunow = 1.*echolc
  
  
  
  #resample to simulate uneven observations 
  #if dtmeanin[i] <= 0 then dont resample just use same timegrid
  dtavenow = dtmeanin[i]
  if (dtavenow > 0):
   datin = np.zeros((ntlc,3))
   datin[:,0] = tlc[idxlo:]
   dout = mrs.myresample('',[''],dtavenow,sampmin=sampmin,sampcode=3,datin=datin)
   tout = dout[:,0]
  else:
   tout = tlc[idxlo:]  
  
  tplotlo = np.min(tout)
  
  ngap = len(gap)
  if (ngap != 0):
   ntout = np.shape(tout)[0]
   idexc = []
   idinc_gap = []
   for ig in range(ngap):
    tlonow = gap[ig][0]
    thinow = gap[ig][1]
    #idlo = np.where((tout - tplotlo > tlonow))[0]
    #idhi = np.where((tout - tplotlo < thinow))[0]
    #a = list(np.concatenate((idlo,idhi)))#list(np.where((tout < tlonow) & (tout > thinow))[0])
    a = np.where((tout - tplotlo > tlonow) & (tout - tplotlo < thinow))[0]
    #print a,tlonow, thinow, np.min(tout), np.max(tout),'bug hunt'
    idexc.append(a)
    idinc_gap.append(np.array([j4 for j4 in range(ntout) if j4 not in a]))
   idexc = np.array([j2 for i2 in idexc for j2 in i2]) 
   #idinc = np.arange(ntout)
   #idinc = np.where(idinc
   
   idinc = np.array([j3 for j3 in range(ntout) if j3 not in idexc])
   
   ninc = np.shape(idinc)[0]
   #raw_input()
   #tout_1stgap = tout[idinc_gap[0]]
   tout = tout[idinc]
   
    
  
  #fout_1stgap = np.interp(tout_1stgap,tlc[idxlo:],fnunow)
  fout = np.interp(tout,tlc[idxlo:],fnunow)
  
  
  
  
  #if have blr light curve, snr specifies variability amplitude / errobar
  #if (wavin[i] == -1.0):
  # fdiskmean = 0
  # sdwav = 1.0
   
  if (wavnow == 0.0):
   fdiskmean = mean_x
   sdwav = sd_x
  
  #scale to the correct flux and variability amplitude
  
  #if (onlyuse1stgap_4_scaling == 1):
  # emean = np.mean(fout_1stgap)
  # esd   = np.mean(fout_1stgap)
  # print 'only use 1st gap 4 scaling', emean,esd,sdwav,fdiskmean
  #else:
  emean = np.mean(fout)
  esd   = np.std(fout)
  
  fnunow = (fout - emean)*sdwav/esd + fdiskmean
    
  
  
  #plot during routine
  if (diag_plot == 1):
   fig = plt.figure()
   ax1 = fig.add_subplot(311)
   ax1.plot(tlcout,xlcout)
   ax2 = fig.add_subplot(312)
   ax2.plot(tout,fnunow)
   ax3 = fig.add_subplot(313)
   ax3.plot(taunow,psinow)
   if dirsave is not None:
    plt.savefig(dirsave+'/fig_diagplot_'+np.str(i)+'.pdf')
  
  
  
  #add noise (if have blr light curve, snr specifies variability amplitude / errobar)
  if (sigforce[i] == -1.0):
   #if (wavin[i] == -1 or wavnow == 0.0):
   sig    = np.zeros(np.shape(fnunow)[0])
   sig[:] = sdwav/snr[i]
   #else:
   # sig    = sdwav/snr[i]#fnunow/snr[i]#fnunow/snr[i]
  else:
   sig = np.ones(np.shape(fnunow)[0])*sigforce[i]
  
  
  #set snr based on sd at each interval rather than globally
  if (noise_gap == 1):
   tlonow = 0.0#tplotlo
   sig = np.ones(np.shape(fnunow)[0])  
   for ig in range(ngap):
    if (ig > 0):
     tlonow = gap[ig-1][1] + tplotlo
    #if (ig == 0):
    thinow = gap[ig][0] + tplotlo
    #elif (ig < ngap - 1):
    # thinow = 
    
    #if (ig < ngap - 1):
    # thinow = gap[ig][0] + tplotlo
    
    
    idincnow = np.where((tout > tlonow) & (tout < thinow))[0]
    stdnow   = np.std(fnunow[idincnow])
    meannow  = np.mean(fnunow[idincnow])
    sig[idincnow] = stdnow/snr[i]
    
    
    
    
  
  
  #print 'before adding noise',fnunow[:10],sigforce[i]
  nfnunow = np.shape(fnunow)[0]
  fnunonoise = 1.*fnunow
  fnunow = mr.normdis(nfnunow, fnunow, sig)
  #for i in range(np.shape(fnunow)[0]):
   #print fnunow[i],fnunonoise[i],np.abs(fnunow[i]-fnunonoise[i])/sig[i]
  #print('sd check n.o points outside error bars (should be ~ 0.32 for gaussian noise)', np.shape( np.where( np.abs(fnunow-fnunonoise)/sig > 1)[0] )[0]/(1.*np.shape(fnunow)[0]))
  sigsave.append(sig)
  t_out.append(tout)
  echosave.append(fnunow)#echosave.append(fnunow)
  
  
 
 
 
 
 #plot result to test
 #print 'making diag_plot', diag_plot
 if (diag_plot == 1):
  col = ['k','purple','b','cyan','green','orange','red']*100
  fig = plt.figure()
  ax1 = fig.add_subplot(211)
  ax2 = fig.add_subplot(212)
  for i in range(nwav):
   if (fakeplotnorm == 1):
    em = np.mean(echosave[i])
    esd = np.std(echosave[i])
    ax1.errorbar(t_out[i],(echosave[i] - em)/esd,sigsave[i]/esd,ls='',color=col[i])
   else:
    ax1.errorbar(t_out[i],echosave[i],sigsave[i],ls='',color=col[i])
   
   ax2.plot(taunow,psi_out[i],color=col[i])
   taumean = np.sum(taunow*psi_out[i])/np.sum(psi_out[i])
   ax2.plot([taumean,taumean],[0,1],color=col[i])
  ax2.set_xlim([0,30])
  ax2.tick_params(axis='x',which='minor',bottom='on')
   #print taunow[:10]
   #print psi_out[i][:10]
  if dirsave is not None:
   plt.savefig(dirsave+'/myfake_test.pdf')
 
 
 
 
 #if dirsave not None on then make a folder and save in the format required by cream
 if dirsave is not None:
  f = open('creamnames.dat','w')
  for i in range(nwav):
   wavnow = wav[i]
   tnow = np.array(t_out[i])
   
   
   fnow = np.array(echosave[i])
   signow = np.array(sigsave[i])
   nt = np.shape(tnow)[0]
   
 
   dat = np.zeros((nt,3))
   dat[:,0] = tnow
   dat[:,1] = fnow
   dat[:,2] = signow
   fname = 'synth_'+np.str(int(wavnow))+'.dat'
   #print fname
   os.getcwd()
   #print "'"+fname+"' "+np.str(wavnow)+"\n"     
   np.savetxt(fname,dat)
   f.write("'"+fname+"' "+np.str(wavnow)+"\n")
  f.close()
 
 
  if (tfx == 1):
   f = open('cream_tfx.par','w')
   f.write(np.str(T0v)+' 0.001 \n')
   f.write(np.str(T0x)+' 0.001 \n')
   f.write(np.str(sv)+' 0.001 \n')
   f.write(np.str(sx)+' 0.002 \n')
   f.close()
   os.system('mv cream_tfx.par ./'+dirsave)
   
  if (affine == 1):
   f = open('cream_affine.par','w')
   for ia in paraffine:
    f.write(np.str(ia)+' \n')
   f.write(np.str(naffine)+' \n')
   f.close()
   os.system('mv cream_affine.par ./'+dirsave)
  
  os.system('mv synth* ./'+dirsave)
  os.system('mv creamnames.dat ./'+dirsave)
  

 echo_lightcurves = []
 #print('nwav',nwav)
 #print('sigsave',sigsave)
 #print('end of sigsave')
 for i in range(nwav):
  t = t_out[i]
  x = echosave[i]
  s = sigsave[i]
  echo_lightcurves.append(np.array((t,x,s)).T)

 output = {'driver time axis':tlcout,
           'driver values':xlcout,
           'echo light curves':echo_lightcurves,
           'response function time axis':taunow,
           'response function values':psi_out,
           'accretion rate':emdot}
 return(output)






 




