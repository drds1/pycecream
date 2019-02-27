import numpy as np
import vortexa_ccf_simple as mccf
import matplotlib.pylab as plt
import os
import matplotlib.gridspec as gridspec
import myconvolve as  mc3
import mylcgen as mlc
import myresample as mrs
import vortexa_mod_randomwalk as mf18
import pandas as pd
import datetime

#enter the time series signals as list of lists. Driver first [[time,values,noise],[time2,values2,noise2]...etc]
#titles for plots
def lag_rwmod(signals,titles=[],
 monte_carlo_iterations=10000,
 extra_frequencies = [],
 random_walk_gridres = 1.0,
 summary_plot='',
 nsim = 500,
 refdate=0,
 laglim=[]):
 
 
 time_model = []
 signal_model = []
 error_model  = []
 model_hi = []
 model_lo = []
 x_iterations = []
 
 gridres = random_walk_gridres
 nsignals = len(signals)
 
 
 #get the lag limits
 xmin = np.min(np.array( [np.min(signals[i][:,0]) for i in range(nsignals)] ))
 xmax = np.max(np.array( [np.max(signals[i][:,0]) for i in range(nsignals)] ))
 laghi = (xmax-xmin)/2
 laglo = -1*laghi
 
 if (laglim == []):
  lagrange = [laglo,laghi]
 else:
  lagrange = laglim

 

  
  
 if (refdate != 0):
  a = np.float(xmax-xmin)
  xmin = pd.datetime.date(refdate)
  xmax = pd.datetime.date(refdate) + pd.Timedelta(days=np.int(a))


 
 
 print('plotlim review',xmin,xmax,laglo,laghi)
 for snow in signals:
  
  tlo = np.min(snow[:,0])
  thi = np.max(snow[:,0])
  dt  = np.mean(snow[1:,0] - snow[:-1,0])
  
  tmod = np.arange(tlo,thi+gridres,gridres)
  nmod = np.shape(tmod)[0]
  x = mf18.fitrw([snow[:,0]],[snow[:,1]],[snow[:,2]],floin=1.0/(thi-tlo),fhiin=1.0/dt,ploton=0,plot_tit='fig_myrwfit',
  dtresin=tmod,nits = monte_carlo_iterations,extra_f=extra_frequencies)
 
 
  s_model = x[4]
  s_noise = x[5]
  npoints = np.shape(snow[:,0])[0]
  
  model_itp = np.interp(snow[:,0],tmod,s_model)
  #print(np.shape(snow[:,0]),np.shape(tmod),np.shape(s_model),np.shape(model_itp))
  refined_noise = np.sqrt(1./npoints* np.sum( (model_itp - snow[:,1])**2) )
  s_noise = np.ones(nmod)*refined_noise
  
  
  
  
  
  
  time_model.append(x[3])
  signal_model.append(s_model)
  error_model.append(s_noise)
  model_lo.append(s_model - s_noise)
  model_hi.append(s_model + s_noise)
  #model_lo.append(x[6])
  #model_hi.append(x[7])
  #error_model.append(x[5])
 
 centroids_model = []
 centroids_data = []
 ccf_save_save=[]
 lag_grid_save=[]
 smodel_drive = np.array([time_model[0],signal_model[0],error_model[0]]).T
 for i in range(nsignals):
  #compute ccf between the model driver and response
  lag_grid,ccf_save,tlags_peak,rlagspeak,tlags_data = mccf.ccf_frrss(signals[0],signals[i],dt=gridres,
  centroid_frac = 0.4,
  fraction_rss=0.8,nsim = nsim)
  centroids_data.append(tlags_data)
  
  
  smodel_now = np.array([time_model[i],signal_model[i],error_model[i]]).T
  
  
  #return(smodel_drive,smodel_now)
  lag_grid,ccf_save,tlags_peak,rlagspeak,tlags_model = mccf.ccf_frrss(smodel_drive,smodel_now,dt=gridres,fraction_rss=0.8,
  centroid_frac = 0.4,
  nsim = nsim,flux_randomisation=1)
  
  #print(i,nsignals,'here',summary_plot)
  #input()
  ##sreturn(smodel_drive,smodel_now)
  
  lag_grid_save.append(lag_grid)
  ccf_save_save.append(ccf_save)
  centroids_model.append(tlags_model)  
 
 

   
 #make figures showing each plot with ccf
 #compare WITH and WITHOUT Random Walk Prior
 if (summary_plot !=''):
  ndown = nsignals
  nalong = 4
  gs1 = gridspec.GridSpec(ndown,nalong)
  gs1.update(left=0.15,right=0.85,wspace=0.1,hspace=0.01,bottom=0.1,top=0.99)
  
  
  
   
 time_model_out = []
 for i in range(nsignals):
  ccf_save = ccf_save_save[i]
  lag_save = lag_grid_save[i]
  if (refdate == 0):
   time_dat = signals[i][:,0]
   time_mod = time_model[i]
  else:
   time_dat = np.array([pd.Timestamp(refdate) + pd.Timedelta(days=signals[i][i2,0]) for i2 in range(np.shape(signals[i][:,0])[0])])
   time_mod = np.array([pd.Timestamp(refdate) + pd.Timedelta(days=time_model[i][i2]) for i2 in range(np.shape(time_model[i])[0])])
  time_model_out.append(time_mod)
   
   
  if (summary_plot !=''): 
   ax1 = plt.subplot(gs1[i,1:])
   axccf = plt.subplot(gs1[i,0]) 
   ax1.errorbar(time_dat,signals[i][:,1],signals[i][:,2],marker='o',ls='',label='data',color='r')
   #ax1.errorbar(signals[i][:,0],signals[i][:,1],signals[i][:,2])
   ax1.plot(time_mod,signal_model[i],label='model',color='b')
   #ax1.annotate(titles[i],(0.1,0.9),textcoords='axes fraction',horizontalalignment='left')
   if (titles == []):
    t = 'signals '+np.str(i)
   else:
    t = titles[i]
    
   ax1.set_ylabel(t)
   #ax1.yaxis.tick_right()
   ax1.set_yticks([])
   ax1.yaxis.set_label_position('right')
   #ax1.fill_between(time_model[i],signal_model[i]-error_model[i],signal_model[i]+error_model[i],alpha=0.2)
   ax1.fill_between(time_mod,model_lo[i],model_hi[i],alpha=0.3,color='b',label='uncertainty region')
   axccf.hist(centroids_model[i],bins=20,alpha=0.4,label='with prior',color='b')
   
   axccf.hist(centroids_data[i],bins=20,alpha=0.4,label='without prior',color='r')
   axccf.set_xlim(lagrange) 
   axccf.set_yticks([])
   
   axccf_2 = axccf.twinx()
   
   
   #return(lag_grid,ccf_save)
   
   for i2 in range(len(ccf_save)):
    axccf_2.plot(lag_save,ccf_save[i2],color='k')

   #[axccf_2.plot(lag_grid,ccf_save[i],color='k') for i in range(len(ccf_save))]
   axccf_2.set_ylim([0,1.1])
   axccf_2.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
   axccf_2.yaxis.tick_left()
   axccf_2.yaxis.set_label_position('left')
   axccf_2.set_ylabel('Correlation coefficient')
   #axccf_2.set_yticklabels(['',''])
   axccf_2.set_xlim(lagrange)
   #axccf_2.set_yticks([])
   #axccf_2.set_yticklabels([])
   ax1.set_xlim([xmin,xmax])
   if (i == 0):
    axccf.set_xticks([])
   
 if (summary_plot != ''): 
  axccf.legend() 
  ax1.legend()
  if (summary_plot=='show'):
   plt.show()
  else:
   plt.savefig(summary_plot)
  
  plt.clf()
   
 
 


 return(centroids_data,centroids_model,time_model_out,signal_model,error_model)
 
 
 
 
 
 
 
 
###test the lag_rwmod code
# 
#
#forecast = 180.0
#
#sddrive = 100.0
#meandrive = 0.0
#noise_drive = 0.8
#
#
#
#lagrange = [-100,100]
#sdecho = 1.0
#meanecho = 0.0
#noise_echo = 0.8 
#tlo = 0.0
#thi = 800.0
#dt  = 30.0
#gridres = 1.0
#lagres = 0.1
#nsim = 500
#
#
#
#
##impulse info
#lagcent = 50.0
#lagwide = 0.1
#lag = np.arange(0,lagrange[1]+lagres,lagres)
#response = np.exp(-0.5*((lag - lagcent)/lagwide)**2)
#
#
#
#
#
#
#
#
#
#datpre      = mlc.mylcgen(tlo=tlo-lagrange[1],thi=thi,dt=lagres,iseed=-1)#132423
#
#
#
#idlo  = np.where(datpre[:,0] > tlo)[0][0]
#npre     = np.shape(datpre[:,0])[0]
#
#
#datstd  = np.std(datpre[:,1])
#snow =  np.ones(npre)*noise_drive*datstd
#dat = mrs.myresample(dir='',fname=[''],dtave=dt,sampmin=20.,sampcode=3,datin=np.array((datpre[idlo:,0],datpre[idlo:,1],snow[idlo:])).T)
#ndat = np.shape(dat[:,0])[0]
#sig = dat[:,2]
#
#stdnow = np.std(dat[:,1])
#dat[:,1] = (dat[:,1] - np.mean(dat[:,1]))/stdnow*sddrive + meandrive
#dat[:,2] = dat[:,2]/stdnow*sddrive*noise_drive
#dat[:,1] = dat[:,1] + np.random.randn(ndat)*dat[:,2]
#
#
#
#x_drive = dat[:,1]
#errors_drive = dat[:,2]
#nx = np.shape(x_drive)[0]
#
#
#
#
#
#
#
#
#
##make a convolved version of the driver
#echo = mc3.mc3(datpre[:,0],datpre[:,1],lag,response)
#necho = np.shape(echo)[0]
#echostd = np.std(echo)
#errors_echo = np.ones(necho)*noise_echo*echostd
#dat_echo = mrs.myresample(dir='',fname=[''],dtave=dt,sampmin=20.,sampcode=3,datin=np.array((datpre[:,0],echo,errors_echo)).T)
#ndat = np.shape(dat_echo[:,0])[0]
#
#echosdnow = np.std(dat_echo[:,1])
#echomeannow = np.mean(dat_echo[:,1])
#
#dat_echo[:,1] = (dat_echo[:,1] - echomeannow)/echosdnow*sdecho + meanecho
#sig_echo = dat_echo[:,2]
#dat_echo[:,2] = sig_echo/echosdnow*sdecho*noise_echo
#dat_echo[:,1] = dat_echo[:,1] + np.random.randn(ndat)*sig_echo
#
#idlo_echo = np.where(dat_echo[:,0] > tlo)[0][0]
#
#
### fit modified random walk model
##nsnow = np.shape(snow)[0]
##d1 = np.zeros((nsnow,3))
##d1[:,:2] = datpre
##d1[:,2] = snow
##n_echo = np.shape(errors_echo)[0]
##d2 = np.zeros((n_echo,3))
##d2[:,0] = datpre[:,0]
##d2[:,1] = echo
##d2[:,2] = errors_echo
##
##
##id1lo = np.where((d1[:,0]>=tlo) & (d1[:,0]<=thi))[0]
##dat_new = 1.*dat
##dat_new[:,0] = dat_new[:,0] - 30.0
##
#
#
#
#
#
#
#denew = np.array(dat_echo[idlo_echo:,:])
##denew[:,0] = denew[:,0] + 50.0
#signals = [dat,denew]#[ dat, dat_new ]#[d1[id1lo,:],d2[id1lo,:]]
#titles  = ['arb','flows']
#
#a = lag_rwmod(signals,titles,summary_plot='show',monte_carlo_iterations=1)
#
#
#
##lag_grid,ccf_save,tlags_peak,rlagspeak,tlags_model = mccf.ccf_frrss(a[0],a[1])
#
##for each iteration of parameter values, compute the ratio of fourier transforms and 
##inverse to get the response function
##nits = np.shape(x_iterations[0][0,:])[0]
##response_save = []
##for it in range(nits):
## drive_ft = np.fft.fft(x_iterations[0][:,it])
## echo_ft  = np.fft.fft(x_iterations[1][:,it])
## response_ft = echo_ft/drive_ft
## response = np.fft.ifft(response_ft)
## response_save.append(response)
## 
##
##fig = plt.figure()
##ax1 = fig.add_subplot(111)
##[ax1.plot(time_model[0],rs) for rs in response_save[:100]]
###plt.savefig('fig_crazyidea_fftresponse.pdf')
###plt.clf() 
##plt.show()
#
#
#
#
# 
###save the sine and cosine value in  2d array
##cwt = np.ones((ndat,nw))
##swt = np.ones((ndat,nw))
##for iw in range(nw):
## cwt[:,iw] = np.cos(w[iw]*time)
## swt[:,iw] = np.sin(w[iw]*time)
##
##
##hnow_cc = np.tensordot(cwt.T/sig2,cwt,axes=1)#/sig2)
##hnow_sc = np.tensordot(swt.T/sig2,cwt,axes=1)#/sig2)
##
##hnow_cs = np.tensordot(cwt.T/sig2,swt,axes=1)#np.sum(cwt[:,iw]*swt[:,iwp]/sig2)
##hnow_ss = np.tensordot(swt.T/sig2,swt,axes=1)#np.sum(swt[:,iw]*swt[:,iwp]/sig2)
# 
# 
##lagcent = 34.0
##lagwide = 0.2
##lag = np.arange(0,lagrange[1]+lagres,lagres)
##response = np.exp(-0.5*((lag - lagcent)/lagwide)**2)
##
##ft = np.fft.fft(response)
##ift = np.fft.ifft(ft)
##plt.plot(np.real(ift),label='new')
##plt.plot(response,label='original')  
##plt.legend()
##plt.show()
#
#
#





#
#
#
#
#