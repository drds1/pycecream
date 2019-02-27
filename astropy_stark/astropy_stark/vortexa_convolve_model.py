import numpy as np
import astropy.convolution as apc
import vortexa_mod_randomwalk as mf18
import mylcgen as mlc
import myresample as mrs
import itertools
import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import pandas as pd
import vortexa_convert_df2lc as c2l

#enter the raw time series data.
#interpolate using random walk model,
#convolve this with kernel,
#interpolate onto response data 
#and evaluate chi-squared 
def rw_convolve(drive,response,kernel,monte_carlo_iterations=1,summary_plot = '',
xgrid_in = [],
tgrid_in = [],
noise_in = []):

 #modify kernel to make compatible with astropy convolve
 laggrid = kernel[:,0]
 dtgrid = np.mean(laggrid[1:] - laggrid[:-1])
 ngrid = np.shape(laggrid)[0]
 id0 = np.where(laggrid >= 0)[0][0]
 n_neg = ngrid - id0 - 1
 n_pad = n_neg - id0
 pad = np.zeros(n_pad)
 kernel_input = np.append(pad,kernel[:,1])
 nki = np.shape(kernel_input)[0]
 
 
 #compute the overlapping time boundaries for variance and interpolation calculation
 tol_min = np.max([np.min(drive[:,0]),np.min(response[:,0])])
 tol_max = np.min([np.max(drive[:,0]),np.max(response[:,0])])
 idxinc  = np.where((response[:,0] > tol_min) & (response[:,0] < tol_max))
 ninc    = np.shape(idxinc)[0]
 
 
 #interpolate driver onto time grid using random walk
 if (tgrid_in == []):
  tlo,thi = np.min(drive[:,0]),np.max(drive[:,0])
  tgrid = np.arange(tlo,thi+dtgrid,dtgrid)
 else:
  tgrid = tgrid_in
  tlo = tgrid[0]
  thi = tgrid[-1]
  print('be careful when manually supplying tgrid \
  must have same spacing as lag grid',tgrid[1]-tgrid[0],dtgrid)
 
 
 ndrive = np.shape(drive[:,0])[0]
 fhi = 1./dtgrid
 flo = 0.5/(thi-tlo)
 nf = (fhi - flo)/flo
 nf= min(nf,1000)
 fhi = nf*flo
 
 print('fitting driver. Number of frequencies',nf,flo,fhi)
 if (xgrid_in == []):
  
  
  x = mf18.fitrw([drive[:,0]],[drive[:,1]],[np.ones(ndrive)*0.1*np.std(drive[:,1])],floin=flo,
  fhiin=fhi,plot_tit='',dtresin=tgrid,nits = monte_carlo_iterations,extra_f=[])
  s_model = x[4]
  s_noise = x[5]
  all_models = x[12]
  response_models = np.array( [apc.convolve(all_models[:,i],kernel_input,boundary = 'extend') for i in range(monte_carlo_iterations)] ).T
  response_mean = np.mean(response_models,axis=1)
 else:
  s_model = xgrid_in
  s_noise = noise_in
  response_mean = apc.convolve(xgrid_in,kernel_input)
 print('done fitting driver',drive[0,0],drive[-1,0],ndrive)
 
 
 #interpolate model response onto data time grid
 rt_itp = np.interp(response[:,0],tgrid,response_mean)
 
 #compute the model accuracy using the raw inputted response time series
 mean_response_model = np.mean(rt_itp[idxinc])
 rms_response_model  = np.std(rt_itp[idxinc])
 mean_response_data  = np.mean(response[:,1])
 rms_response_data    = np.std(response[:,1])
 
 
 
 response_transform  = (response_mean - mean_response_model)*rms_response_data/rms_response_model + mean_response_data
 
 n_response = np.shape(response[:,0])[0]
 variance = np.sum((rt_itp[idxinc]-response[idxinc,1])**2)/ninc
 #print('model cisquared/N: ',cisq)
 
 
 
 
 
 
 
 #summary plot
 if (summary_plot != ''):
  import matplotlib.pylab as plt
  fig = plt.figure()
  ax1 = fig.add_subplot(311)
  ax1.plot(drive[:,0],drive[:,1],marker='o',ls='')
  ax1.plot(tgrid,s_model,color='b')
  #ax1.fill_between(tgrid,s_model-s_noise,s_model+s_noise,alpha=0.3,color='b')
  
  ax2 = fig.add_subplot(312)
  ax2.plot(response[:,0],response[:,1],color='r',marker='o',ls='')
  ax2.plot(tgrid,response_transform)
  #[ax2.plot(tgrid,response_models[:,i],color='r') for i in range(monte_carlo_iterations)]
  
  ax3 = fig.add_subplot(313)
  ax3.plot(kernel[:,0],kernel[:,1],color='r')
  
  ax1.set_xlim([tlo,thi])
  ax2.set_xlim([tlo,thi])
  if (summary_plot == 'show'):
   plt.show()
  else:
   plt.savefig(summary_plot)
  plt.clf()
 
 
 return(response_transform,s_model,variance,s_noise)
















#design wrapper for above function to insert into fixtures code
def rw_convolve_gridsearch(lc1,lc2,
kernel_model='gaussian',
pars=[np.arange(0.0,55.0,5.0),np.arange(1.0,5.0,1.0)],
laggrid = np.arange(0.0,50.0,0.2),
tgrid_in = []):

 #specify custom time grid
 if (tgrid_in == []):
  tlo = np.min(np.array([np.min(lc1[:,0]),np.min(lc2[:,0])]))
  thi = np.max(np.array([np.max(lc1[:,0]),np.max(lc2[:,0])]))
  dtgrid = laggrid[1] - laggrid[0]
  tgrid = np.arange(tlo,thi+dtgrid,dtgrid)
  lg = np.array(laggrid)
 else:
  tgrid = tgrid_in
  dtgrid = tgrid[1]-tgrid[0]
  lg = np.arange(laggrid[0],laggrid[-1]+dtgrid,dtgrid)
  
  
 npars = len(pars)
 pars_unique = list(itertools.product(*pars))
 npars_unique = len(pars_unique)
 var_list = []
 par_list = []
 model_list = [] 
 for ip in range(npars_unique):
  pnow = list(pars_unique[ip])
  if (kernel_model == 'gaussian'):
   centroid = pnow[0]
   width    = pnow[1]
   kern   = np.exp(-0.5*((lg - centroid)/width)**2)
   kernel = np.array([lg,kern]).T
  else:
   raise Exception('please enter a valid kernel function')
   
  #only need to calculate random walk model on first iteration (I think) 
  if (ip == 0):
   x = rw_convolve(lc1,lc2,kernel,monte_carlo_iterations=1,summary_plot = '',tgrid_in = tgrid)
  else:
   x = rw_convolve(lc1,lc2,kernel,monte_carlo_iterations=1,summary_plot = '',
   tgrid_in = tgrid,xgrid_in = s_model,noise_in = s_noise )
   
  response,s_model,variance,s_noise = x
   
  var_list.append( variance )
  par_list.append( pnow )
  model_list.append( response ) 

 var_list = np.array(var_list)
 par_list = np.array(par_list)
 #model_list = model_list
 
 
 
 
 idx_best = np.argmin(var_list)
 model_best = model_list[idx_best]
 var_best   = var_list[idx_best]
 par_best   = par_list[idx_best]




 return(var_best,par_best,model_best,var_list,par_list,model_list,s_model,s_noise) 










#make syntax correct for data frame input
def mass_convolve(meta_df,dateref=-1,forecast_period = 30,figure_title=''):

 
 nframes = len(meta_df)
 labels  = ['time series '+np.str(i) for i in range(nframes)]
 color   = ['k','r','b','cyan','purple','green']*nframes

 #identify earliest date from all the time axis and use this as the first time by default
 datemax_save = []
 if (dateref == -1):
  dates = pd.concat([meta_df[i].iloc[:,0] for i in range(nframes)])
  datemax_save.append([meta_df[i].iloc[-1,0] for i in range(nframes)])
  datemin = min(dates)
  datemax = max(dates) + pd.Timedelta(days=forecast_period)
  date_trunc = np.min(np.array(datemax_save))
  print ('Earliest dat in sample...',datemin)
 else:
  datemin = dateref[0]
  datemax = dateref[1] + pd.Timedelta(days=forecast_period)
  date_trunc = dateref[1]



 #convert data frames to light curve format
 meta_lc = []
 for i in range(nframes):
  x = c2l.convert_df_2_lightcurve(meta_df[i],reftime=datemin,noise = 0.01,normalise = 0)
  meta_lc.append(x)
 



 #fit the convolved rw model to all the light curves relative to the flows
 tlo = 0
 thi = (datemax - datemin).days + forecast_period
 dtgrid = 1.0
 timegrid = np.arange(tlo,thi+dtgrid,dtgrid)

 timegrid_dates = np.array([pd.Timestamp(datemin) + pd.Timedelta(days=tgnow) for tgnow in timegrid])     
 ntgrid = np.shape(timegrid)[0]
 meta_model = []
 laggrid = np.arange(0.0,50.0,0.2)
 for i in range(nframes):
  x = rw_convolve_gridsearch(meta_lc[i],meta_lc[0],
    kernel_model='gaussian',
    pars=[np.arange(0.0,55.0,5.0),np.arange(1.0,5.0,1.0)],
    laggrid = laggrid,
    tgrid_in = timegrid)
  var_best,par_best,model_best,var_list,par_list,model_list,s_model,s_noise = x
  noise = np.ones(ntgrid)*np.sqrt(var_best)
  a_s = [timegrid_dates,s_model,s_noise,model_best,noise,var_best,par_best]
  
  print('model means')
  print(np.mean(s_model),np.mean(model_best),np.mean(meta_lc[i][:,1]),np.mean(meta_df[i].values[:,1]))
  print('model std')
  print(np.std(s_model),np.std(model_best),np.std(meta_lc[i][:,1]),np.std(meta_df[i].values[:,1]))
  #input()
  

  meta_model.append(a_s)



 
 #produce figures for each fit
 if (figure_title != ''):
  ndown = nframes + 2
  nalong = 3
  fig = plt.figure()
  gs1 = gridspec.GridSpec(ndown,nalong)
  gs1.update(left=0.15,right=0.85,wspace=0.05,hspace=0.0,bottom=0.1,top=0.99)
  
  
  for i in range(nframes-1):
   ax1 = plt.subplot(gs1[i,:])
   #ax1.plot(meta_df[i].values[:,0],meta_df[i].values[:,1],marker='o',ls='',label=labels[i],color=color[i])
   ax1.plot(meta_df[0].values[:,0],meta_df[0].values[:,1],marker='o',ls='',label=labels[0],color=color[0])
   
   ax1.set_xlim([timegrid_dates[0],timegrid_dates[-1]])
   ax1.tick_params(axis='x',rotation=30)
   driver_time   = meta_model[i+1][0]
   driver_signal = meta_model[i+1][1]
   driver_noise  = meta_model[i+1][2]
   response_signal = meta_model[i+1][3]
   response_noise  = meta_model[i+1][4]
   if (i < nframes - 1):
    ax1.set_xticks([])
   
   ax1.plot(driver_time,response_signal,label='convolve model '+labels[i+1],color=color[i])
   
   #ax1.fill_between(driver_time,response_signal - response_noise, response_signal + response_noise,label='convolve model '+labels[i+1],color=color[i])
   
   #ax1.fill_between(driver_time,driver_signal - driver_noise, driver_signal + driver_noise,label='convolve model '+labels[i+1],color=color[i])
   print(driver_time)
   print(np.mean(driver_signal))
   print(np.mean(driver_noise))
   
   
   #if (i == 0):
   # for i2 in range(1,nframes):
     
   #  rs = meta_model[i2][3]
   #  rn = meta_model[i2][4]
   #  ax1.fill_between(driver_time,rs-rn,rs+rn,alpha=0.4) 
  
  
  
  #construct bar plot of variances for each driving model 
  ax1 = plt.subplot(gs1[ndown-1,0])
  labels = ['time series '+np.str(i) for i in range(nframes)]
  xpos = np.arange(nframes)
  vars = [meta_model[i][5] for i in range(nframes)]
  [ax1.barh([xpos[i]],[vars[i]],label=None,color=color[i]) for i in range(nframes)]
  ax1.set_title('model variances')
  ax1.set_yticks(xpos)
  ax1.set_yticklabels(labels)
  ax1.tick_params(axis='y',rotation=0)
  

  
  #construct a plot of each response function model weighted by variances
  ax1 = plt.subplot(gs1[ndown-1,2])
  ax1.set_xlabel('lag')
  ax1.set_ylabel('response')
  for i in range(nframes):
   centroid,width = meta_model[i][6]
   resp = np.exp(-0.5*((laggrid - centroid)/width)**2)
   ax1.plot( laggrid, resp, label = labels[i],color=color[i] )
  ax1.legend()
  
 
  if (figure_title == 'show'):
   plt.show()
  else:
   plt.savefig(figure_title)


 return()



#
#
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
#lagwide = 2.0
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
##dat[:,1] = dat[:,1] + np.random.randn(ndat)*dat[:,2]
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
#import myconvolve as  mc3
#import matplotlib.pylab as plt
#echo = mc3.mc3(datpre[:,0],datpre[:,1],lag,response)
#
#
##plt.plot(datpre[:,0],datpre[:,1])
##plt.plot(datpre[:,0],echo)
##plt.show()
#
#necho = np.shape(echo)[0]
#echostd = np.std(echo)
#errors_echo = np.ones(necho)*noise_echo*echostd
#dat_echo = mrs.myresample(dir='',fname=[''],dtave=dt,sampmin=2.,sampcode=3,datin=np.array((datpre[:,0],echo,errors_echo)).T)
#ndat = np.shape(dat_echo[:,0])[0]
#
#echosdnow = np.std(dat_echo[:,1])
#echomeannow = np.mean(dat_echo[:,1])
#
#dat_echo[:,1] = (dat_echo[:,1] - echomeannow)/echosdnow*sdecho + meanecho
#sig_echo = dat_echo[:,2]
#dat_echo[:,2] = sig_echo/echosdnow*sdecho*noise_echo
#dat_echo[:,1] = dat_echo[:,1] #+ np.random.randn(ndat)*sig_echo
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
#
#kernel = np.array([lag,response]).T
#
#rw_convolve(dat[:,:2],denew[:,:2],kernel,monte_carlo_iterations=1,summary_plot = '')
#
#
#
#
#
##try iterating for several different reponse functions
##time test the code
#
#import time
##impulse info
#lagcent = np.arange(0,100,10)
#lagwide = np.arange(1,10,1)
#lag = np.arange(0,lagrange[1]+lagres,lagres)
#
#icount = 0
#
#for lc in lagcent:
# for lw in lagwide:
#  response = np.exp(-0.5*((lag - 1.*lc)/1.*lw)**2)
#  rin = np.array([lag,response]).T
#  
#  tlo = time.time()
#  if (icount == 0):
#   x = rw_convolve(dat[:,:2],denew[:,:2],rin,monte_carlo_iterations=1,summary_plot = '')
#   xgrid_in = x[1]
#  else:
#   x = rw_convolve(dat[:,:2],denew[:,:2],rin,monte_carlo_iterations=1,summary_plot = '',xgrid_in = xgrid_in)
#  thi = time.time()
#  print('simulation ',icount,' run time',thi-tlo,' centroid',lc,' width',lw,' cisq',x[2]) 
#  print()
#  icount = icount + 1
#lagcent = 50.0
#lagwide = 0.1
#
#
#
#



