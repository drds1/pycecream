import json, os, sys, logging
import pandas as pd, numpy as np
import sqlalchemy
from tqdm import tqdm
from pandas.io.sql import SQLTable
import matplotlib.pylab as plt
import os
import matplotlib.gridspec as gridspec
import vortexa_ccf_simple as mccf
import vortexa_convert_df2lc as c2l
import my_kateccf_py3 as mc3
import vortexa_forecast_tools as vs
import datetime
import vortexa_shift_time as vst
import vortexa_combine_signal as vcs
import vortexa_improve_lag_idea as rwlag

#format should be a list of data frames where each frame
#has 2 or 3 axes (1st is the time as a date time object, 2nd are the flux values,
#3rd (optional) are the uncertainties . Good to have these for CCF otherwise a 
#certain amount of guess work will ensue.
#20/11/2018 added forecast period. If forecast_period = , do no forecasting
#...else perform a forecast of 'forecast_period' days ahead

def vortexa_lag_analysis(df_in,df_titles='',figure_title = 'vortexa_lags_summary.pdf',
dateref = -1,laglim=[-50,50],gridfig=[16,12],
figures=[['driving time series','','Flow (tonnes)'],
['response time series','Date','Arbitrage'],
['ccf simple','',''],
['mean lags','Arbitrage forecast month','lag (days)'],
['ccf frrss','lag (days)',''],
['power spectrum','P(f)','frequency (cycles/day)']
],
figure_positions_x=[[0,12],[0,12],[0,3],[4,8],[9,12],[0,3]],
figure_positions_y=[[0,3],[4,7],[9,12],[9,12],[9,12],[13,16]],
sigmode = 0.2,
lagstat='peak',
nsim = 5000,
histbins = 20,
global_titles = ['']*100,
forecast_period = 0.0,
datemin_force = '',
datemax_force = '',
normalise = False
):
 
 figheads = [figures[i][0] for i in range(len(figures))]
 try:
  id_drive = figheads[0].index('driving time series')
 except:
  id_drive = -1
  
 try:
  id_echo = figheads.index('response time series')
 except:
  id_echo = -1
 
 try:
  id_ccf = figheads.index('ccf simple')
 except:
  id_ccf = -1
 
 try:
  id_lag = figheads.index('mean lags') 
 except:
  id_lag = -1
  
 try: 
  id_ccf_frrss = figheads.index('ccf frrss') 
 except:
  id_ccf_frrss = -1

 try: 
  id_pspec = figheads.index('power spectrum') 
 except:
  id_pspec = -1
  

  
  
  
  
  
 #initialise lists to store the uncertainties from the forecasts
 forecast_times = []
 forecast_uncerts = []
 forecast_uncerts_incarb = []
   
   
 #make plots of OB vs time and flows vs time for each region
 ndown  = gridfig[0] 
 nalong = gridfig[1]
 
  
 #insert list of data frames where each data frame has a time,value, uncertainty axis
 #measure lags relative to first
 nframes = len(df_in)
 necho   = nframes - 1
 color=['k','r','b','cyan','purple','orange','gree','skyplue','magenta']*nframes

 #must have at least 2 data frames to perform lag analysis
 if (nframes < 2):
  print('Must have at least 2 data frames to perform lag analysis in vortexa_lag_figures.py')
  return()

  
 #convert dates to datetime format
 df = []#list(df_in)
 for i in range(nframes):
  time = pd.to_datetime(df_in[i].iloc[:,0])
  dat  = df_in[i].iloc[:,1]
  
  if (normalise == True):
   dat = (dat - np.mean(dat))/np.std(dat)
  
  df.append( pd.concat([time,dat],axis=1) )
  #df[i].iloc[:,0] = pd.to_datetime(df[i].iloc[:,0])
  if (( datemin_force != datetime.datetime) & (datemax_force !='')):
   df[i]= df[i][(df[i].iloc[:,0] > pd.Timestamp(datemin_force)) & (df[i].iloc[:,0] < pd.Timestamp(datemax_force))]
   print(min(df[i].iloc[:,0]),max(df[i].iloc[:,0]),'min max forced dates')
 #save the combine forecast array using time shift
 
 
 #print(df[0].iloc[:10,:])
 #print(df[1].iloc[:10,:])
 #input('test break here')
 #print(df_in[0].iloc[:10,:])
 #print(df_in[1].iloc[:10,:])
 #input('test break here') 

 combine_x = []
 combine_noise = []
 combine_r = [1.]
 mape_save = []
  
 #perform forecasting if required
 lc_forecast = [[] for i in range(nframes)]
 if (forecast_period > 0):
  idx = 0
  for dfnow in df_in:
   
   
   a_s = vs.sarima(dfnow,orderin=(1,1,0),seasonal_order = 'auto',trend='c',enforce_invertibility=False,
   alpha_confidence=0.32,time_forecast=forecast_period,diagnostic_plot='')
   
   
   as0 = a_s[0]
   
   a_s = vs.polyfit(dfnow,order='test',time_forecast=forecast_period,alpha_confidence=0.32,frac_check=0.9)
   print(dfnow.values[-1,0],'last time')
   print(as0,'sarima')
   print(a_s[0],'poly')
   #input()
   mape_save.append(a_s[4])
   a = a_s[:4]
   lc_forecast[idx] = a
   
   idx = idx + 1

   
 #identify earliest date from all the time axis and use this as the first time by default
 datemax_save = []
 if (dateref == -1):
  dates = pd.concat([df[i].iloc[:,0] for i in range(nframes)])
  datemax_save.append([df[i].iloc[-1,0] for i in range(nframes)])
  datemin = min(dates)
  datemax = max(dates) + datetime.timedelta(days=forecast_period)
  date_trunc = np.min(np.array(datemax_save))
  print ('Earliest dat in sample...',datemin)
 else:
  datemin = dateref[0]
  datemax = dateref[1] + datetime.timedelta(days=forecast_period)
  date_trunc = dateref[1]
 #truncate the figures to make both the arb and flows run to the same time
 

 


 
 #if no titles supplied then make up here
 if (df_titles == ''):
  dataframe_titles = [np.str('Dataframe '+np.str(i)) for i in range(nframes)]
 else:
  dataframe_titles = df_titles
 
 #initialise figures
 gs1 = gridspec.GridSpec(ndown,nalong)
 gs1.update(left=0.15,right=0.85,wspace=0.05,hspace=0.0,bottom=0.1,top=0.99)
  
 
 #deal with the driving time series
 df_drive_time = df[0].iloc[:,0]
 df_drive_x    = df[0].iloc[:,1]
 if (len(df_drive_time)==0):
  print('first list element (driving time series) of input data frame cannot be empty')
  return()
  
 #return(df)
 lc_drive = c2l.convert_df_2_lightcurve(pd.concat([df_drive_time,df_drive_x],axis=1),reftime=datemin,noise = 0.01,normalise = 1)
 lc_drive = np.array(lc_drive,dtype=float)
 
 
 if (id_drive > -1):
  fig_x = slice(figure_positions_x[id_drive][0],figure_positions_x[id_drive][1],1)
  fig_y = slice(figure_positions_y[id_drive][0],figure_positions_y[id_drive][1],1)
  ax1 = plt.subplot(gs1[fig_y,fig_x])
  #plot the error bars if present
  try: 
   df_sig_x = df[0].iloc[:,2]
   ax1.errorbar(df_drive_time,df_drive_x,df_sig_x,
   label=None,color=color[0])
  except:
   pass
  
  
  ax1.plot(np.array(df_drive_time),np.array(df_drive_x),label=None,color=color[0],marker='o',ls='',markersize=3)
  print('min max  drive info',np.min(df_drive_x),np.max(df_drive_x))
  #plot the forecast if present
  fc = lc_forecast[0]
  if (len(fc) > 0):
   time_out,ymean,ylo,yhi = fc
   ax1.plot(time_out,ymean,label=None,ls='--',color=color[0])
   ax1.fill_between(time_out,ylo,yhi,alpha=0.2,label='forecast (historic flows only)',color=color[0])
   #for it in range(np.shape(df_drive_time)[0]):
   # print(it,df_drive_time[it],df_drive_x[it])
   #for it in range(np.shape(time_out)[0]):
   # print(it,'fc',time_out[it],ymean[it])
   #input()
   for it in range(np.shape(time_out)[0]):
     print('fc drive',it,time_out[it],ylo[it],ymean[it],yhi[it])
   
   ##indicate start of forecast period
   #ax1_b = ax1.twiny()
   #xl1 = list(ax1.get_xlim())
   yl1 = list(ax1.get_ylim())
   #ax1_b.set_xlim(xl1)
   #ax1_b.set_ylim(yl1)
   #ax1_b.set_xticklabels(['forecast start'])
   #ax1_b.set_yticklabels([])
   #ax1_b.set_yticks([])
   #ax1_b.set_xticks([time_out[0]])
   #ax1_b.plot([time_out[0]]*2,yl1,label=None,color='k',ls='--')
   #ax1_b.xaxis.set_tick_params(labelbottom='on',labeltop='off')
   ax1.plot([time_out[0]]*2,yl1,label=None,color='k',ls='--') 
  
  
  ax1.set_xlabel(figures[id_drive][1])
  ax1.set_ylabel(figures[id_drive][2])
  ax1.legend(fontsize='x-small')
  #set the global title
  ax1.set_title(global_titles[0])
 
 
 
 
 
 
 
 
 
 #deal with each of the response time series'
 if (id_echo > -1):
  fig_x = slice(figure_positions_x[id_echo][0],figure_positions_x[id_echo][1],1)
  fig_y = slice(figure_positions_y[id_echo][0],figure_positions_y[id_echo][1],1)
  ax2 = plt.subplot(gs1[fig_y,fig_x])
  ax2.set_xlabel(figures[id_echo][1])
  ax2.set_ylabel(figures[id_echo][2])
  ax2.set_title(global_titles[id_echo])

 
 #initialise ccf example
 if (id_ccf > -1):
  fig_x = slice(figure_positions_x[id_ccf][0],figure_positions_x[id_ccf][1],1)
  fig_y = slice(figure_positions_y[id_ccf][0],figure_positions_y[id_ccf][1],1)
  ax3 = plt.subplot(gs1[fig_y,fig_x])
  #ax3.set_xlabel('lag (days)')
  #ax3.set_ylabel('cross correlation function') 
  ax3.set_xlabel(figures[id_ccf][1])
  ax3.set_ylabel(figures[id_ccf][2])

  ax3.set_title(global_titles[id_ccf])

 #initialise mean lag plot
 if (id_lag > -1):
  fig_x = slice(figure_positions_x[id_lag][0],figure_positions_x[id_lag][1],1)
  fig_y = slice(figure_positions_y[id_lag][0],figure_positions_y[id_lag][1],1)  
  ax4 = plt.subplot(gs1[fig_y,fig_x])
  ax4.set_xlabel('Arbitrage forecast month')
  ax4.set_ylabel('lag (days)')
  ax4.set_title(global_titles[id_lag])
  ax4.set_xlabel(figures[id_lag][1])
  ax4.set_ylabel(figures[id_lag][2])

 

 #initialise ccf frrrss histogram
 if (id_ccf_frrss > -1):
  fig_x = slice(figure_positions_x[id_ccf_frrss][0],figure_positions_x[id_ccf_frrss][1],1)
  fig_y = slice(figure_positions_y[id_ccf_frrss][0],figure_positions_y[id_ccf_frrss][1],1)
  ax_frrss = plt.subplot(gs1[fig_y,fig_x]) 
  ax_frrss.set_title(global_titles[id_ccf_frrss])
  ax_frrss.set_xlabel(figures[id_ccf_frrss][1])
  ax_frrss.set_ylabel(figures[id_ccf_frrss][2])

 
 
 #define lists to save the lag info
 tsave = []
 lagsave = []
 r_coefs = []
 
  
 lc_echo_save = []
 rwlag_save=[np.array([np.array(lc_drive[:,0]),np.array(df_drive_x),np.ones(np.shape(df_drive_x)[0])]).T]
 for i in range(necho):
  df_echo_time = df[i+1].iloc[:,0]
  df_echo_x = df[i+1].iloc[:,1]
  
  if (id_echo > -1):
   ax2.plot(np.array(df_echo_time),np.array(df_echo_x),color=color[i+1],label=None,ls='',marker='o',markersize=3) 
   try:#look for errorbars in echo light curves
    df_echo_sig = df[i+1].iloc[:,2]
    ax1.errorbar(df_echo_time,df_echo_x,df_echo_sig,
    label=None,color=color[i+1])
   except:
    pass
  
   #plot the forecast if present
   fc = lc_forecast[i+1]
   if (len(fc) > 0):
    time_out,ymean,ylo,yhi = fc
    ax2.plot(time_out,ymean,label=None,ls='--',color=color[i+1])
    ax2.fill_between(time_out,ylo,yhi,alpha=0.2,label=None,color=color[i+1])
    for it in range(np.shape(time_out)[0]):
     print('fc',i,it,time_out[it],ylo[it],ymean[it],yhi[it])
   
  #perform the ccf function and plot to the figure
  lc_echo = c2l.convert_df_2_lightcurve(df[i+1],reftime=datemin,noise = 0.01,normalise = 1)
  lc_echo = np.array(lc_echo,dtype=float)
  lc_echo_save.append(lc_echo)
  
  try:
   tccf,ccf,lo,mean,hi,idlo,idmax,idhi = mccf.ccf_simp(lc_drive,lc_echo,dt=1.0,r_crit = 0.8,laglim=laglim)
  except:
   print('problem with ccf opperation')
   tccf = np.array([0])
   ccf  = np.array([0])
   lo = 0
   mean = 0
   idlo = 0
   idhi = 0
   idmax = 0
   
   
  
  #attempt the frrss method (first try fake errorbars)
  if (id_ccf_frrss > -1):
   cadence_drive = np.mean(lc_drive[1:,0]-lc_drive[0:-1,0])
   cadence_std_drive = np.std(lc_drive[1:,0]-lc_drive[0:-1,0])
   cadence_echo  = np.mean(lc_echo[1:,0]-lc_echo[0:-1,0])
   cadence_std_echo = np.std(lc_echo[1:,0]-lc_echo[0:-1,0]),np.std(lc_echo[1:,0]-lc_echo[0:-1,0])
   print('frrss cadence drive (mean rms) ',cadence_drive, cadence_std_echo)
   print('frrss cadence response (mean rms) ',cadence_echo, cadence_std_echo)
   
   overlap_min = max(np.min(lc_drive[:,0]),np.min(lc_echo[:,0]))
   overlap_max = min(np.max(lc_drive[:,0]),np.max(lc_echo[:,0]))
   overlap_period = overlap_max-overlap_min
   
   lolim =  -overlap_period/4
   hilim =  overlap_period/4
   lag_range = [lolim,hilim]
   gridres = 0.5*(cadence_drive + cadence_echo)/2
   print('frrss overlap period (days) ',overlap_period,
   '    lower and upper lag limits ',lolim,hilim,
   '     grid resolution ',gridres)
   
   
   try:
    lag_grid,ccf_save,tlags_peak,rlagspeak,tlags_centroid = mccf.ccf_frrss(lc_drive,lc_echo,dt=gridres,fraction_rss=0.8,nsim = nsim)
   except:
    tlags_centroid = np.array([0])
    ccf_save = np.array([0])
    tlags_peak = np.array([0])
    rlagspeak = np.array([0])
    tlags_centroid = np.array([0])
    print('frrss ccf failed')
   
   
   if (lagstat == 'peak'):
    tlags_centroid = tlags_peak
    
   lagrange = hilim-lolim
   lag_margain = 0.1*lagrange
   idinc = np.where((tlags_centroid > lolim + lag_margain) & (tlags_centroid < hilim - lag_margain))[0]
   tlags_centroid = tlags_centroid[idinc]  
   r_coefs.append(np.array(rlagspeak[idinc]))
   print (np.shape(tlags_centroid),'cent shape')
   

   hist = ax_frrss.hist(tlags_centroid,bins=histbins,alpha=1.0,color=color[i+1],label=None,normed=False,histtype='step')
   histmax_ccf = hist[0].max()
   
   print(np.shape(lc_echo),np.shape(lc_drive),np.shape(tlags_centroid))
   if (np.shape(tlags_centroid)[0] == 0):
    print('no centroids calculated - porbably no overlap between time series')
    lo,med,hi = [0,0,0]
   else:
    lo,med,hi = np.percentile(tlags_centroid,[15.865,50,84.135])
   ax_frrss.plot([lolim+lag_margain]*2,[0,histmax_ccf],color=color[0],ls='--',lw=4,label=None)
   ax_frrss.plot([hilim-lag_margain]*2,[0,histmax_ccf],color=color[0],ls='--',lw=4,label=None)
   

   rwlag_save.append( np.array([np.array(lc_echo[:,0]),np.array(df_echo_x),np.ones(np.shape(df_echo_x)[0])]).T )
  
  
  
  
  
  
  
  
  
   
   
   #combine the driving time series and arbitrage weighted by the correlation coefficient
   custom_mean = [np.mean(df[0].iloc[:,1]),np.mean(df[i+1].iloc[:,1])]
   custom_rms =  [np.std(df[0].iloc[:,1]),np.std(df[i+1].iloc[:,1])]
   print('combination_mean...',custom_mean)
   print('combination_rms...',custom_rms)
   fc = lc_forecast[i+1]
   if (len(fc) > 0):
    time_out,ymean,ylo,yhi = lc_forecast[0]
    t0   = pd.to_datetime(time_out)
    y0   = ymean
    sig0 = (yhi - ylo)/2
    t1,y1,ylo1,yhi1 = fc
    sig1 = (yhi1 - ylo1)/2
    df1 = pd.DataFrame([t0,y0,sig0]).T
    df2 = pd.DataFrame([t1,y1,sig1]).T
    crap,x_combine,noise_combine = vst.shift_time(df1,df2,shift=0,renorm = 1,custom_mean = custom_mean,custom_rms = custom_rms)#return(lc2_op,x_combine,noise_combine)
    combine_r.append( np.mean(np.array(rlagspeak[idinc])) )
    if (i + 1 == 1):
     combine_x.append(x_combine[:,0])
     combine_noise.append(noise_combine[:,0])
     forecast_uncerts.append(noise_combine[:,0])
    combine_x.append(x_combine[:,1])
    combine_noise.append(noise_combine[:,1])
    
 
    if (i == necho - 1):
     
     #add the updated forecast to the driving light curve plot
     combine_x = np.array(combine_x).T
     combine_noise = np.array(combine_noise).T
     combine_r = np.array(combine_r)
     signal_out,signal_noise = vcs.comb_signal(combine_x,combine_noise,combine_r)#return(signal_out,signal_noise)
     forecast_uncerts_incarb.append(signal_noise)
     forecast_times = pd.to_datetime(t0)  

   
   if (id_ccf > -1):
    ax3.plot(tccf,ccf,color=color[i+1],label = dataframe_titles[i+1]+' ccf')
    ax3.plot([tccf[idlo],tccf[idhi]],[ccf[idlo],ccf[idhi]],label=None,color=color[i+1])
  
   #make summary plot of mean lags
   if (id_lag > -1):
    ax4.errorbar([i],[med],[[med- lo],[hi - med]],marker='o',label = dataframe_titles[i+1],color=color[i+1])
    tsave.append(i)
    lagsave.append(med)
  
  
  
  
  
  
  
  
  
  
  
 #apply the random walk lag detector
 #return(centroids_data,centroids_model,time_model,signal_model,error_model)
 
 try:
  rwlag_op = rwlag.lag_rwmod(rwlag_save,titles=[],
  monte_carlo_iterations=1,
  extra_frequencies = [],
  random_walk_gridres = 1.0,
  summary_plot = '',#figure_title[:-4]+'_rwlag.pdf',
  refdate=datemin)
  centroids_model,time_model,signal_model,error_model = rwlag_op[1],rwlag_op[2],rwlag_op[3],rwlag_op[4]
 except:
  centroids_model = [ np.array([0,0]) for it in range(nframes) ]
  time_model = [ np.array([datemin,datemax]) for it in range(nframes) ]
  signal_model = [ np.array([0,0]) for it in range(nframes) ]
  error_model = [ np.array([0,0]) for it in range(nframes) ]
  
  

 
 for i in range(nframes):
  tm_now = time_model[i]
  sm_now = signal_model[i]
  no_now = error_model[i]
  if (i == 0):
   ax1.plot(tm_now,sm_now,color=color[i],label=dataframe_titles[i]+' model')
   ax1.fill_between(tm_now,sm_now-no_now,sm_now+no_now,color=color[i],label=None,alpha=0.3)
  else:
   ax2.plot(tm_now,sm_now,color=color[i],label=dataframe_titles[i]+' model')
   ax2.fill_between(tm_now,sm_now-no_now,sm_now+no_now,color=color[i],label=None,alpha=0.3)
 for i in range(necho): 
  print('model info')
  print(tm_now)
  print(sm_now)
  print(no_now)
  try:
   ax_frrss.hist(centroids_model[i],bins=histbins,alpha=1.0,color=color[i+1],
   histtype='step',label=dataframe_titles[i],normed=False)
  except:
   print('problem computing frrss histogram')

    
  
  
  
 
 if (id_lag > -1):
  #ax4.plot(tsave,lagsave,color=color[0],label=None)
  ax4.plot(tsave,-30*np.array(tsave),color=color[0],ls='--',label='arbitrage expected lag')
  ax4.legend(fontsize='xx-small')
 
 
 if (id_ccf_frrss > -1):
  ax_frrss.legend(fontsize='xx-small')
  yl = list(ax_frrss.get_ylim())
  idx = 0
  ax_frrss.tick_params(labelleft=False)
  for lagnow in lagsave:
   #ax_frrss.plot([lagnow]*2,yl,color=color[idx+1],label=None)
   idx=idx+1
 
 #compute power spectrum
 tx = lc_drive[:,0]
 tmean = np.mean(tx[1:]-tx[:-1])
 x  = lc_drive[:,1]
 nx = np.shape(tx)[0]
 ft = np.fft.fft(x)
 pspec = np.real(np.conj(ft)*ft)
 freq = np.fft.fftfreq(nx,tmean)
 nf = np.shape(freq)[0]
 nf2 = np.int(nf/2)
 fhi = freq[nf2]
 flo = 1.e-10
 idinc = np.where((freq>flo) & (freq<fhi) & (pspec>1.e-10))[0]
 
 if (id_pspec > -1):
  fig_x = slice(figure_positions_x[id_pspec][0],figure_positions_x[id_pspec][1],1)
  fig_y = slice(figure_positions_y[id_pspec][0],figure_positions_y[id_pspec][1],1)
  axps = plt.subplot(gs1[fig_y,fig_x]) 
  #axps.set_xlabel('frequency (cycles/day)')
  #axps.set_ylabel('P(f) Power spectrum')
  axps.set_xlabel(figures[id_pspec][1])
  axps.set_ylabel(figures[id_pspec][2])

  
  
  axps.plot(freq[idinc],pspec[idinc],color=color[0],label=dataframe_titles[0])
  
  for i in range(necho):
   lcn = lc_echo_save[i]
   dtmean = np.mean(lcn[1:,0]-lcn[:-1,0])
   nechotemp = np.shape(lcn[:,0])[0]
   freq = np.fft.fftfreq(nechotemp,dtmean)
   idinc = np.where((freq>flo) & (freq<fhi))[0]
   nf = np.shape(freq)[0]
   nf2 = np.int(nf/2)
   fhi = freq[nf2]
   ft = np.fft.fft(lcn[:,1])
   pspec = np.real(np.conj(ft)*ft)[:nf2]

   
   axps.plot(freq[idinc],pspec[idinc],color=color[i+1],label=dataframe_titles[i+1])
  
  

  
  ylim = list(axps.get_ylim())
  xlim = list(axps.get_xlim())
  ps_sort = np.sort(pspec)
  
  ylim[0] =ps_sort[2]#max([10**-10,ylim[0]])


  axps.set_ylim(ylim)



  xrefs = [1./365,1./31.,1./7]
  xreflab =['1 year','1 month','1 week']
  #
  axps2 = axps.twiny()
  axps2.set_xlim(xlim)
  axps2.set_ylim(ylim)
  
  axps2.set_xticks(xrefs)
  axps2.set_xticklabels(xreflab)
  
  [axps2.plot([xn]*2,ylim,ls='--',color='k',lw=1,label=None,zorder=0) for xn in xrefs]
  axps.set_xscale('log')
  axps.set_yscale('log')
  #axps2.set_xscale('log')
  
  #try:
  # axps.set_xscale('log')
  # axps.set_yscale('log')
  # axps2.set_xscale('log')
  #except:
  # print('cannot log-scale power spectrum')
 
 if (id_drive > -1):
  if (forecast_period > 0):
   ax1.plot(t0,signal_out,color='purple',label=None)
   ax1.fill_between(t0,signal_out-signal_noise,signal_out+signal_noise,color='purple',label='forecast (inc arb model)',alpha=0.3)
  ax1.legend()
  ax1.set_xlim([datemin,datemax])
 if (id_echo > -1):
  if (forecast_period > 0):
   yl2 = list(ax2.get_ylim())
   ax2.plot([datemax - datetime.timedelta(days=forecast_period)]*2,yl2,label=None,color='k',ls='--')
  ax2.legend()
  ax2.tick_params(axis='x',rotation=30)
  ax2.set_xlim([datemin,datemax])
  if (forecast_period > 0):
   ax2_b = ax2.twiny()
   ax2_b.plot([datemax - datetime.timedelta(days=forecast_period)]*2,yl2,label=None,color='k',ls='--')
   ax2_b.set_xlim([datemin,datemax])
   ax2_b.set_ylim(yl2)
   ax2_b.set_xticklabels(['forecast start'])
   #ax2_b.set_yticklabels([])
   #ax2_b.set_yticks([])
   ax2_b.set_xticks([datemax - datetime.timedelta(days=forecast_period)])
   ax2_b.xaxis.set_tick_params(labelbottom='off',labeltop='on')
   ax2_b.set_xlim([datemin,datemax])
   ax2.set_xlim([datemin,datemax])
   ax1.set_xlim([datemin,datemax])
   ax1.plot([datemax - datetime.timedelta(days=forecast_period)]*2,yl1,label=None,color='k',ls='--')
   #ax1.plot([time_out[0]]*2,yl1,label=None,color='k',ls='--') 
  ax1.set_xticklabels([])
  ax1.set_xlabel('')
  #print(datemax,forecast_period)
  #print(ax2_b.get_xlim())
  #print(ax2.get_xlim())
  #print(ax1.get_xlim())
  #input()
 if (id_ccf > -1):
  pass
  #ax3.legend()
 
 if (id_lag > -1):
  #ax4.legend()
  pass
 
 
 
 
 
 plt.savefig(figure_title)
 
 #print(global_titles[0])
 #print( np.array(forecast_uncerts), 'mean',np.mean(forecast_uncerts) )
 #print( np.array(forecast_uncerts_incarb),'mean',np.mean(forecast_uncerts_incarb) )
 #input()
 ax1.set_xlim([datemin,datemax])
 ax2.set_xlim([datemin,datemax])
 if (forecast_period > 0):
  ax2_b.set_xlim([datemin,datemax])
 return(r_coefs,forecast_times,np.array(forecast_uncerts),np.array(forecast_uncerts_incarb),
 np.array(mape_save))
 #





 
 #optional output each light curve as a text file for further more specific analysis
 #...
 
#  
#import mylcgen as mlc
#flows = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=1100,dt=30.0,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
#flowmean = np.mean(flows[:,1])
#flowrms  = np.std(flows[:,1])
#time = flows[:,0]
#
#tref = pd.Timestamp(year=2015,month=1,day=1)
#timeflows = [tref + pd.Timedelta(days=fl) for fl in time]
#
#timelag =42.0
#timeflows2 = [tref + pd.Timedelta(days=fl+timelag) for fl in time]
#
#
#
#meta_df = [ pd.concat([pd.DataFrame(timeflows),pd.DataFrame(flows[:,1])],axis=1) ]
#meta_df.append( pd.concat([pd.DataFrame(timeflows2),pd.DataFrame(flows[:,1])],axis=1)  )
#a_op = vortexa_lag_analysis(meta_df,figure_title = './figures/test_fake.pdf')
#
# 
  
   
 
 