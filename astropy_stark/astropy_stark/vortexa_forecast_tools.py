#python project to load in numbers spreadsheets for each 
#exp_xxxx.numbers file in /Users/david/projects/expenses_data
#seasonal arima model article https://www.analyticsvidhya.com/blog/2018/02/time-series-forecasting-methods/

import numpy as np
import matplotlib.pylab as plt
import glob
import pandas as pd
import datetime
import statsmodels.api as sm
import vortexa_convert_df2lc as df2lc
import vortexa_polyfit as vp

#a few notes on the final argument of seasonal order
#e.g 120. This can tak evalues up to the size of the input time series
# Higher numbers will encourage the model to search longer into the past for
# seasonal variability at the expense of computation time
# if too small a number is chosen here, the model (although faster) may fail to capture
# seasonal variations

def sarima(lc1,orderin=(1,1,0),seasonal_order = 'auto',trend='c',enforce_invertibility=False,
alpha_confidence=0.32,time_forecast=10.0,diagnostic_plot=''):




 if (type(lc1) == np.ndarray):
  time   = lc1[:,0]
  signal = lc1[:,1]
  print ('sarima input is numpy array')
 elif (type(lc1) == pd.core.frame.DataFrame):
  print('sarima input is a data frame')
  a = df2lc.convert_df_2_lightcurve(lc1,reftime='auto',noise = 0.1,normalise = 0)
  print(a)
  time = a[:,0]
  signal = a[:,1]
 
 print(type(lc1))
 
 
 dt = np.mean(time[1:]-time[:-1])
 nforecast = np.int(time_forecast/dt)
 ntime = np.shape(time)[0]
 tfclo = time[ntime-1]
 tfchi = tfclo + nforecast*dt
 idx_forecast_lo = ntime - 1
 idx_forecast_hi = ntime + nforecast
 
 
 #set up the model parameters
 if (seasonal_order == 'auto'):
  so_final = np.int(0.70*ntime)
  so = (0,1,0,so_final)
 else:
  so = seasonal_order
 
 #fit the sarima model
 model   = sm.tsa.statespace.SARIMAX(endog=signal,
 order=orderin,
 seasonal_order=so,
 trend='c',
 enforce_invertibility=False)
 
 for i in range(ntime):
  print(i,time[i],signal[i])
 print('idx fc ',idx_forecast_lo,idx_forecast_hi)
 results = model.fit()
 pred    = results.get_prediction(start = ntime, end= idx_forecast_hi )
 ps      = pred.summary_frame(alpha=alpha_confidence)
 pslo    = np.array(ps['mean_ci_lower'])
 pshi    = np.array(ps['mean_ci_upper'])
 npred   = np.shape(pslo)[0]
 


 #output the forecast time series
 ymean = np.array(ps['mean'])
 ylo  = np.array(ps['mean_ci_lower'])
 yhi  = np.array(ps['mean_ci_upper'])
 nym = np.shape(ymean)[0]
 xmod = np.linspace(tfclo,tfchi,nym)

  

 #plot the time series 9can also use dweek.plot() for simple option but less customisation
 if (diagnostic_plot != ''):
  x = time
  y = signal
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.plot(x,y,label='data')
  ax1.plot(xmod,ymean,label='forecast')
  ax1.set_title(lab)
  ax1.fill_between(xmod,ylo,yhi,alpha = 0.3,label='uncertainty')
  ax1.set_xlabel('Time')
  ax1.set_ylabel('light curve')
  plt.legend(fontsize='xx-small')
  if (disgnostic_plot != 'show'):
   plt.savefig(diagnostic_plot)
   plt.clf()
  else:
   plt.show()
   
 #if the input was a data frame convert the time axis back into datetime format
 if (type(lc1) == pd.core.frame.DataFrame): 
  date_start = lc1.values[0,0]
  date_final = lc1.values[-1,0]
  time_out = np.array( [date_final + datetime.timedelta(days=xm-xmod[0]) for xm in xmod] )
 else:
  time_out = xmod 
   
   
   
   
   
   
 #compute the model Mean Absolute Percentage Error (MAPE) using cross validation
 #extract a smaller time series (smaller by length equal to the forecast)
 #use cross validation to work out the MAPE.
 frac_check = 0.9
 id_check = np.int(frac_check*ntime)
 model_check = sm.tsa.statespace.SARIMAX(endog=signal[:id_check],
 order=orderin,
 seasonal_order=so,
 trend='c',
 enforce_invertibility=False,
 enforce_stationarity=False)
 results_check = model_check.fit()
 pred_check = results_check.get_prediction(start = 0, end=ntime - 1 )
 ps_check      = pred_check.summary_frame(alpha=alpha_confidence)
 y_check = np.array(ps_check['mean'])
 
 print(id_check,ntime)
 print(y_check)
 print(signal[id_check:])
 
 mape = np.sum( np.abs(y_check[id_check:] - signal[id_check:])/signal[id_check:] )/(ntime-id_check)
 

 
 #import matplotlib.pylab as plt
 
 #plt.clf()
 #plt.plot(signal)
 #plt.plot(y_check)
 #plt.show()
 #print(mape,ntime-id_check)
 #input() 
 
 
 return(time_out,ymean,ylo,yhi,mape)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#  sarima(lc1,orderin=(1,1,0),seasonal_order = 'auto',trend='c',enforce_invertibility=False,
#alpha_confidence=0.32,time_forecast=10.0,diagnostic_plot=''):
#time_forecast argument should be in days not units
 

def polyfit(lc1,order='test',time_forecast=10.0,alpha_confidence=0.32,frac_check=0.9):
 
 
 if (type(lc1) == np.ndarray):
  time   = lc1[:,0]
  signal = lc1[:,1]
  print ('sarima input is numpy array')
 elif (type(lc1) == pd.core.frame.DataFrame):
  print('sarima input is a data frame')
  a = df2lc.convert_df_2_lightcurve(lc1,reftime='auto',noise = 0.1,normalise = 0)
  print(a)
  time = a[:,0]
  signal = a[:,1]
 
 print(type(lc1))
 
 
 dt = np.mean(time[1:]-time[:-1])
 nforecast = np.int(time_forecast/dt)
 ntime = np.shape(time)[0]
 tfclo = time[ntime-1]
 tfchi = tfclo + nforecast*dt
 idx_forecast_lo = ntime - 1
 idx_forecast_hi = ntime + nforecast
  
 #grid search to find the optimum number of parameters
 if (order == 'test'):
  best_cisqred,best_aic,best_bic = vp.fit_search(time,signal,maxorder=6)
  idx = best_cisqred
 else:
  idx = order
  
 
 #perform the fit
 time_plus_fc = np.arange(time[0],time[-1]+time_forecast+dt,dt)
 yg_med,yg_lo,yg_hi,cov,cisq,cisq_red,bic,aic = \
 vp.fit(time,signal,np.ones(ntime),order=3,xgrid=time_plus_fc,confidence=alpha_confidence,nits=2000)
 
 
 
 
 #if the input was a data frame convert the time axis back into datetime format
 if (type(lc1) == pd.core.frame.DataFrame): 
  date_start = lc1.values[0,0]
  date_final = lc1.values[-1,0]
  time_out = np.array( [date_final + datetime.timedelta(days=xm-time[-1]) for xm in time_plus_fc[ntime:]] )
 else:
  time_out = time_plus_fc[ntime:]
   
    
 y_out   = yg_med[ntime:]
 ylo_out = yg_lo[ntime:]
 yhi_out = yg_hi[ntime:]
 
 
 
 
 
 #evaluate the mape
 ncheck = np.int(frac_check*ntime)
 time_check_fc = np.arange(time[0],time[-1]+dt,dt)
 check_med,check_lo,check_hi,check_cov,check_cisq,check_cisq_red,check_bic,check_aic = \
 vp.fit(time[:ncheck],signal[:ncheck],np.ones(ncheck),order=idx,xgrid=time,confidence=alpha_confidence,nits=2000)
 mape = np.sum( np.abs( (signal[ncheck:] - check_med[ncheck:])/signal[ncheck:]) )/ (ntime-ncheck + 1) 
  
  
 
 
 
 
 
 return(time_out,y_out,ylo_out,yhi_out,mape)
 












