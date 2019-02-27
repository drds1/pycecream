import numpy as np
import matplotlib.pylab as plt
import glob
import pandas as pd
import datetime
import statsmodels.api as sm



def sarima(lc1,orderin=(1,1,0),seasonal_order = (0,1,0,120),trend='c',enforce_invertibility=False,
alpha_confidence=0.32,time_forecast=10.0,diagnostic_plot=''):

 time   = lc1[:,0]
 signal = lc1[:,1]
 dt = np.mean(time[1:]-time[:-1])
 nforecast = np.int(time_forecast/dt)
 ntime = np.shape(time)[0]
 tfclo = time[ntime-1]
 tfchi = time_forecast_lo + nforecast*dt
 idx_forecast_lo = ntime
 idx_forecast_hi = ntime + nforecast
 
 #fit the sarima model
 model   = sm.tsa.statespace.SARIMAX(endog=signal,
 order=order,
 seasonal_order=seasonal_order,
 trend='c',
 enforce_invertibility=False)
 results = model.fit()
 pred    = results.get_prediction(start = idx_forecast_lo, end= idx_forecast_hi )
 ps      = pred.summary_frame(alpha=alpha_confindence)
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
   
   
 return(xmod,ymean,ylo,yhi)