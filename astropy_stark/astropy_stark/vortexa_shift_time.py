import numpy as np
import pandas as pd








#shift time series 2 in time by 'shift' and 
#map onto lc1 by adding, multiplying by the appropriate means and rms
def shift_time(lc1,lc2,shift=0,renorm = 1,custom_mean = [],custom_rms = []):
 
 
 if (type(lc1) == pd.core.frame.DataFrame):
  x1         = lc1.values[:,1] 
  time1in    = pd.to_datetime(lc1.values[:,0]) 
  reftime    = time1in[0]
  time1      = np.array( [(t1p-reftime).days for t1p in time1in],dtype='float' )
  try:
   sig1 = lc1.values[:,2]
  except:
   print('no uncertainties given')
   sig1 = x1*0 + 1.0
 else:
  x1         = lc1[:,1]
  time1in    = lc1[:,0]
  reftime    = lc1[0,0]
  time1    = np.array(time1in - reftime,dtype='float')
  try:
   sig1 = lc1[:,2]
  except:
   print('no uncertainties given')
   sig1 = x1*0 + 1.0
 lc1_rms  = np.std(x1)
 lc1_mean = np.mean(x1) 
 
 
 
 
 #transform and interpolate second light curve onto first
 if (type(lc2) == pd.core.frame.DataFrame):
  lc2_rms  = np.std(lc2.values[:,1])
  lc2_mean = np.mean(lc2.values[:,1])
 else:
  lc2_rms  = np.std(lc2[:,1])
  lc2_mean = np.mean(lc2[:,1])
  
 #shift 2 by the input time lag 'shift' and transform onto lc1 scale
 if (type(lc2) == pd.core.frame.DataFrame):
  time2pre = pd.to_datetime(lc2.values[:,0])
  a = time2pre - reftime
  time2 = np.array([an.days for an in a],dtype='float') + shift
  x2 = lc2.values[:,1]
  try:
   sig2 = lc2.values[:,2]
  except:
   print('no uncertainties given')
   
 else:
  time2pre = lc2[:,0]
  time2 = np.array(time2pre + shift - reftime,dtype='float')
  x2    = lc2[:,1]
  try:
   sig2 = lc2.values[:,2]
  except:
   print('no uncertainties given')
  
 if (custom_mean != []):  
  x2new = np.array( (x2 - custom_mean[1])*custom_rms[0]/custom_rms[1] + custom_mean[0] ,dtype='float')
 else: 
  x2new = np.array( (x2 - lc2_mean)*lc1_rms/lc2_rms + lc1_mean ,dtype='float')
 

 
 #print(type(time1))
 #print(type(time2))
 #print(type(x2new))
 #print(np.mean(x2new),np.std(x2new),'mean and std',lc1_mean,lc1_rms)
 #import matplotlib.pylab as plt
 #plt.plot(x1)
 #plt.show()
 #input()
 
 x2op = np.interp(time1,time2,x2new)
 

 #try interpolating errors if present
 if (custom_mean != []):
  rms1 = custom_rms[0]
  rms2 = custom_rms[1]
 else:
  rms1 = lc1_rms
  rms2 = lc2_rms
   
 try:
  sig2op = np.interp(np.array(time1,dtype='float'),np.array(time2,dtype='float'),np.array(sig2,dtype='float'))*rms1/rms2
 except:
  print('no uncertainties given')
  sig2op = x2op*0 + 1.0
 
 
 
 
 
 
 
 #return the transformed and interpolated light curve in the same format it entered 
 #data frame or numpy array
 if (type(lc2) == pd.core.frame.DataFrame):
  lc2_op = pd.DataFrame(lc2)
  lc2_op.values[:,0] = time1in
  lc2_op.values[:,1] = x2op
  try:
   lc2_op.values[:,2] = sig2op
  except:
   print('no uncertainties given')
   
 else:
  lc2_op = np.array(lc2)
  lc2_op[:,0] = lc1[:,0]
  lc2_op[:,1] = x2op
  try:
   lc2_op[:,2] = sig2op
  except:
   print('no uncertainties given')
   sig2op = x2op*0 + 1.0
   
   
   
   
   
   
   
 #also return an array in a format
 #ready to insert into 'vortexa_combine_signal.py function'  
 n1 = np.shape(x1)[0]
 x_combine = np.zeros((n1,2))
 noise_combine = np.zeros((n1,2))
 x_combine[:,0] = x1
 noise_combine[:,0] = sig1
 x_combine[:,1] = x2op
 noise_combine[:,1] = sig2
 
  
 return(lc2_op,x_combine,noise_combine)
 
 
 
 
 
 
 
 
 
 
 