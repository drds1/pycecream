import numpy as np
import pandas as pd


#input dataframe with time axis as datetime object,
#2nd column should be the y values of the timeseries
#noise argument assigns error bars as a function of 
#light curve standard deviation (I dont know if this
#is allowed but should be fine)

def convert_df_2_lightcurve(df,reftime='auto',noise = 0.1,
normalise = 1):
 
 
 #print(df.iloc[:10,0])
 #print(df.iloc[:10,1])
 #input()
 ny,nalong = np.shape(df.values)
 
 
 if (reftime == 'auto'):
  print('setting reference time as the earliest elemnt in the time axis in day units')
  t0 = pd.Timestamp(year=2015,month=1,day=1)
 else:
  t0 = reftime 
  
 datemin = pd.Timestamp(year=2015,month=1,day=1)
 #try:
 
 cols = df.columns[0]
 dfsort = df.sort_values(by=cols)
 
 

 
 timesort = dfsort.iloc[:,0]
 datsort  = dfsort.iloc[:,1]
 
 #df = df.sort_index(axis=0)
 print(type(timesort),type(t0))
 dt = timesort - t0
 
 #t  = np.array([dtn.days() for dtn in dt])
 #something wrong above. Manually convert to days instead
 dt = np.array(dt,dtype='timedelta64[s]').astype('timedelta64[D]')
 t = dt.astype('float')

 #for i in range(ny):
 # print(t[i],df.iloc[i,1])
 y  = np.array(datsort)
 ystd = np.std(y)
 
 #if dataframe has 3 columns then one is the errors and we dont have to manufacture these
 if (nalong > 2):
  ysig = df.values[:,2]
 else:
  ysig = np.zeros(ny) + noise*ystd
 

 if (normalise == 1):
  y = (y - np.mean(y))/ystd
  ysig = ysig/ystd
 
 
  
#except:
#  raise Exception('first column must be time axis in datetime format')
  
  
 #return ny X 3 array with ny columns and rows corresponding to time, y, errorbar
 return( np.array([t,y,ysig]).T )
 
 
 