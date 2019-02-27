import numpy as np
import pandas as pd





#input list of pandas data frames (timestamp, quantity) with different timestamps
#rebin onto single axis and convert to numpy array

#merge dates of all other frames onto first (if no date entry matching up with first 
#date frame then input nan

def df2np(df_list,datemin=-1):
 nlist = len(df_list)
 df_combine = df_list[0]
 df_0 = df_list[0]
 n0 = len(df_0)
 output_array = np.zeros((n0,nlist))
 output_array[:,0] = df_0.values[:,1]
 
 for i in range(nlist):
  df_now = df_list[i]
  columns = df_now.columns
  cdate = columns[0]
  columns0 = df_combine.columns
  cdate0 = columns0[0]  
  df_0[cdate0] = pd.to_datetime(df_0[cdate0])
  df_now[cdate] = pd.to_datetime(df_now[cdate])
  
  #print(cdate,cdate0)
  #print(df_now)
  
  df_new = pd.merge(df_0, df_now, left_on = cdate0, right_on = cdate,how='left')
  #c_new = df_new.columns
  vals_new = df_new.values[:,-1]
  #df_combine = pd.merge(df_combine, df_new, left_on = cdate0, right_on = cdate,how='left')
  output_array[:,i] = vals_new

 
 dates = df_0.values[:,0]#df_combine.values[:,0] 
 output_data = output_array#df_combine.values[:,1:]
 
 #print(dates)
 
 #print(df_combine.columns)
 #print(output_data[:10,:])
 #print(output_data.shape)
 sums = np.sum(output_data,axis=1)
 idgood = np.where(sums == sums)[0]
 dates = dates[idgood]
 output_data = output_data[idgood,:]
 
 return(dates,output_data)