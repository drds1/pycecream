#script to merge light curves from two or more telescopes using drw or mean/rms approximation

import numpy as np
import glob
import os
import sys
import re
import myfitrw_current as mf

ifit = 1

#/Users/ds207/Documents/standrews/sta/python/mymerge_test/ps/rm005_gi_cream_notag_ps.dat

#link to the directory containing light curves (light curves will be merged relative to the first entry
dir_data = ['/Users/ds207/Documents/standrews/sta/python/mymerge_test/sdss/',
'/Users/ds207/Documents/standrews/sta/python/mymerge_test/ps/'] 

#search string to load all light curves
data_search = ['*rm*','*rm*']

#light curve file names should have an i.d number indicating to which target they belong.
#script below scans for these and casts them into a 3 digit number (if more than 1000 targets to merge,
#change this to 4 digit)









nmerge = len(dir_data)

#lc = [[] for i in range(nmerge)]
lc  = [None]*nmerge

for i in range(nmerge):
 flist   = glob.glob(dir_data[i]+data_search[i])
 nflist  = len(flist)
 id_list = [map(int, re.findall(r'\d+', master_list[i]))[-1] for i in range(nflist)]#"%03d"%a
 
 if (i == 0):
  master_list = list(flist)
  id_master = list(id_list)
  nflist_master = len(id_master)
 
 #fit data either using mean_rms or using rw fitting or iterated RM RMS_extra fitting
 for idx in range(nflist): 
  dat = np.loadtxt(flist[idx])
  tdat = dat[:,0]
  xdat = dat[:,1]
  sdat = dat[:,2]
  
  if (ifit == 1):#rms_mean method
   stretch = np.std(xdat)
   os      = np.mean(xdat)
  elif (ifit == 2):#fit rw model to light curve before computing mean and rms of model
   mod = mf.fitrw(tdat,xdat,sdat,floin=-1,fhiin=-1,ploton=0,dtresin=0,nits = 1)
   xmod = 
   stretch = np.std(
# nfile = np
 