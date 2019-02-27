#return the uncertainties from cream output parameters (default returns mdot, deginc and tr slope)

import numpy as np
import glob
import myreadlines as mr


def gcp(dir,dirsub='',paridx=[2,3,4],idbi = 2./3,conlims=[15.865,50,84.135],deginc=1,idcos = 3):
 
 
 if (dirsub == ''):
  #print 'get_cream_pars   looking for cream output in '+dir+'/output_20*'
  dirnow = sorted(glob.glob(dir+'/output_20*'))[-1]
 else:
  dirnow = dir+'/'+dirsub
 
 print 'get_cream_pars   looking for cream output in '+dirnow
 
 try:
  pload = np.loadtxt(dirnow+'/outputpars.dat')
 except:
  pload = mr.myreadlines(dirnow+'/outputpars.dat')[0]
 
 ndat = np.shape(pload)[0]
 
 itlo = np.int(ndat*idbi)
 pload = pload[itlo:,paridx]
 
 
 iddeg = np.where(np.array(paridx) == idcos)[0]
 if (np.shape(iddeg)[0] > 0):
  pload[:,iddeg[0]] = np.arccos(pload[:,iddeg[0]])*180/np.pi 
 
 oplo,opmed,ophi = np.percentile(pload,conlims,axis=0)
 
 return(oplo,opmed,ophi)
 