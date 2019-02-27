#!/usr/bin/python
#-*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt
import scipy
import CCF_interp_py3 as myccf
from scipy import stats 

""" The following is a demonstration of how to use the CCF_interp.
It requires the two light curves (also provided) sample_lc1.dat
and sample_lc2.dat. It will output three data files containing the
CCF, CCCD, and CCPD, and one plot showing the light curves, CCF,
CCCD, and CCPD, all into the current directory. 
"""

########################################
###Read in two light curves
########################################
def kateccf(lc1,lc2,lag_range=[-100,100],interp = 2.0, nsim = 500, mcmode = 0, sigmode = 0.2, plotop = '', fileop = 0,output_file_folder = './',filecent='sample_centtab.dat',filepeak='sample_peaktab.dat',filesample='sample_peaktab.dat'):

 if type(lc1) is np.ndarray:
  mjd1,flux1,err1 = lc1[:,0], lc1[:,1], lc1[:,2]
 else:
  mjd1, flux1, err1 =  np.loadtxt(lc1, unpack = True, usecols = [0, 1, 2])
 
 if type(lc2) is np.ndarray:
  mjd2,flux2,err2 = lc2[:,0], lc2[:,1], lc2[:,2]
 else: 
  mjd2, flux2, err2 =  np.loadtxt(lc2, unpack = True, usecols = [0, 1, 2])
  
 #########################################
 ##Set Interpolation settings, user-specified
 #########################################
 #lag_range = [-100, 100]  #Time lag range to consider in the CCF (days). Must be small enough that there is some overlap between light curves at that shift (i.e., if the light curves span 80 days, these values must be less than 80 days)
 #interp = 2. #Interpolation time step (days). Must be less than the average cadence of the observations, but too small will introduce noise.
 #nsim = 500  #Number of Monte Carlo iterations for calculation of uncertainties
 #mcmode = 0  #Do both FR/RSS sampling (1 = RSS only, 2 = FR only) 
 #sigmode = 0.2  #Choose the threshold for considering a measurement "significant". sigmode = 0.2 will consider all CCFs with r_max <= 0.2 as "failed". See code for different sigmodes.
 
 ##########################################
 #Calculate lag with python CCF program
 ##########################################
 tlag_peak, status_peak, tlag_centroid, status_centroid, ccf_pack, max_rval, status_rval, pval = myccf.peakcent(mjd1, flux1, mjd2, flux2, lag_range[0], lag_range[1], interp)
 tlags_peak, tlags_centroid, nsuccess_peak, nfail_peak, nsuccess_centroid, nfail_centroid, max_rvals, nfail_rvals, pvals = myccf.xcor_mc(mjd1, flux1, abs(err1), mjd2, flux2, abs(err2), lag_range[0], lag_range[1], interp, nsim = nsim, mcmode=mcmode, sigmode = sigmode)
 
 lag = ccf_pack[1]
 r = ccf_pack[0]
 
 perclim = 84.1344746    
 
 ###Calculate the best peak and centroid and their uncertainties using the median of the
 ##distributions. 
 centau = stats.scoreatpercentile(tlags_centroid, 50)
 centau_uperr = (stats.scoreatpercentile(tlags_centroid, perclim))-centau
 centau_loerr = centau-(stats.scoreatpercentile(tlags_centroid, (100.-perclim)))
 print('Centroid, error: %10.3f  (+%10.3f -%10.3f)'%(centau, centau_loerr, centau_uperr))
 
 peaktau = stats.scoreatpercentile(tlags_peak, 50)
 peaktau_uperr = (stats.scoreatpercentile(tlags_peak, perclim))-centau
 peaktau_loerr = centau-(stats.scoreatpercentile(tlags_peak, (100.-perclim)))
 print('Peak, errors: %10.3f  (+%10.3f -%10.3f)'%(peaktau, peaktau_uperr, peaktau_loerr))
         
 
 ##########################################
 #Write results out to a file in case we want them later.
 ##########################################
 if (fileop == 1):
  centfile = open(output_file_folder+filecent, 'w')
  peakfile = open(output_file_folder+filepeak, 'w')
  ccf_file = open(output_file_folder+filesample, 'w')
  for m in xrange(0, np.size(tlags_centroid)):
      centfile.write('%5.5f    \n'%(tlags_centroid[m]))
  centfile.close()
  for m in xrange(0, np.size(tlags_peak)):
      peakfile.write('%5.5f    \n'%(tlags_peak[m]))
  peakfile.close()
  for m in xrange(0, np.size(lag)):
      ccf_file.write('%5.5f    %5.5f  \n'%(lag[m], r[m]))
  ccf_file.close() 
 
 
 ##########################################
 #Plot the Light curves, CCF, CCCD, and CCPD
 ##########################################
 if (plotop != ''):
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.2, wspace = 0.1)
  
  #Plot lightcurves
  ax1 = fig.add_subplot(3, 1, 1)
  ax1.errorbar(mjd1, flux1, yerr = err1, marker = '.', linestyle = None, color = 'k', label = 'LC 1 (Continuum)')
  ax1_2 = fig.add_subplot(3, 1, 2, sharex = ax1)
  ax1_2.errorbar(mjd2, flux2, yerr = err2, marker = '.', linestyle = None, color = 'k', label = 'LC 2 (Emission Line)')
  
  ax1.text(0.025, 0.825, lc1, fontsize = 15, transform = ax1.transAxes)
  ax1_2.text(0.025, 0.825, lc2, fontsize = 15, transform = ax1_2.transAxes)
  ax1.set_ylabel('LC 1 Flux')
  ax1_2.set_ylabel('LC 2 Flux')
  ax1_2.set_xlabel('MJD')
  
  #Plot CCF Information
  xmin, xmax = -99, 99
  ax2 = fig.add_subplot(3, 3, 7)
  ax2.set_ylabel('CCF r')
  ax2.text(0.2, 0.85, 'CCF ', horizontalalignment = 'center', verticalalignment = 'center', transform = ax2.transAxes, fontsize = 16)
  ax2.set_xlim(xmin, xmax)
  ax2.set_ylim(-1.0, 1.0)
  ax2.plot(lag, r, color = 'k')
  
  ax3 = fig.add_subplot(3, 3, 8, sharex = ax2)
  #ax3.set_xlim(xmin, xmax)
  ax3.axes.get_yaxis().set_ticks([])
  ax3.set_xlabel('Centroid Lag: %5.1f (+%5.1f -%5.1f) days'%(centau, centau_uperr, centau_loerr), fontsize = 15) 
  ax3.text(0.2, 0.85, 'CCCD ', horizontalalignment = 'center', verticalalignment = 'center', transform = ax3.transAxes, fontsize = 16)
  n, bins, etc = ax3.hist(tlags_centroid, bins = 50, color = 'b')
  
  ax4 = fig.add_subplot(3, 3, 9, sharex = ax2)
  ax4.set_ylabel('N')
  ax4.yaxis.tick_right()
  ax4.yaxis.set_label_position('right') 
  #ax4.set_xlabel('Lag (days)')
  #ax4.set_xlim(xmin, xmax)
  ax4.text(0.2, 0.85, 'CCPD ', horizontalalignment = 'center', verticalalignment = 'center', transform = ax4.transAxes, fontsize = 16)
  ax4.hist(tlags_peak, bins = bins, color = 'b')
  
  plt.savefig(output_file_folder+plotop, orientation = 'landscape', bbox_inches = 'tight') 
  plt.close(plotop)
          
 return(tlags_centroid,tlags_peak,lag,r,centau,centau_uperr,centau_loerr,peaktau,peaktau_uperr,peaktau_loerr)
 
 


 
 