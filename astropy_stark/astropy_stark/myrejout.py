import numpy as np


#Kates outlier rejection criteria. Need to apply this to sdss light curves prior to merging
def myrejout_kate(tin,yin,sigyin,itmax = 5,tel='spec',epochdelete=[]):
 
 t = 1.*tin
 y = 1.*yin
 sigy = 1.*sigyin
 #deal with indicees to delete
 if (len(epochdelete) > 0):
  edrev = sorted(epochdelete,reverse=True)
  for idx in edrev:
   t = np.delete(t,idx)
   y = np.delete(y,idx)
   sigy = np.delete(sigy,idx)
 
 
 
 if (tel == 'phot'):
  ynew = 1.*y
  tnew = 1.*t
  signew = 1.*sigy
  nold = np.shape(ynew)[0]
  nnew = -666.
  it = 1
  
  
  while ((nnew != nold) and (it <= itmax)):
   
   nold = np.shape(ynew)[0]
   print nold
   
   #remove too-tiny fluxes
   idxinc_lo = np.where(np.abs(ynew) > 1.)[0]
   print len(idxinc_lo),'len(idxinc_lo)'
   #recalculate nmad, calculate median errors and stuff
   median = np.median(ynew)
   diff = np.abs(ynew - median)
   mad = np.median(diff)
   nmad = mad/0.6745
   med_err = np.median(signew)
   #print 'median,mad,nmad,med_err'
   #print median,mad,nmad,med_err
   #remove too-big fluxes
   idxinc_big = np.where(np.abs(ynew) < 5.*nmad)[0]
   print len(idxinc_big), 'len(idxinc_big)'
   
   #remove points with tiny errors
   idxinc_losig = np.where(signew > 0.2*med_err)[0]
   print len(idxinc_losig),'len(idxinc_losig)'
   
   #remove points with too big errors
   idxinc_hisig = np.where(signew < 5.*med_err)[0]
   print len(idxinc_hisig),'len(idxinc_hisig)'
   
   i1 = np.intersect1d(idxinc_lo,idxinc_big)
   i2 = np.intersect1d(i1,idxinc_losig)
   i3 = np.intersect1d(i2,idxinc_hisig)
   
   nnew = np.shape(i3)[0]
   ynew = y[i3]
   tnew = t[i3]
   signew = sigy[i3]
   it = it + 1
  
  
  #not sure how to do this exclude
  #remove entire light curve if too big offset dispersion
  #if (tel == 'bok'):
  # bound = 5000
  #elif (tel == 'cfht'):
  # bound = 600
  
  
  
  
 elif (tel == 'spec'):
  print np.shape(y)
  print np.shape(sigy)

  idxkeep = np.where((np.abs(y) > 0) & (sigy > 0))[0]
  tnew = t[idxkeep]
  ynew = y[idxkeep]
  signew = sigy[idxkeep]
 
 
 return(tnew,ynew,signew)
 
 
 
 
 
 
 
###test routine on kates old light curve bs new light curve to make sure you get the same points out
##rm_005 a good test
#
#file1 = '/Users/ds207/Documents/standrews/sta/lightcurves/sdss_all_lc_oct2017/updated_quasar_lcs/cream_save/rm_005/spec_i_rm005_i.dat'
##file1 = '/Users/ds207/Google Drive/rm_campaign_cream/katelc_group_1febmerge_0/rm005_gimerge/rm005_i_spec.dat'
#file2 = '/Users/ds207/Documents/standrews/sta/fort/fortcode/lucky_results/rm_005/spec_i_rm005_i.dat'
##file1 = '/Users/ds207/Google Drive/rm_campaign_cream/katelc_group_1febmerge_0/rm005_gimerge/rm009_i_cfht_a2i.dat'
##file2 = '/Users/ds207/Documents/standrews/sta/fort/fortcode/lucky_results/rm_005/cfht_i_rm009_i_rs.dat'

#file1 = '/Users/ds207/Documents/standrews/sta/lightcurves/sdss_all_lc_oct2017/updated_quasar_lcs/cfht_i/rm608_i_rs.dat'

#
#
#
#dat1 = np.loadtxt(file1)
#dat2 = np.loadtxt(file2)
#n1 = np.shape(dat1[:,0])[0]
#n2 = np.shape(dat2[:,0])[0]
#
#dat1new = myrejout_kate(dat1[:,0],dat1[:,1],dat1[:,2],tel='phot',epochdelete=[])
#dat2new = myrejout_kate(dat2[:,0],dat2[:,1],dat2[:,2],tel='spec',epochdelete=[18,13])
#
##dat1new = myrejout_kate(dat1[:,0],dat1[:,1],dat1[:,2],tel='phot',epochdelete=[13])
##dat2new = myrejout_kate(dat2[:,0],dat2[:,1],dat2[:,2],tel='phot',epochdelete=[18,13])
#
#n1new = np.shape(dat1new[0])[0]
#n2new = np.shape(dat2new[0])[0]
#
#
#
#print 'old lc before after',n1,n1new
#print 'newlc before after',n2,n2new
#
#
#
#



 
 
 