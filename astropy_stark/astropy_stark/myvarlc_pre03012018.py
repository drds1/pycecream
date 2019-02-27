import scipy.optimize
import numpy as np
import myrandom as mr
import astropy.stats as ast
import matplotlib.pylab as plt
#jon trumps implementation of keiths intrinsic variability code
#
##function getrms, fluxes, errors
##
##;; Keith Horne's method to estimate intrinsic variability.
##;;
##;; V = sum( (xi-u)^2 gi^2 ) / sum(gi)
##;;   V = intrinsic variance
##;;   xi = flux measurement
##;;   u = mean flux
##;;   gi = 1 / (1+e^2/V)
##;;   e = flux error
##;;
##;; Iterate to solve for V.  Note that e should be a combination of
##;; the standard flux error plus the spectrophotometry error.
##
##  if stddev(fluxes) eq 0 or n_elements(fluxes) le 1 then return,[0,-1]
##
##  meanflux = median(fluxes)
##;;   errors = sqrt(errors^2 + (0.04*meanflux)^2)  ;add 0.04% specphoterror
##
##  varguess = (stddev(fluxes)^2 - median(errors)^2)>0.05
##
##  repeat begin
##     varold = varguess
##     gi = 1 / (1+errors^2/varold)
##     varnew = total( (fluxes-meanflux)^2 * gi^2 ) / total(gi)
##     varguess = varnew
##  endrep until abs(varold-varnew)/varnew le 1e-5 or varnew lt 1e-12
##
##  varerror = 2d / ( 2*total( (fluxes-meanflux)^2 / (varnew+errors^2)^3 ) $
##                            - total(1d/(varnew+errors^2)^2) )
##
##  return, [sqrt(varnew>0), sqrt(varerror>0)/sqrt(varnew>1e-4)]
##
##;  varerror2 = 2d * ( 2d*total((fluxes-meanflux)^2 * gi^3) $
##;                    / total((fluxes-meanflux)^2 * gi^2) * total(gi) - total(gi^2) )
##;
##;  return, [sqrt(varnew>0), sqrt(varnew>0)/sqrt(varerror2>1e-4)]
##
##end


#converted from Keiths fortran routine
#* =================================================
#	subroutine avgrmsx( n, dat, sig, good, avg, rms, sigavg, sigrms )
#* Maximum likelihood estimator of the mean and intrinsic rms
#* given independent samples with error bars
#* In:	n	i4 number of data
#*	dat(n)	r4 data
#*	sig(n)	r4 1-sigma error bars (<0 to skip)
#* Out:	good	r4 number of "good" data
#*	avg	r4 mean value
#*	sigavg	r4 1-sigma error bar
#*	rms	r4 excess rms
#*	sigrms	r4 1-sigma of excess rms
#* 2014 Aug Keith Horne @ St Andrews
#* 2014 Nov KDH @ StA - improve estimate of sigrms
#* 2016 Mar KDH @ StA - catch crash if no variance
#* 2016 Mar KDH @ StA - scale to <1/sig^2>=1. avoids vmin=0.


def myvarlc(fluxin,sigin):
 
 n = np.shape(fluxin)[0]

 #trap no data
 good = 0.
 avg = 0.
 rms = 0.
 sigavg = -1.
 sigrms = -1.
 if( n <= 0 ):
  print 'flux,sig arrays empty, no data',n
  return(good, avg, rms, sigavg, sigrms)

 
 #only use data with non zero error bar
 idinc = np.where(sigin > 0)[0]
 flux = fluxin[idinc]
 sig  = sigin[idinc]
 ng = np.shape(sig)[0]
 #trap no valid data
 if( ng <= 0 ):
  print '** ERROR in avgrmsx. n', n, ' ng', ng
  return(good, avg, rms, sigavg, sigrms)	
  
 #average positive error bars
 e1 = np.mean(sig)
 #KDH : V1 MAY VANISH IF E1 TOO SMALL
 v1 = e1 * e1
 
 
 #initial estimate
 x = e1/sig
 x2 = x*x
 sum1 = np.sum(flux*x2)
 sum2 = np.sum(x2)
 varguess = np.std(flux) 


 #optimal average and its variance
 avg = sum1 / sum2
 sigavg = e1 / np.sqrt( sum2 )
 v0 = e1 * ( sigavg / sum2 )

 #scale factor ( makes <1/sig^2>=1 )
 v1 = ng * e1 * ( e1 / sum2 )
 e1 = np.sqrt( v1 )
 v0 = v0 / v1
 
 
 #convergence threshold
 tiny = 2.e-5
 #lower limit for extra variance
 vmin = tiny * v0







#max-likelihood estimate of mean and extra varian
 nloop = 1000
 for loop in range(nloop):
  #stow for convergence test
  oldavg = avg
  oldrms = rms

  #data loop
  a = avg / e1
  
  d = flux/e1
  e = sig/e1
  #weight
  w = 1./(v0 + e*e)
  #"goodness" of the data 
  g = v0 * w
  x = g * (d - a)
  xx = x*x
  
  sum1 = np.sum(g*d)
  sumg = np.sum(g)
  sum  = np.sum(g*g)
  sum2 = np.sum(xx)
  sum3 = np.sum(g*xx)


  #"good" data points ( e.g. with sig < rms )
  good = sumg

  if( sumg < 0.0 ):
   print '** ERROR in AVGRMSX. NON-POSITIVE SUM(G)=', sumg
   print '** loop', loop, ' ndat', n, ' ngood', ng
   print '** sumg', sumg, ' sum1', sum1, ' sum', sum
   print '** tiny', tiny, ' v0', v0, ' vmin', vmin
   print '** dat(:)', flux[:]
   print '** sig(:)', sig[:]



  a = sum1 / sumg
  v0 = sum2 / sumg
  v0 = max( vmin, v0 )
  va = v0 / sumg
  #new avg and rms
  avg = a
  rms = np.sqrt( v0 )
  #error bars on avg and rms
  sigavg = np.sqrt( va )
  g = 2.0 * ( 2.0 * sum3 / sum2 * sumg - sum )
  
  
  
  sigrms = rms
  if( g > 0.0 ): 
   sigrms = rms / np.sqrt( g )

   
  #restore scaling
  avg = a * e1
  rms = rms * e1
  sigavg = sigavg * e1
  sigrms = sigrms * e1




  #KDH: CATCH SUM2=0 (NO VARIANCE)
  if( sum2 < 0.0):
   rms = 0.
   sigrms = -1.
   print 'problem CATCH SUM2=0 (NO VARIANCE)'
   return(good, avg, rms, sigavg, sigrms)

  #converge when test < 1
  if( loop > 1 ):
   safe = 0.9
   avg = oldavg * ( 1. - safe ) + safe * avg
   rms = oldrms * ( 1. - safe ) + safe * rms
   chiavg = ( avg - oldavg ) / sigavg
   chirms = ( rms - oldrms ) / sigrms
   test = max( abs( chiavg ), abs( chirms ) ) / tiny
   #report on last 5 iterations
   if( loop > nloop - 3 ):
	print 'Loop', loop, ' of', nloop, ' in AVGRMSX'
	print 'Ndat', n, ' Ngood', good, ' Neff', g
	print ' avg', avg, ' rms', rms
	print ' +/-', sigavg, ' +/-', sigrms
	print ' chiavg', chiavg, ' chirms', chirms, ' test', test
	if( test < 1. ):
	 print '** CONVERGED ** :))'

   #converged
   if( test < 1. ):
    return(good, avg, rms, sigavg, sigrms)

  #quit if v0 vanishes
  if( v0 <= 0. ): 
   return(good, avg, rms, sigavg, sigrms)

  #next loop


 #failed to converge
 print '** AVGRMSX FAILED TO CONVERGE :((((('
 return(good, avg, rms, sigavg, sigrms)






#python adaptation to jon trump iqrerror code
def iqrerror(flux,error):
 niter = 1000
 nv = np.shape(flux)[0]
 iqr = np.zeros(niter)

 for ee in range(niter):
  #for some reason jon trumps code doesn't use all the points
  #for the monte carlo resampling. Maybe to save time
  # i turned this off but you can put it back on if you want
  ###ind = long(randomu(seed,nv+1)*(nv+1))<nv
  
  flux1 = flux + np.random.randn(nv)*error
  #fluxsort = sort(flux1)

  x = np.percentile(flux1,[25,75])
  iqr[ee] = x[1]-x[0]
  #iqr[ee] = flux1[fluxsort[3*nv/4]] - flux1[fluxsort[nv/4]]
 return(np.std(0.74*iqr))









#linear fit
def func(x,a,b): 
 return(a+b*x)

#finction to fit linear trend to light curve and use above maximum likelihood estimator to 
#get the snr
#iterates linear fit over time using iq code to caclualte intrinsic errors
#iterativelt and refits linear line
#returns Keiths ML var and uncertainty
#snr in kates paper are ML/sigML from Keiths code i.e var[0]/var[1] (ML avgrms variance / uncertainty)
def linvarfit(time,flux,err,diag = 0,frac_intersect=0):
 
 iqr=0.0
 iter=0
 maxiter=10
 
 nv = np.shape(flux)[0]-1
 
 mederr = np.median(err)
 mederr2 = mederr*mederr
 iqr2 = iqr*iqr
 err2 = err*err
 toterr2 = err2+iqr2
 toterr  = np.sqrt(toterr2)
 x = np.median(np.abs(toterr2 - (err2 + iqr2))/flux) 
 c = ast.median_absolute_deviation(err2)
 c2 = c*c
 cc = c2/4/nv 
 
 #i dont think there is any need to iterate here 
 #but go ahead and see what happens
 while ((x < 1.e-4) and (iter <= maxiter)):
  
  
  #linear fit
  fit_coef=scipy.optimize.curve_fit(func,time,flux,sigma=toterr)
  npar = np.shape(fit_coef[1])[0]
  sig_coef = np.sqrt([[fit_coef[1][i,i]] for i in range(npar)])[:,0]
  fit_coef = fit_coef[0]  
  modelfit = fit_coef[0] + fit_coef[1]*time
  flux0 = flux - modelfit
  if (frac_intersect==1):
   idgood = np.where(np.abs(flux0) - toterr < 0)[0]
   ngood = np.shape(idgood)[0]
   npcgood = 1.*ngood/nv
  else:
   npcgood = 1.
  #make sure linear fit working properly
  if (diag==1):
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   ax1.errorbar(time,flux,toterr,ls='',capsize=2,label='toterr')
   ax1.errorbar(time,flux,err,ls='',capsize=1,label='errorbars')
   ax1.plot(time,modelfit,color='k',label='fit')
   ax1.set_title('iteration '+np.str(iter))
   plt.legend()
   plt.show()

  #ML variance estimator
  var_keith = myvarlc(flux0,err)
  var_trump = [var_keith[2],var_keith[2]/var_keith[4]]
  
  #IQR range of 25-75% of cumulative distribution,
  # with correction (x0.74) to equal sigma for a Gaussian
  x = np.percentile(flux0,[25,75])
  
  vsort = np.sort(flux0)
  iqr_obs = 0.74 * (x[1] - x[0])
  iqrerr_obs = iqrerror(flux0,err)
  iqr = np.sqrt(np.max(iqr_obs**2 - mederr2,0))
  iqr2 = iqr*iqr
 
  iqrerr = 1/iqr * np.sqrt( (iqr_obs*iqrerr_obs)**2 + cc)
  toterr2 = err2+iqr2
  toterr  = np.sqrt(toterr2)
  x = np.median(np.abs(toterr2 - (err2 + iqr2))/flux) 

#the below are another way of getting the snr still coded in c need to trnslate to 
#python one below half done dont need final one for routine

#AD: absolute average deviation (e.g. vanden Berk et al 2004)
#   aad = getaad(flux0,err)
#
  #print 'median error', np.median(err)
  print 'keith ml var and var/var_error',var_trump
  #print 'iqr and iqrerr', iqr,iqrerr
  
  iter = iter + 1

 if (frac_intersect==1):
  return(var_trump,npcgood)
 else:
  return(var_trump)










