import numpy as np
import scipy
import scipy.optimize


#ingest set of data find median and uncertainties
#input dat array
# op dmed median of dat, dlo,dup lower and upper sig confidence limits (i.e for 1 sigma, set sig = 0.68)
# disg = (dup-dlo)/2


def stat_uncert(dat,sig=0.68):
 nd = len(dat)#.shape()
 
 dat_sort=np.sort(dat)
 dmed = np.median(dat_sort)
 
 f = (1.-sig/2)
 nup = int(np.ceil((1.-f)*nd))
 nlo = int(np.ceil(f*nd))
 
 
 dlo = dat_sort[nlo]
 dup = dat_sort[nup]
 dsig = np.abs((dup-dlo)/2)
 
 return(dmed,dsig,dlo,dup)
 
 



#return uncertainty in theta give cos(theta) and sig_(costheta)



def sigtheta(costheta, sd_ct):

 ct_mean = np.mean(costheta)
 
 a = np.arccos(ct_mean)
 degmean = 180./np.pi*a
 sdmean =  1./np.sin(a) * sd_ct * 180/np.pi
 
 return(degmean,sdmean)
 
 





# fit polynomial to data by orthogonalizing parameters

#n order of poly
#x,y,sig array of data

def polyfit(x,y,sig,n):
 
 sum1oversig = np.sum(1./sig)
 
 xorth = np.sum(x**n / sig**n) / sum1oversig  
 
 xsuborth = x - xorth
 
 
 
 #fit a nth order poly to xsuborth data
 
 if (n == 1):
  def func(x,a,b): 
   return(a+b*x)
 elif (n==2):
  def func(x,a,b,c):
   return(a+b*x+c*x*x)
 elif (n==3):
  def func(x,a,b,c,d):
   return(a + b*x+c*x*x + d*x*x*x)
 elif (n==4):
  def func(x,a,b,c,d):
   return(a + b*x+c*x*x + d*x*x*x + e*x*x*x*x*x)
 else:
  raise ValueError('This function only works with polynomials up to n=4')
 
 
 
 fit_coef=scipy.optimize.curve_fit(func,np.array(x)-xorth,np.array(y),sigma=np.array(sig))
 sig_coef = np.array(())
 
 for i in range(n+1):
  sig_coef=np.append(sig_coef,np.sqrt([fit_coef[1][i,i]]))
  
 fit_coef = fit_coef[0]
 return(fit_coef,sig_coef,xorth)
  
 
 
 
 

#bayesian information criterion data x, model y,sig, k free params
def bic(x,y,sig,k):
 
 xsuby = x-y
 xsuby2 = xsuby*xsuby
 sig2 = sig*sig
 n=np.shape(x)[0]
 
 chisq = np.sum(xsuby2/sig2)
 
 bic = chisq + k*np.log(n) + n*np.log(2*np.pi) + np.sum(np.log(sig2))
 
 return(bic)
 

 
 
 
#convert uncertainty in  log or cos to real parameter and back again
#convert mean parameter (aveold) and uncertainty sdold and mean both from (1) and to (-1) log par and uncertainty, 
#and from (2) and to (-2) cosine par with uncertainty



def myconvlogsig(iversion,aveold,sdold):


 rln10 = np.log(10.)
 deg2rad = np.pi/180


#if going from log to real
 if (iversion == 1):
  a = 10**aveold
  avenew = a*(1. + 0.5*sdold*sdold*rln10*rln10)
  sdnew  = rln10*a * sdold



 #!! real to log
 elif (iversion == -1):

  a = alog10(aveold)
  avenew = a - 0.5*sdold*sdold/(rln10*aveold*aveold)
  sdnew  = 1./(rln10*aveold) * sdold







#!! cos to real
 elif (iversion ==2):
  print 'mystats doing fuck all'






# real to cos
 elif (iversion == -2):
  print 'mystats doing fuck all'
 
 

 else:
  print 'mystats.f90: iversion',iversion,' is not an option!!'
  exit
  
 
 a = cos(aveold*deg2rad)

 avenew = a*(1. - 0.5*sdold*sdold)
 sdnew  = sin(aveold*deg2rad)*sdold






 return(avenew,sdnew)
 
 
 
 
#function to return the fwhm of a distribution if funny shape just use ends 
def fwhm(t,x):
 xmin = np.min(x)
 xmax = np.max(x) - xmin
 xnew = (x - xmin)/xmax
 idx  = np.where(xnew > 0.5)[0]
 tlo  = t[idx[0]]
 thi  = t[idx[-1]]
 return(thi - tlo)
 
def fwhmall(t,x):
 xmin = np.min(x)
 xmax = np.max(x) - xmin
 xnew = (x - xmin)/xmax
 idx  = np.where(xnew > 0.5)[0]
 tlo  = t[idx[0]]
 thi  = t[idx[-1]]
 return(tlo,thi)
 
 
 
 
 