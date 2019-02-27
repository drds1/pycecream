import numpy as np
import matplotlib.pylab as plt



#iteratively fit polynomials up to order, using BIC,AIC and reduced chi_squared to determine 
#the optimum number of parameters.
def fit_search(x,y,maxorder=8,xgrid=[]):
 
 
 
 yg_med   = []
 yg_lo    = []
 yg_hi    = []
 cov      = []
 cisq     = []
 cisq_red = []
 bic      = []
 aic      = []
 
 sig = 0.*x + 1.
 for i in range(1,maxorder):
  a = fit(x,y,sig,order=i,xgrid=xgrid)
  yg_med.append(a[0])
  yg_lo.append(a[1])
  yg_hi.append(a[2])
  cov.append(a[3])
  cisq.append(a[4])
  cisq_red.append(a[5])
  bic.append(a[6])
  aic.append(a[7])

 #sort by increasing fit metrics
 id_aic     = np.argsort(aic)
 id_bic     = np.argsort(bic)
 id_cisqred = np.argsort(np.abs(np.array(cisq_red) - 1.0))#np.argsort(cisq_red)
           
 print('best cisqred fit order:',id_cisqred[0])
 print('best aic fit order:',id_aic[0])
 print('best bic fit order:',id_bic[0])
 return(id_cisqred[0]+1,id_aic[0]+1,id_bic[0]+1)














#fit polynomial of a given order conf=0.05 for 95 % , 0.3173 for 1 sigma
def fit(x,y,sig=0,order=3,xgrid=[],confidence=0.3173,nits=20000,figure_title=''):

 #evaluate model on arbitrary time grid
 if (xgrid != []):
  xg = np.array(xgrid)
 else:
  xg = np.array(x)
 nxg = np.shape(xg)[0]

 #cannot fit trend with fewer points than polynomial order
 nx = np.shape(x)[0]
 oi = min(nx-4,order)
 if (oi < 0): 
  print('fit is degenerate, just using straight line')
  try:
   yg_med = np.ones(nxg)*y[0]
  except:
   yg_med = np.zeros(nxg)
  yg_lo  = yg_med*1
  yg_hi  = yg_med*1
  cov    = np.zeros((1,1))
  disq   = 0
  cisq = 0
  cisq_red = 0
  bic = 0
  aic = 0
  rmd = 0
  return(yg_med,yg_lo,yg_hi,cov,cisq,cisq_red,bic,aic,rmd)

 #if sig = [] then dont know error bars assign equal weight
 if (type(sig) == np.ndarray):
  s = 1.*sig
  noerrors = 0
 else:
  s = np.ones(nx)
  noerrors = 1
 
 
 

 #compute and subtract pivot point (get better fit with reduced correlated parameter uncertainties)
 xp = np.sum(x/s**2)/np.sum(1./s**2)
 
 
 coefs,cov = np.polyfit(x - xp,y,w=1./s,deg = oi,full=False,cov=True)
 
 if (noerrors == 1):
  
  ymod_itp = np.sum(np.array([coefs[-1-ip]*(x-xp)**ip for ip in range(oi+1)]),axis=0)
  
  print(np.shape(ymod_itp))
  #input()
  
  var =   np.sum((ymod_itp - y)**2)/nx
  
  #print('var',var)
  #cov = cov*var
 
 

 #monte carlo sample error snake
 c_mvn  = np.random.multivariate_normal(coefs,cov,size=nits)
 yg = np.zeros((nxg,nits))
 
 
 
 xgxp_grid = np.array([(xg-xp)**ip for ip in range(oi+1)])

 
 
 for iteration in range(nits):
  cnow = c_mvn[iteration,:]
  yg_temp = np.dot(xgxp_grid.T,cnow[::-1]) 
  #yg_temp = np.sum( np.array([cnow[order-ip]*(xg-xp)**ip for ip in range(order+1)] ), axis = 0)  

  yg[:,iteration] = yg_temp


 #compute final model and uncertainties
 #ygrid = np.mean(yg,axis=1)
 
 pc = [confidence/2 *100, 50. ,(1.-confidence/2)*100 ]
 yg_med_conf = np.percentile(yg,pc,axis=1).T
 
 yg_lo  = yg_med_conf[:,0]
 yg_med = yg_med_conf[:,1]
 yg_hi  = yg_med_conf[:,2]
   
 

 #interpolate model onto data
 yg_itp  = np.interp(x,xg,yg_med)
 
 
 
 
 
 

 
 #modify uncertainty to suit confidence
 ymed    = np.median(yg_med)
 diffs   = y - yg_itp
 id_hi     = np.where(diffs > 0)[0]
 diff_hi   = np.sort(diffs[id_hi])
 nhi       = np.shape(diff_hi)[0] 
 try:
  frac_hi   = diff_hi[np.int(nhi*(1.-confidence))] / np.median(yg_hi - yg_med)
 except:
  frac_hi = 1
 
 yg_hi_mod = (yg_hi - yg_med)*frac_hi + yg_med
 
 
 
 # print(nhi,'herere',order)
 # print(coefs)
 #
 # for it in range(np.shape(yg_med)[0]):
 #  print(xg[it],yg_lo[it],yg_med[it],yg_hi[it])
 # print()
 # for it in range(np.shape(yg_itp)[0]):
 #  print(x[it],y[it],yg_itp[it])  




 
 
 
 #for i in range(np.shape(yg_lo)[0]):
 # print(i,xg[i],yg_lo[i],yg_med[i],yg_hi[i])
 #print('cov')
 #print(cov)
 #print('pars')
 #print(coefs)
 #print('mvg')
 #print(c_mvn)
 #
 #print(frac_hi,diff_hi[np.int(nhi*(1.-confidence))],np.median(yg_hi - yg_med))
 #print()
 
 
 
 id_lo     = np.where(diffs < 0)[0]
 diff_lo   = np.sort(np.abs(diffs[id_lo]))
 nlo       = np.shape(diff_lo)[0]
 try:
  frac_lo   = diff_lo[np.int(nlo*(1.-confidence))] / np.median(yg_med - yg_lo)
 except:
  frac_lo = 1.0
 yg_lo_mod = yg_med - (yg_med-yg_lo)*frac_lo
 
 #print(frac_lo,diff_lo[np.int(nlo*(confidence))],np.median(yg_med - yg_lo))
 #print() 
 #input()
  
  
  
  
 #absdiff_sort = np.sort(np.abs(diff))#np.where(np.abs(diff)
 #frac = absdiff_sort[np.int( nxg*(1.-confidence) )]/absdiff_sort[np.int(nxg/2)]
 
  
  
 #interpolate model uncertainties onto data
 hi_itp  = np.interp(x,xg,yg_hi_mod)
 lo_itp  = np.interp(x,xg,yg_lo_mod)
 sig_itp = (hi_itp - lo_itp)/2
 
 
 #compute model evaluation stats (big numbers for lots of parameters or poor fits)
 #reduced chisquared
 cisq = np.sum(((yg_itp - y)/sig_itp)**2)
 ndof = nx - oi + 1
 cisq_red = cisq/ndof
 
 #autocoorelation
 rmd = np.corrcoef(yg_itp,y)[0,1]
 
 #Bayesian information criterion 
 bic = cisq + (oi+1)*np.log(nx)
 
 
 #Akaike information criterion
 aic = cisq + 2*(oi+1)
 
 
 #results summary
 print('fit results for poilynomial order: ',oi)
 print('')
 print('Smaller numbers better')
 print('reduced chisquare: ',cisq_red)
 print('aic: ',aic)
 print('bic: ',bic)
 print('chisquare: ',cisq)
 print('correlation coefficient: ',rmd)
 print()
 print()
 

 #fig = plt.figure()
 #ax1 = fig.add_subplot(111)
 #ax1.scatter(x,y)
 #ax1.plot(xg,yg_med,label='model')
 #ax1.plot(x,yg_itp,label='interpolated')
 #ax1.fill_between(xg,yg_lo,yg_hi,alpha=0.3,label='old uncertainties')
 #ax1.fill_between(xg,yg_lo_mod,yg_hi_mod,alpha=0.3,label='new uncertainties')
 #plt.legend()
 #plt.show()
 
 

 if (figure_title != ''):
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.scatter(x,y)
  ax1.plot(xg,yg_med,label='model')
  ax1.plot(x,yg_itp,label='interpolated')
  ax1.fill_between(xg,yg_lo,yg_hi,alpha=0.3,label='old uncertainties')
  ax1.fill_between(xg,yg_lo_mod,yg_hi_mod,alpha=0.3,label='new uncertainties')
  plt.legend()
  if (figure_title == 'show'):
   plt.show()
  else:
   plt.savefig(figure_title)
  
 return(yg_med,yg_lo,yg_hi,cov,cisq,cisq_red,bic,aic,rmd)
 
 
 
 
 
 
 
 
 
 
##test with fake data
#
#import numpy as np
#import matplotlib.pylab as plt
#import glob
#import pandas as pd
#import datetime
#import statsmodels.api as sm
#
#combine = 0
#frequency = 1./30. #fake 6 month signal
#f2 = 1./60
#tlo = 0.0
#thi = 365
#tref = 0.0#thi/2
#dt = 1.0
#lab = ''
#amp1 = 5.0
#amp2 = 10.0
#amps = [amp1,amp2]
#freqs = [frequency,f2]
#poly = [0,0,5.e-5,0.0]
#tfclo = thi
#tfchi = thi + 365
##generate synthetic time series for test
#t = np.arange(tlo,thi,dt)
#x = 5*np.sin(2 * np.pi * frequency * t) + 10*np.sin(2 * np.pi * f2/2 * t)
##add trend final -10.0 indicates amplitude at the end of the time sequence
#grad = 1.0
#y1   = 10.0
#noiseamp = 0.5
#
#pc_forecast = 0.4 #forecast ahead 40% the original length of the time series
#
##specify the confidence interval with the alpha argument 
##e.g alpha = 0.05 is the 95pc confidence interval (2 sigma)
##alpha = 0.32 is the 1 sigma confidence interval (68pc)
#alpha_conf = 0.32#0.05
#
#
#
#def trend(t,poly,tref,freqs,amps):
# nt = np.shape(t)[0]
# npoly = len(poly)
# namps = len(amps)
# xtot = []
# for i in range(nt):
#  trend = np.sum([poly[ip]*(t[i]-tref)**ip for ip in range(npoly)] )
#  seasonal = np.sum([amps[ip]*np.sin(2*np.pi*freqs[ip]*t[i]) for ip in range(namps)] )
#  xtot.append(trend + seasonal)
#
# return(np.array(xtot))
# 
# 
#

#signal = trend(t,poly,tref,freqs,amps) 
##grad*(x-t[-1]) + y1
#nt = np.shape(t)[0]
#noise = np.random.randn(nt)*noiseamp
#signal = signal + noise
#
#t = np.arange(100)
#x = np.random.randn(100)*1.0 + 5.0 #+ 2.3*t + 3*t**2
#
#
#tgrid = np.arange(-20,120,0.1)
#fit_search(t,x,maxorder=8,xgrid=tgrid)

##signal = trend(t,poly,tref,freqs,amps) 
###grad*(x-t[-1]) + y1
##nt = np.shape(t)[0]
##noise = np.random.randn(nt)*noiseamp
##signal = signal + noise
##
##t = np.arange(100)
##x = np.random.randn(100)*1.0 + 5.0 #+ 2.3*t + 3*t**2
##
##
##tgrid = np.arange(-20,120,0.1)
##fit_search(t,x,maxorder=8,xgrid=tgrid)
