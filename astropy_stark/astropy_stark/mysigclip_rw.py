#David Starkey UIUC 08/02/2018
#code to use sigma clipping algorithm to reject outliers from a randomw walk model fit
#import myfitrw_current as mfw
import myfitrw_simp as mfw
#import astropy.stats as apys
from pylab import *
from mylcgen import *

def mscrw(x,y,sig=[],niteration=10,sigmaclip = 5):

 xin = 1.*x
 nx = np.shape(x)[0]
 idxinclude = np.arange(nx)
 tdat = x[idxinclude]
 ydat = y[idxinclude]
 
 #if no error bars just treat everything equally
 if (sig == []):
  sigdat = np.ones(nx)
 else:
  sigdat = sig[idxinclude]
 
 #this bit tells the code to use the errors on the model rather than the data errorbars if
 #your data has no error bars
 if (sig == []):
  nits = 100
 else:
  nits = 1
  
  
 for i in range(niteration):
  
  nxold = np.shape(sigdat)[0]
  op = mfw.fitrw([tdat],[ydat],[sigdat],floin=-1,fhiin=-1,ploton=1,dtresin=-1,nits = nits)
  ymod = op[-2]
  ymodsd = op[-1]
  yin = np.abs(ydat - ymod)
  
  if (sig == []):
   idinc = np.where(yin < sigmaclip*ymodsd)
  else:
   idinc = np.where(yin < sigmaclip*sigdat)
   
  #clf()
  #
  #
  #errorbar(tdat,ydat,sigdat,ls=None)
  #plot(tdat,ymod)
  #title('sdfsdfs')
  #print idinc
  for i in range(np.shape(yin)[0]):
   print ydat[i],ymod[i],ymodsd[i]
  show()
  
  #a = apys.sigma_clip(yin, sigma=sigmaclip, sigma_lower=None, sigma_upper=None, iters=1, axis=None, copy=True)
  #idinc = ~a.mask
  tdat = tdat[idinc]
  ydat = ydat[idinc]
  sigdat = sigdat[idinc]
  
  nxnew = np.shape(sigdat)[0]
 
  print 'sigclip report'
  print 'nold',nxold
  print 'nnew',nxnew
  print 'rejected',nxold-nxnew,'points'
  if (nxold == nxnew):
   break
 return(tdat,ydat,sigdat)
 
 
 
 
 

 
 
 
 
 
 
##!!!!!!!!!!!!!!!!!!! test the code using this program
#
#





#generate test light curve and add noise
cadence = 3.0
timeall = []
sigall  = []
yall    = []
datpre      = mylcgen(tlo=0,thi=100,dt=cadence,iseed=132423)
npre     = np.shape(datpre[:,0])[0]
sig = np.std(datpre[:,1])/10
#add noise
sigpre = np.ones(npre)*sig
datpre[:,1] = np.random.randn(npre)*sig + datpre[:,1]
#add bad data to test sigma clipping
datpre[10,1] = datpre[10,1] + 100*sig
datpre[5,1] = datpre[5,1] - 40*sig

timeall = datpre[:,0] + np.random.randn(npre)*0.1#npre,datpre[:,0],0.1)
yall = datpre[:,1]
sigall = sigpre





#if your data has error bars, put these in as the sigma argument
#a = mscrw(timeall,yall,sig = sigall,niteration=3,sigmaclip = 2)

#if your data has no error bars use the error snakes on the fitted model in the sigma clip
a = mscrw(timeall,yall,niteration=1,sigmaclip = 5)




#make a seperate plot of old and new data
clf()
scatter(timeall,yall)
scatter(a[0],a[1])

show()
datmean  = np.std(datpre[:,1])
snow =  np.ones(npre)/10*datmean 
#dat = #mrs.myresample(dir='',fname=[''],dtave=2.0,sampmin=0.8,sampcode=3,datin=np.array((datpre[:,0],datpre[:,1],snow)).T)
ndat = np.shape(dat[:,0])[0]
sig = dat[:,2]
for i in range(ndat):
 dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
sigall.append( sig )
yall.append( dat[:,1] +50 )
timeall.append( dat[:,0] )



a = fitrw(timeall,yall,sigall,ploton=1) 



#
#     
 