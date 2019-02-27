import numpy as np
import scipy
import scipy.stats


#randomly sample nsamp from a distribution in bins of width binwid
#if binwidin -ve the absolute value specifies number of bins rather than the width
#if specbw 1, then alter bin widths to ensure equal samples in each bin
def myransamp(dist,nsamp,binwidin,dloin='',dhiin='',specbw = 1):
 
 if (specbw == 1):
  print 'myransamp.py adaptive bin size ON'
 
 idxout = []
 binout = []
 if (dloin == ''):
  dlo = np.min(dist)
 else:
  dlo = dloin
  
 if (dhiin == ''):
  dhi = np.max(dist)
 else:
  dhi = dhiin

 blosave = []
 bwsave = []
#one option specifies binwidth the other specifies bin number with option to equalize number of objects in each bin 
 if (binwidin < 0):
  nbin = np.abs(binwidin)
  for i in range(nbin):
   if (specbw == 1):
    dsort = np.sort(dist)
    ndist = np.shape(dist)[0]
    nb = int(np.floor(ndist/nbin)) 
    blonow = dsort[i*nb]
    blosave.append(blonow)
    if (i < nbin - 1):
     blonext = dsort[(i+1)*nb]
     bwsave.append(blonext - blonow)
    else:
     bwsave.append(dhi - blonow)   
   else:
    binwid = (dhi - dlo)/(nbin - 1)
    blosave.append(dlo + i*binwid)   
    bwsave.append(binwid)
    
 else:
  binwid = binwidin
  nbin = np.floor(dhi-dlo)/binwid + 1
  for i in range(nbin):
   blosave.append(dlo + binwid*i)
   bwsave.append(binwid)
    
  
 nperbin = nsamp /nbin 
 
 for i in range(nbin):
  blonow = blosave[i]#dlo + i*binwidth
  bhinow = blonow + bwsave[i]#blonow + binwidth
  idxsamp = np.where((dist > blonow) & (dist < bhinow))[0]
  
  nss = np.shape(idxsamp)[0]
  nnow = min(nss,nperbin)
  if (nnow < nperbin):
   print 'bins too fine for number of requested samples per bin, picking smaller sample size',nnow,nperbin
  idbin = np.random.choice(idxsamp,size=nnow,replace=False)
  idxout.append(idbin)
  binout.append(blonow)
  
 return(blosave,bwsave,idxout)
  
  
  
  
  
  
  
  
  
#this code will bin the data into a 2D histogram of approximately equal 
#bin sizes, then randomly sample nsamp  
#minbin roughly how many samples you want per bin
#datout should be a list of approximately uniformly sampled dat

#if nbinin > 0 then this specifies the number of linearly spaced bins in each dimension
#from the lowest to highest data point

#dlimin = [[dlo,dhi],[dlo,dhi].....] for each dimension limits for each bin
#... else if -1 use minimum and maximum point 
#...else use dlminin*nsd about mean either side

def myransampnd(dat,nsamp,minbin=10,nbinin = -1,dlminin=-1):
 ndim = np.shape(dat[0,:])[0]
 ndat = np.shape(dat[:,0])[0]
 datout = []
 if (isinstance(dlminin,list)):
  dlo = [dlminin[i][0] for i in range(ndim)]
  dhi = [dlminin[i][1] for i in range(ndim)]
 elif (dlminin == -1):
  dlo = [np.min(dat[:,i]) for i in range(ndim)] 
  dhi = [np.max(dat[:,i]) for i in range(ndim)] 
 elif (dlminin > 0):
  dlo = [np.mean(dat[:,i]) - dlminin*np.std(dat[:,i]) for i in range(ndim)] 
  dhi = [np.mean(dat[:,i]) + dlminin*np.std(dat[:,i]) for i in range(ndim)] 

 #need fewer bins than samples I think this should give about minbin samples per bin
 if (nbinin == -1):
  nbin = nsamp**(1./ndim)/minbin
 else:
  nbin = nbinin
  
 
 #decide bin sizes
 bins = []
 idxselect = []
 for i in range(ndim):
  if (nbin == -1):
   ddim = dat[:,i]
   dsort = np.sort(ddim)
   dbin = dsort[::minbin]
  else:
   dbin = np.linspace(dlo[i],dhi[i],nbin)
   
  bins.append(dbin) 
  #idxselect = 
  
 bins = np.array(bins)
 print np.shape(bins)
 print np.shape(dat)
 dhist = scipy.stats.binned_statistic_dd(dat,np.ones(ndat),bins=bins)
 idbin = dhist[-1]
 idbinu = np.unique(idbin)
 
 
 print np.shape(idbinu)
 print np.shape(idbin)
 for idnow in idbinu:
  try:
   idxu = np.where(idbin == idnow)[0]
   nnow = 1.0#min(np.shape(idxu)[0],minbin)
   idpick = np.random.choice(idxu,size=nnow,replace=False)
   datout.append(dat[idpick,:])
  except:
   continue
 
 return(datout) 
  
 #np.histogramdd(dat,bins=bins)



#myransamp(dist,nsamp,binwidin,dloin='',dhiin='',specbw = 1):
