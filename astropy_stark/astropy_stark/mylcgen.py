#generate a light curve with a bbroken power law power spectrum according to the prescription of timmer and kronig 1995




import numpy as np
#import myrandom as mr
import matplotlib.pylab as plt



def mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=0,thi=100,dt=0.125,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0):
 '''
 Generate random walk light curve. Can also customize the powr spectrum slope using a and b arguments
 :param datfile:
 :param p0:
 :param f0:
 :param a:
 :param b:
 :param tlo:
 :param thi:
 :param dt:
 :param ploton:
 :param iseed:
 :param meannorm:
 :param sdnorm:
 :return:
 '''
 #!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!
 #Input parameters
 
 #use clock to get different random light curve each call if -1, else iseed specifies light curve
 if (iseed > 0):
  np.random.seed(int(iseed))
 else:
  np.random.seed()
 
 flo = 0.5/(thi-tlo)
 fhi = 1./dt
 #datfile ='mylcggen_op.dat'
 #
 #p0 = 1.0      # scale factor 
 #f0 = 0.1      #break frequency cycles per day
 #a  = -2.0
 #b  = -2.0     # second and first slopes of the broken power law
 #tlo = 0.0
 #thi = 100.0
 #dt  = 0.1
 #
 #flo = 0.01
 #fhi = 4.0     # low and hi fourier frequency limits 
 df = 1.*flo
 #
 #ploton = 1
 
 
 #!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!
 
 
 
 
 
 
 
 
 
 
 
 
 
 #!!!!!!!!!#!!!!!!!!!#!!! Make and save the light curve !!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!
 
 time = np.arange(tlo,thi+dt,dt)
 nt = np.shape(time)[0]
 x    = np.zeros(nt)
 dat  = np.zeros((nt,2))
 
 nf = int(np.ceil((fhi - flo)/df + 1))
 freq = np.arange(flo,fhi+df,df)
 w    = 2*np.pi*freq
 
 a   = np.sqrt( p0*(freq/f0)**a / (1 + (freq/f0))**(a-b))
 sk  = np.random.randn(nf)*a#mr.normdis(nf,0,1)*a
 ck  = np.random.randn(nf)*a#mr.normdis(nf,0,1)*a
 
 dps = np.zeros((nf,3))
 dps[:,0] = freq
 dps[:,1] = sk
 dps[:,2] = ck
 
 #print 'mylcgen iseed...',iseed
 for i in range(nt):
  tnow = time[i]
  x[i] = np.sum(sk*np.sin(w*tnow) + ck*np.cos(w*tnow))
  #print nt, i,'here'
 dat[:,0] = time
 dat[:,1] = x
 
 if (len(datfile) > 0):
  np.savetxt(datfile,dat)
  np.savetxt('fft_'+datfile,dps)
 #!!!!!!!!!#!!!!!!!!!#!!!#!!!!!!!!!#!!!!!!!!!#!!!#!!!!!!!!!#!!!!!!!!!#!!!#!!!!!!!!!#!!!!!!!!!#!!!
 
 
 
 #normalise if option turned on
 if (sdnorm > 0):
  datsd = np.std(dat[:,1])
  datmean = np.mean(dat[:,1])
  dat[:,1] = (dat[:,1] - datmean)*sdnorm/datsd + datmean
 
 if (meannorm >= 0):
  dat[:,1] = dat[:,1] - dat[:,1].mean() + meannorm
 
 
 #!!!!!!!!!#!!!!!!!!!#!!! If ploton make power spectrum and plot#!!!!!!!!!#!!!!!!!!!#!!!#!!!!!!!!!#!!!!!!!!!#!!!
 if (ploton == 1):
  fig = plt.figure()
  ax1=fig.add_subplot(211)
  ax1.plot(time,x)
  
  ax2=fig.add_subplot(212)
  ax2.plot(freq,sk**2+ck**2,ls='',marker='o')
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  plt.show()
 
 return(dat)





