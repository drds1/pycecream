import numpy as np




#combine signals with specific uncertainties and weights r
#should all be on same time axis prior to inserting in here
#
#signal is list of arrays each array is N X 3 with 1,2 and 3 columns are the time, flux 
#and uncertainty axis
#r2 are a list of weights from 0 to 1 on how much to regard each series 
#can use correlation coefficients for this
#first axis indicates time, 2nd indicates separate signals

def comb_signal(x,noise,r):

 n2 = noise**2

 top = np.sum(r*x/n2,axis=1)
 bot = np.sum(r/n2,axis=1)
 signal_out   = top/bot
 signal_noise = np.sqrt(1./bot)
 
 
 #will crash if zeros used for uncertainties
 if 0 in noise:
  raise Exception('cannot have 0 value for noise array in vortexa_combine_signal.py')
  
 
 #nx,ny = np.shape(x)
 #for i in range(nx):
 # print(i,noise[i,:])
 #print('mean',np.mean(x,axis=1))
 #print('rms',np.std(x,axis=1))
 #print('top',top)
 #print('bot',bot)  
 #input()

 return(signal_out,signal_noise)
 
 
  
 
 
 
 
 
 
 
 
 
 
 
 
 
##test combine signal code\
#import mylcgen as mlc
#import matplotlib.pylab as plt
# 
# 
# 
# 
#anoise = 0.1
#bnoise = 0.1
#tlo = 0.0
#thi = 100.0
#dt = 1.0
#
#a = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=tlo,thi=thi,dt=dt,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
#na = np.shape(a[:,0])[0]
#asd = np.std(a[:,1])
#a[:,1] = np.random.randn(na)*anoise*asd + a[:,1]
#
#
#b = mlc.mylcgen(datfile='',p0=1.0,f0=0.1,a=-2,b=-2,tlo=tlo,thi=thi,dt=dt,ploton=0,iseed=-1,meannorm = -1., sdnorm = -1.0)
#nb = np.shape(a[:,0])[0]
#bsd = np.std(a[:,1])
#b[:,1] = np.random.randn(nb)*bnoise*bsd + b[:,1]
#
#
#
#
#
#noise = np.zeros((na,2))
#noise[:,0] = anoise*asd
#noise[:,1] = bnoise*bsd
#
#signal = np.zeros((na,2))
#signal[:,0] = a[:,1]
#signal[:,1] = b[:,1]
#
#
#r  = [1.0,1.0]
#cb = comb_signal(signal,noise,r)
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.fill_between(a[:,0],signal[:,0]-noise[:,0],signal[:,0]+noise[:,0],label='signal 1',alpha=0.6)
#ax1.fill_between(a[:,0],signal[:,1]-noise[:,1],signal[:,1]+noise[:,1],label='signal 2',alpha=0.6)
#ax1.fill_between(a[-50:,0],cb[0][-50:]-cb[1][-50:],cb[0][-50:]+cb[1][-50:],label='combined signal',alpha=0.6)
#ax1.legend()
#plt.show()
#
#






