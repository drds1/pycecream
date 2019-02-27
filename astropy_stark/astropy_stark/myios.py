#this should be here already

def ioslinfit(dat,sig,modin,cspec=0.0,sspec=0.0, specgen = 1.0, rms_pri = 0.0):
 
 #cspec=1 fix offset to mean of data
 #.......=2 ............. median
 #sspec = 1 ...force s to rms of data
 #sspec = 2 force sspec to be 1
 
 
 #if (cspec == 1):
  
 if (specgen == 1):
  mod = modin - np.mean(modin)
 else:
  mod = 1*modin
 
 n=100
 s = 1.0
 c = 1.0
 for i in range(n):
  
  top = np.sum( (dat-s*mod)/sig**2 )
  bot = np.sum( 1./sig**2 )
  c = top/bot
  #print s, c, 'myios.iosinterp'
  sigc = 1/np.sqrt(bot)
  if (cspec==2):
   c = np.median(dat)
  elif (cspec ==1):
   c = np.mean(dat)
  elif (cspec < 0):
   c = np.abs(cspec)
  
  top = np.sum( (dat-c)*mod/sig**2 )
  bot = np.sum( mod**2/sig**2 )
  s = top/bot
  if (sspec==2):
   s = np.abs(s)
  elif (sspec ==1):
   s = np.std(dat)/np.std(mod)
  elif (sspec < 0):
   s = np.abs(sspec)
  
  sigs = 1/np.sqrt(bot)
  
  if (rms_pri == 1):
   modnew = s*mod + c
   s = s*np.std(dat)/np.std(modnew)
  
  
  #print i, s, c,'asdkjsadja'
  
 return(s,sigs,c,sigc) 
  



### python implementation of iterated optimal scaling
#inputs: t: time array
#        x: data array
#     xsig: error bar array
#        w: angular frequency array
#  sk1 ck1: best guess for inut sine and cosine 

#outputs: p[NW,2] sine and cosine amplitudes
#      sigp[NW,2] uncertainties on sine and cosine amplitudes  


import numpy as np
from pylab import *


def ios(t,x,xsig,w,sk1,ck1):
    overflow=1.e10

    ck=1.*ck1
    sk=1.*sk1
    nt=x.shape[0]
    nw=w.shape[0]
    xmean=np.mean(x)
    sigsq=1.*xsig*xsig
    p=np.zeros((nw,2))
    sigp=np.zeros((nw,2))
    nits=100
##### subtract the mean and return this as a constant offset
    xdup=x-xmean
    cw=np.zeros((nt,nw))
    sw=np.zeros((nt,nw))
    for it in range(nt):
        cw[it,:]=np.cos(t[it]*w[:])
        sw[it,:]=np.sin(t[it]*w[:])
    
   
## sin parameter    
    for iw in range(nw):
        for iteration in range(nits):
            top=np.sum((xdup[:]-ck*cw[:,iw])*sw[:,iw]/sigsq[:])
            bot=np.sum(sw[:,iw]*sw[:,iw]/sigsq[:])
            
            if (top == 0):
             sk = 0
            else:
             sk=top/bot
            
            sigsk=np.sqrt(1./bot)
            #print top,bot,xdup.sum(),ck
            #print iteration,iw
            #print cw[:10,iw]
            #print sw[:10,iw]
            #print xdup[:10]
            #print sigsq[:10]
            #print np.sum((xdup[:]-ck*cw[:,iw])*sw[:,iw]/sigsq[:])
            #print sigsq[:]
            #raw_input()
### cos parameter
            top=np.sum((xdup[:]-sk*sw[:,iw])*cw[:,iw]/sigsq[:])
            bot=np.sum(cw[:,iw]*cw[:,iw]/sigsq[:])
            ck=top/bot
            sigck=np.sqrt(1./bot)
            
            sumold=ck+sk
            if iteration > 1:
                convtest=abs((sk+ck - sumold)/sumold)
                if convtest < 0.0001:
                    break
        
        xdup=xdup-sk*sw[:,iw]-ck*cw[:,iw]   # subtract next frequency             
        if abs(sk) > overflow:
            sk=0
       
        if abs(ck) > overflow:
            ck=0
            
        #print iw,w[iw],sk,ck,top
        #print ''
        #raw_input('space')
        p[iw,0]=sk
        p[iw,1]=ck
        sigp[iw,0]=sigsk
        sigp[iw,1]=sigck        
        
        p[-1.:]=0.0
        sigp[-1,:]=0.0
        #print iw, p[iw,:],sigp[iw,:],'myios.py',np.sqrt(sigsq[0]/np.sum(np.sin(w[iw]*t[:])**2)),w[iw],np.sum(np.sin(w[iw]*t[:])**2)/nt           
    return(p,sigp)
    
    
    
    


#!!!!!!!!!!!!!!!! test below
#twopi=np.pi*2
#w0=2*np.pi*0.5
#tlo=0.0
#thi=30.0
#dt=0.1
#t=np.arange(tlo,thi+dt,dt)
#nt=t.shape[0]
#
#tmean=15.0
#tsd=5.0
#x=1/np.sqrt(twopi)/tsd* exp(-0.5*((t-tmean)/tsd)**2)
#
##3.2*np.sin(w0*t)
#
#
#xmean=x.mean()
#xsig=np.ones(nt)*1.0
#sk1=1.0
#ck1=1.0
#wlo=twopi*1./(t[-1]-t[0])
#whi=twopi*1./(2.*dt)
#nw=int(np.ceil(0.5*nt))
#dw=(whi-wlo)/(nw-1)
#w=wlo+np.arange(nw)*dw
#p=np.zeros((nw,2))
#sigp=np.zeros((nw,2))
#
#p,sigp=ios(t,x,xsig,w,sk1,ck1)
#
#tmodlo=t[0]
#tmodhi=t[-1]
#dtmod=dt/10
#ntmod=int(np.ceil((tmodhi-tmodlo)/(dtmod))+1)
#tmod=tmodlo+np.arange(ntmod)*dtmod
#xmod=np.zeros(ntmod)
#for i in range(ntmod):
#    xmod[i]=xmean+np.sum(p[:,0]*np.sin(w[:]*tmod[i])+p[:,1]*np.cos(w[:]*tmod[i]))
#
#plot(t,x,linestyle='none',marker='o')
#plot(tmod,xmod)
#show()
#