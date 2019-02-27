# Function to calculate a skewed gaussian with the mean, sig, skew and kurt inputted
# INPUT: tau[1:Ntau],tau0,sig0,pskew,pkurt
# OUTPUT: psitau[1:Ntau]
import numpy as np

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  m3 function
def m3(y):
 yy=y*y
 yyy=yy*y
 m3=1./np.sqrt(3.)*(2*yyy-3*y)

 return(m3)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  m4 function
def m4(y):

 yy=y*y
 yyyy=yy*yy

 m4=1./np.sqrt(24.) * (4*yyyy - 12*yy +3*y)
 
 return(m4)
 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! tfskewgauss function
def tfskewgauss(tau,p):#tau0,sig0,pskew,pkurt):
 
 psi0 = p[0]
 tau0 = p[1]
 sig0 = p[2]
 pskew = p[3]
 pkurt = p[4]
 
 y = (tau - tau0)/sig0

 psitau = psi0 * np.exp(-y*y/2) * (1. + pskew * m3(y) + pkurt*m4(y) / ( sig0 *np.sqrt(2*np.pi) ) )

 return (psitau)



## tested works 2/10/2014
#
## test the above function below
#
#p0 = 1
#tau0 = 10.0
#sig = 3.0
#skew = 0.2
#kurt = 0.0
#
#tau = np.arange(0,30,0.01)
#psitau = tfskewgauss(tau, (p0,tau0, sig, skew, kurt))
#dtau = tau[1] - tau[0]
#psiarea = np.sum(psitau*dtau)
#top = np.sum(psitau*tau)
#bot = np.sum(psitau)
#taumean = top/bot
#
#from pylab import *
#plot(tau,psitau)
#show()
#
##!! calculate the moments and compare
#import mystatmoments as ms
#
#opsig2 = ms.mystatmoments(tau,psitau/psiarea,taumean,2)
#opskew = ms.mystatmoments(tau,psitau/psiarea,taumean,3)
#opkurt = ms.mystatmoments(tau,psitau/psiarea,taumean,4)


