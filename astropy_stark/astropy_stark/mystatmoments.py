# routine to calculate the nth cenral moment of an input distribution
# Bug fixed 15/10/14 incorrect calculation of skew and kurt fixed now
import numpy as np
import pylab as plt

def mysig2(xdist,ydist,x0,area,a):
 
 sum = np.sum(ydist*a*a)
 return(sum/area)


def mystatmoments(xdist,ydist,x0,nth):
 
 a = xdist - x0
 ndist = xdist.shape[0]
 dx = (xdist[-1] - xdist[0])/(ndist-1)
 area = np.sum(ydist)
 
 sig2 = mysig2(xdist,ydist,x0,area,a)
 if (nth == 2): 
  answer = sig2
 elif (nth == 3):
  aa = a*a
  sum = np.sum(ydist*aa*a)
  answer = sum/area/(sig2**1.5)  #before 15/10/14 just divided by sig2 and not sig^n
 elif (nth ==4):
  aa= a*a
  sum = np.sum(ydist*aa*aa)
  answer = sum/area/(sig2*sig2)
 
 return(answer)
 
 
##!!!!!! test here
#x = np.arange(0,100,0.00001)
#
#x0 = 15.0
#sig = 4.0
#
#a = 10/(np.sqrt(2*np.pi*sig*sig))
#
#y = a * np.exp(-0.5*(x-x0)**2/(sig*sig))
##
#
#opsig = mystatmoments(x,y,x0,2)
#opskew = mystatmoments(x,y,x0,3)
#opkurt = mystatmoments(x,y,x0,4)
##
#print opsig, opskew, opkurt, x.shape[0]



 