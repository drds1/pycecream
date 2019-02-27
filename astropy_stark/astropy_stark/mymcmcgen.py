# tested and works 7/10/14 
# general mcmc code python
# IP data: x,y,sigy, func_mod: the model function, nparm: the number of parameters, pstep: best guess of teh step sizes, p: best starting guesses of the parameters, psafe: 0 if step uniformly, 1 if step in log

import numpy as np
import os
import pylab as plt




def mymcmc(x,y,ysig,func_mod,nparm, p, pstep, psafe, p_primean, p_prisig, positivity):
 
 nits = 5000
 checknum = 10
 checkfrac = 0.25
 ploton = 0

 pold = np.zeros(nparm)
 pAT = np.zeros(nparm)
 cisqsave=np.zeros(0)
 parmsave=np.zeros((0,nparm))
 

 checkarange = np.arange(nparm)
 
#optional for plotting purposes (the ion allows the plot to be updated live)
 if (ploton == 1):
  plt.ion()
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x,y)
  yplot = y*0.0
  line1, = ax.plot(x,yplot)
# 
 
# for the bof 
 #idx_pri_log = np.where((p_primean > 0) and ( psafe == 1) )
 #we cant do the above as np.where doesn't work for two arrays
 #instead do it like this
 a = np.logical_and((p_primean > 0), (psafe == 1))
 idx_pri_log=np.where(a == 1)[0]
 
 a = np.logical_and((p_primean > 0), (psafe == 0))
 idx_pri_lin = np.where(a == 1)[0]
 
 
 #cant_start_safe_parms_at_0
 idx_safe_check = np.where(psafe == 1)[0]
 if ( np.any(p[idx_safe_check] <= 0) == 1):
  print 'You cannot make starting values == 0 for parameters stepped in log!'
  raw_input()

 
 
#start the iteration loop 
 
 for iteration in range(nits):


# change step sizes
# keep step sizes reasonable
  reset_check = np.mod(iteration,checknum)
  if (reset_check == 0):
   idxup = np.where(pAT/checknum > checkfrac)[0] #these steps are being doubled
   idxdown = np.where(checkarange != idxup)[0] #these steps are bing halved
   pstep[idxup] = 2*pstep[idxup]
   pstep[idxdown] = 0.5*pstep[idxdown]
   pAT[:] = 0
# finished changing step sizes

#step in the parameter or log of the parameter
  for ip in range(nparm):
   pip = p[ip]
   pipstep = pstep[ip]
   pold[ip] = pip

# step in the parameter (or the of the parameter if psafe =1)  
   a = np.random.randn(1)*pipstep 
   if (psafe[ip] ==1):
    logp = np.log10(pip)
    pip = 10**(logp + a)
   else:
    pip = pip + a
  
   #if (ip == 1):
    #print pip,a,logp, iteration, 'psafe bug'
    #raw_input()
    
  #print 'parm test', pip, a, p[ip]
   p[ip] = pip
  
  
# evaluate the model with the new parameters. p is the 1d array of parameters
   ymod = func_mod(x,p)
   #print ymod
   if (positivity == 1):
    ymod = np.exp(ymod) - 1
   #print ymod
   #print ''
# calculate the cisq of the new model wrt data
   a = y - ymod
   cisq = np.sum( a*a/(ysig*ysig) )

# add on any priors   
   a = (p[idx_pri_lin]- p_primean[idx_pri_lin])/p_prisig[idx_pri_lin]
   bof = np.sum(a*a)
   
   #print a
   #print bof
   #print 'bof bug'
   #raw_input()
   
   a = (np.log10(p[idx_pri_log]) - p_primean[idx_pri_log])/p_prisig[idx_pri_log]
   bof = bof +np.sum(a*a) 
   
   #print p[idx_pri_log], p_primean[idx_pri_log], p_prisig[idx_pri_log]
   #print a
   #print bof
   #print 'bof bug 2'
   #raw_input()
    
   bof = bof + cisq
   
   
# decide to oaccept
   cisqnew = cisq
   bofnew = bof
   if ( (iteration == 0) and (ip == 0) ):
    cisqold = cisq
    bofold = bof
# accepts 
  #print cisqnew,cisqold,'cisqtest' 
  #raw_input()
   if ( (bofnew < bofold) or (iteration == 1) or (np.exp(-0.5*(bofnew-bofold)) > np.random.rand(1)) ):
    pAT[ip] = pAT[ip]+1
    cisqold = cisqnew
    bofold  = bofnew
# reject   
   else:
    p[ip] = pold[ip]

  cisqsave = np.append(cisqsave,cisqold) #save the bof
  parmsave = np.vstack( [parmsave, [p[:]] ] )

 
  print iteration
  print p
  print pstep
  print bofold
  #raw_input()
# optional bit (plot the result of each iteration for comparison
  if (ploton == 1):
   yplot = func_mod(x,p)
   if (positivity == 1):
    yplot=np.exp(yplot) - 1
   #print 'poditive'
   #xmodtest = np.arange(-5,5,0.1)
   #ymodtest=np.exp(func_mod(xmodtest
   line1.set_ydata(yplot)
   fig.canvas.draw()
  #print yplot[:3]
  #print y[:3]
  #print x[:3]
  #raw_input()
  
 return(cisqsave,parmsave)
 


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
#test inputs (these will be turned into inputs for a definition
#x = np.arange(10)
#y = 2*x + 3
#nx = x.shape[0]
#ysig = np.ones(nx)*x.mean()/10
#
#
#
#
#def f(xin,a,b):
# return(a*xin+b)
#
#nparm = 2
#pstep = np.ones(nparm)
#p=np.zeros(nparm)
#psafe = np.zeros(nparm)
#
#op = mymcmc(x,y,ysig,f,nparm, p, pstep, psafe)
#