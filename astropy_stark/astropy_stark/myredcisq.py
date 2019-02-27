import numpy as np


def myredcisq(tdat,xdat,sigdat,tmod,xmod,npar = 0):

 ndat = np.shape(tdat)[0]
 xint = np.interp(tdat,tmod,xmod)
 cisq = np.sum( ((xint - xdat)/sigdat)**2)
 redcisq = cisq/(ndat-npar)
 
 return(cisq,redcisq)