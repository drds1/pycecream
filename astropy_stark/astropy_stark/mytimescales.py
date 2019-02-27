import numpy as np



def tlc(embh,rld):
 rm = 2.58e13*rld
 rsc = 3000*embh
 c   = 2.9979e8
 tday = rm/c

 return(tday)




def tdyn(embh,rld):
 gnewt = 6.673e-11
 rm    = rld*2.58e13
 rm3 = rm**3
 
 emkg = 2.e30*embh
 return(np.sqrt(rm3/(gnewt*emkg)))

def tth(embh,rld,alpha=0.03):
 tdy = tdyn(embh,rld)
 return(tdy/alpha)
 

def tvisc(embh,rld,alpha=0.03):
 
 rm = 2.58e13*rld
 rcm = rm*100
 r18 = rcm/1.e18
 secyr = 24*3600*365.25
 return(secyr * 10**10 * (alpha/0.03) * r18**1.5)
 
 
 
 
 
 
 
 

embh = 1.e8
rld = 1
secday = 86400

print 'tlc, tdyn, tth, tvisc'

tl = tlc(embh,rld)/secday
tt = tth(embh,rld)/secday
tv = tvisc(embh,rld)/secday
td  = tdyn(embh,rld)/secday

print tl,td,tt,tv



