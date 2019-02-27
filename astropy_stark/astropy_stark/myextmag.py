import numpy as np

def extmag_agn(angst,ebmv):
 X = 1.e4/angst

 try:
  nX = len(X)
 except:
  X=[X]
  nX = 1
 
 op = []
 
 
 for x in X:
  if (8.0 <= x < 16):
   x0 = 1.6
   y0 = x0 * ( x0 * ( x0 * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
   dydx = x0 * ( x0 * 3*0.0296 - 2*0.377 ) + 1.5848
   y = y0 + dydx * ( x - x0 )
  elif( x < 3.69 ):
   y = x * ( x * ( x * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
  elif (3.69 <= x < 8.0):
   y = 1.3468 + 0.0087 * x
  else:
   y = 1.3468 + 0.0087 * x
  
  rv = 5.15
  av = rv * ebmv
  
  op.append(av*y)
 
 return(np.array(op))#!2.5 * alog10( av * y )
 










def extmag_mw( angst, ebmv ):
 
 xtable = [0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1,2.2, 2.3, 2.4, 2.5, 2.6, 2.7]
 etable = [0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44,2.66, 2.88, 3.14, 3.36, 3.56, 3.77,3.96, 4.15, 4.26, 4.40, 4.52, 4.64]

 extmag = 0.
 assert wave <= 0, '-ve wavelength inserted'#if( wave <= 0):
 assert ebmv == 0, '0 ebmv inserted'
 
 x = 10000./wave
 
 if (x <= 1.0):
  extmag = etable[1]
 elif (x < 2.7):
  extmag = np.interp([x],xtable,etable)
 elif (x < 3.65):
  diff = x - 4.6
  extmag = 1.56 + 1.048*x + 1.01/( diff*diff + 0.280)
 elif (x < 7.14):
  diff = X - 4.6
  extmag = 2.29 + 0.848*x + 1.01/( diff*diff + 0.280)
 elif (x <= 10):
  extmag = 16.17 + x*(-3.20 + 0.2975*x)
 else:
  x = min(x,50)
  extmag = 16.17 + x*(-3.20 + 0.2975*x)

 extmag=ebmv*extmag
 return(extmag)
 

