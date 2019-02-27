import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import myscatterhist
import mymedquart
import os

# hardwired parameters  !!!!!!!!!!!!!!!!!!


def myscatterhistbasic(dirmain,subfile,parmlab,idx1,idx2,chop,subfilelab):

 nbins = 100

 col = ['b','r','c','m','g','y','k','0.75','b','r','c','m','g']
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 os.chdir(dirmain)


 nplots = len(subfile)
 parm1 = []#np.zeros((0,nplots))
 parm2 = []#np.zeros((0,nplots))


# Load the data first to decide what plot limits to use
 for idx in range(nplots):
  
  dat = np.loadtxt(subfile[idx])
  parm1.append(dat[chop:,idx1])
  parm2.append(dat[chop:,idx2])
#
 
  print 'running scatterplot code',idx, len(parm1[0])
 #parm1 = np.random.randn(200000) + 10
 #parm2 = np.random.randn(200000)

# determine plot limits
 xloall=[]
 xhiall=[]
 yloall=[]
 yhiall=[]
 opx=[]
 opy=[]
 for idx in range(nplots):
  medloupx = mymedquart.mymedquart(parm1[idx],0.68)

  
  medloupy = mymedquart.mymedquart(parm2[idx],0.68)
  
  print medloupx
  print medloupy
  print 'med loup above are the medloup'
  raw_input()
  xlo = medloupx[1]
  xhi = medloupx[2]
  xmed = medloupx[0]
  xloall.append(xlo - np.abs(xlo)/10)
  xhiall.append(xhi + np.abs(xhi)/10)
  
  ylo = medloupy[1]
  yhi = medloupy[2]
  ymed = medloupy[0] 
  yloall.append(ylo - np.abs(ylo)/10)
  yhiall.append(yhi + np.abs(yhi)/10) 
  opx.append([xmed,xlo,xhi])
  opy.append([ymed,ylo,yhi])
 
 
 xlo = min(xloall)
 xhi = max(xhiall)
 ylo = min(yloall)
 yhi = max(yhiall)
 print xloall
 print xhiall
 print yloall
 print yhiall
 print 'above are the lowest pointsof all the data'
 
 
 raw_input()

 for idx in range(nplots):
  axHistx = myscatterhist.myscatterhist(parm1[idx],parm2[idx],col[idx],np.str(subfilelab[idx]),parmlab[idx1],parmlab[idx2],xlo,xhi,ylo,yhi,nbins)
 
  print 'finished', subfile[idx],subfilelab[idx]
#print plot
 axHistx.legend(prop={'size':10})
 
 


#os.chdir('../..')


 print 'limits', xlo,xhi,ylo,yhi
 plt.xticks(rotation = 90)

 plt.show()



#dirmain = './inclinations'
#subfile = ['op_tfbop_0.720','op_tfbop_0.730','op_tfbop_0.740','op_tfbop_0.750','op_tfbop_0.760','op_tfbop_0.770','op_tfbop_0.780']#,'11_of_12']
#parmlab = ['P0','w0','alpha','d_alpha','Error Bar Expansion Factor','Offset Level']


