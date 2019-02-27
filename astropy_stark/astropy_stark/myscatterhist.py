import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import mymedquart


def myscatterhist(x,y,col='b',l='',xlab='',ylab='', xlo = 666,xhi = 666,ylo = 666,yhi = 666,nbins = 100):
 
 paperplot = 1
 print 'inside function', l
 nullfmt   = NullFormatter()         # no labels
 
 # definitions for the axes 
 left, width = 0.1, 0.25
 bottom, height = 0.1, 0.25
 bottom_h = left_h = left+width+0.02
 
 

 rect_scatter = [left, bottom, width, height]
 rect_histx = [left, bottom_h, width, 0.2]
 rect_histy = [left_h, bottom, 0.2, height]
 

 # start with a rectangular Figure
 #plt.figure(1, figsize=(8,8))
 
 axScatter = plt.axes(rect_scatter)
 axHistx = plt.axes(rect_histx)
 axHisty = plt.axes(rect_histy)
 
 # no labels
 axHistx.xaxis.set_major_formatter(nullfmt)
 axHisty.yaxis.set_major_formatter(nullfmt)
 
 # the scatter plot:
 axScatter.scatter(x, y, color = col, s=0.5)#np.ones(x.shape[0])*1)
 # add vertical lines to the plot to indicate the mean and sd of the parameter
 medloupx = mymedquart.mymedquart(x,0.68)
 lo = medloupx[1]
 up = medloupx[2]
 med = medloupx[0]
 #axScatter.vlines(med,ylo,yhi,color=col,linewidth=3)
 #axScatter.vlines(lo,ylo,yhi,color=col)
 #axScatter.vlines(up,ylo,yhi,color=col)
 
 medloupy = mymedquart.mymedquart(y,0.68)
 lo = medloupy[1]
 up = medloupy[2]
 med = medloupy[0]
 #axScatter.hlines(med,xlo,xhi,color=col,linewidth=3)
 #axScatter.hlines(lo,xlo,xhi,color=col)
 #axScatter.hlines(up,xlo,xhi,color=col)
 
 
 # now determine nice limits by hand:
 if (xlo == 666 & xhi == 666):
  xhi = x.max()
  xlo = x.min()
 
 if (ylo == 666 & yhi == 666):
  yhi = y.max()
  ylo = y.min()
 
 binwidthx = (xhi - xlo)/nbins
 binwidthy = (yhi - ylo)/nbins
 
 
 axScatter.set_xlim( (xlo, xhi) )
 axScatter.set_ylim( (ylo, yhi) )
 axScatter.set_xlabel(xlab)
 axScatter.set_ylabel(ylab)
 
 #we need a sensible way to define the number of bins. The routine automatically tries to include everything
 binsx = np.arange(xlo, xhi + binwidthx, binwidthx)
 binsy = np.arange(ylo, yhi + binwidthy, binwidthy)
 
 #print '10',binsx.shape[0],binsy.shape[0],xlo,xhi,ylo,yhi,binwidthx, binwidthy
 #raw_input()
 axHistx.hist(x, bins=binsx, color = col, label=l)
 #axHistx.legend(l)
 
 
 #print '11'
 axHisty.hist(y, bins=binsy, orientation='horizontal', color = col)
 
 axHistx.set_xlim( axScatter.get_xlim() )
 axHisty.set_ylim( axScatter.get_ylim() )

 return axHistx


col =['r','b','k','p','y']
#fig = plt.figure()

plt.figure(1, figsize=(8,8))
for i in range(1):
 #ax1 = fig.add_subplot(2,2,i+1)
 a = np.random.randn(1000) + 25*i
 b = np.random.randn(1000) + 2*i

 
 ax1 = myscatterhist(a,b,ifig,nalong,nvert)

plt.show()

