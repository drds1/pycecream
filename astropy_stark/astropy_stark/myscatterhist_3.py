#My scatter hist 2 allows Multiple 2d parmparm plots to be plotted on the same window
#Useful for methods paper
# need to make better to allow for par file input write this in different subroutine
# Need to upgrade to allow for multiple plots to be shown on same plot
# allow for forcing x and y lims to be same

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import mymedquart
import scipy
import scipy.stats
import myreadlines
import mystats as ms

def myscatterhist_2(xr,y,ialong,ivert,nalong,nvert,col='b',l='',xlab='',ylab='', xlo = 666,xhi = 666,ylo = 666,yhi = 666,nbins = 100, labstyle='bl', translucency = 0.0, lineplot=[0,0,0,0,0,0,0],title=[''],special =[[0],[0]],sigconts=[], truth = [666,666],showconts=[],scatterskip=1):
 
 
# plot log axis if necessary 
 if (special[0][0] == -2):
  x = np.log10(xr)
 else:
  x = xr
 
 #print x.shape,x.min(),x.max(), special
 #raw_input()
 #x = x[:-500:]
 #y = y[:-500:]
 
 paperplot = 1
 #print 'inside function', l
 #raw_input()
 
 nullfmt   = NullFormatter()         # no labels
 lineplot = np.array(lineplot)#plot the comparison with ccf and cream in anotherwindow
 plotlp   = np.mean(lineplot)
 
 # definitions for the axes 
 
 left0 = 0.1
 right0 = 0.9
 bottom0 = 0.1
 top0 = 0.9
 dspace = 0.05
 scatterfrac = 0.5
 dlineplot = 0.04
 
 onesubsf    = 1.-scatterfrac
 width = (1. - left0 - (1.-right0) - (nalong-1)*dspace)/nalong
 height =( 1. - bottom0 - (1.-top0) -(nvert-1)*dspace)/nvert
 
 left   = left0 + ialong*(dspace + width)
 bottom = top0 - (ivert+1)* height - (ivert)*dspace
 
 width_s   = scatterfrac*width
 height_s  = scatterfrac*height 
 left_hy   = left   + width_s
 width_hy  = width - width_s
 
 bottom_hx = bottom + height_s
 height_hy = height - height_s

 

 rect_scatter = [left, bottom, width_s, height_s]
 rect_histx = [left, bottom_hx, width_s, onesubsf*height]
 rect_histy = [left_hy, bottom, onesubsf*width, height_s]
 rect_line  = [left_hy + dlineplot, bottom_hx + dlineplot, onesubsf*width - dlineplot, onesubsf*height - dlineplot]
 
 #print rect_histx
 #print rect_scatter
 #print width, height#,(1. - left0 - right0 - (nalong-1)*dspace),dspace,1. - left0 - right0
 #print bottom,left,ivert,dspace
 #print height,height_s,height_hy
 #raw_input()

 # start with a rectangular Figure
 #plt.figure(1, figsize=(8,8))
 
 axScatter = plt.axes(rect_scatter)
 axHistx = plt.axes(rect_histx)
 axHisty = plt.axes(rect_histy)

 
 # no labels
 #axHistx.xaxis.set_major_formatter(nullfmt)
 #axHisty.yaxis.set_major_formatter(nullfmt)
 
 # if labstyle = 'bl' turn off tick labeling for all but x for lower plots and y for left plots
 axHistx.set_yticklabels(' ')
 axHistx.set_yticks([])
 axHisty.set_yticklabels(' ')
 axHistx.set_xticklabels(' ')
 axHisty.set_xticklabels(' ')
 axHisty.set_xticks([])
 
 if (labstyle == 'bl'):
  if (ivert != nvert - 1):
   axScatter.set_xticklabels('')
  if (ialong > 0):
   axScatter.set_yticklabels('')
   
 
 
 # the scatter plot:
 axScatter.scatter(x[::scatterskip], y[::scatterskip], color = col, s=0.000015)#np.ones(x.shape[0])*1)

 
 
 #plot the true values if necessary
 if (truth[0] != 666):
  axScatter.vlines(truth[0],ymin=ylo,ymax=yhi,color='k')
 if (truth[1] != 666):
  axScatter.hlines(truth[1],xmin=xlo,xmax=xhi,color='k')
 
 
 
 #print 'scatter plotted', xlo, xhi, ylo, yhi, x.mean(), y.mean()
 # add vertical lines to the plot to indicate the mean and sd of the parameter
 #print x
 #raw_input()
 medloupx = mymedquart.mymedquart(x,0.68)
 
 #print medloupx

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

 if (int(xlo) == 666 & int(xhi) == 666):
  xhi = x.max()
  xlo = x.min()
  
 if (int(ylo) == 666 & int(yhi) == 666):
  yhi = y.max()
  ylo = y.min()
 
 binwidthx = 1.*(xhi - xlo)/nbins
 binwidthy = 1.*(yhi - ylo)/nbins
 

 axScatter.set_xlim( (xlo, xhi) )
 axScatter.set_ylim( (ylo, yhi) )
 axScatter.set_xlabel(xlab)
 axScatter.set_ylabel(ylab)
 
 if (len(special[1]) > 0):
  temp = special[1]
  if (temp[0] == -3):
   axScatter.set_yticks(np.cos(np.pi*np.array(temp[1:])/180))
   axScatter.set_yticklabels(temp[1:])

 if (len(special[0]) > 0):
  temp = special[0]
  if (temp[0] == -2):
   axScatter.set_xticks(temp[1:])
   axScatter.set_xticklabels(temp[1:])


## plot contours if necessary

 nconts = len(sigconts)
 #print 'OOOOAHHHHHHHHHH!!!!!!!',nconts, xlo, xhi, binwidthx,(xhi-xlo)/binwidthx, ylo, yhi, binwidthy,(yhi-ylo)/binwidthy
 if (nconts > 0):
  # fit a KDE to the data
  
  
  ntemp =  len(x)
  #print ntemp,len(y)
  
  xtemp = np.array(x)
  ytemp = np.array(y)
  idx_xnan = np.isnan(xtemp)
  idx_x0= x == 0
  idx_y0= y == 0
  idx_xinf= abs(x) == np.inf
  idx_yinf= abs(y) == np.inf
  idx_ynan = np.isnan(ytemp)
  xtmean = np.median(xtemp)
  ytmean = np.median(ytemp)  
  xtemp[idx_xnan] = xtmean
  xtemp[idx_x0] = xtmean
  ytemp[idx_y0] = ytmean
  ytemp[idx_ynan] = ytmean
  ytemp[idx_yinf] = ytmean
  xtemp[idx_xinf] = xtmean
  x = np.array(xtemp)
  y = np.array(ytemp)
  
  
  

  both = np.array((x,y)).T
  u_x = ms.stat_uncert(x,sig=0.32)
  u_y = ms.stat_uncert(y,sig=0.32)
  
  u_y_2=ms.sigtheta(y,u_y[1])
  
  print 'UNCERT mmdot', '%.2f' % u_x[0], '%.2f' % u_x[1]
  print 'UNCERT inc', '%.2f' % u_y_2[0], '%.2f' % u_y_2[1]
  print ''
  try:
   pdf1=scipy.stats.kde.gaussian_kde(both.T)
  except:
   for i in range(ntemp):
    print x[i],y[i]
    if ( (x[i] == np.nan) or (x[i] == -np.inf) or (y[i] == np.nan) or ( y[i] ==  -np.inf)): 
     print ' It broke!! ', i, x[i], y[i]
     raw_input()
  
  # create a grid over which we can evaluate pdf
  q,w=np.meshgrid(np.arange(xlo,xhi,binwidthx), np.arange(ylo,yhi,binwidthy))
  r1=pdf1([q.flatten(),w.flatten()])
  # sample the pdf and find the value at the 95th percentile
  s1=[]
  for ic in range(nconts):
   s1.append(scipy.stats.scoreatpercentile(pdf1(pdf1.resample(1000)), sigconts[ic]))
   
  # reshape back to 2d
  r1.shape=(q.shape[0],q.shape[1])

  for ic in range(nconts):
   #print 'OOOOAHHHHHHHHHH!!!!!!!', s1[ic],[sigconts[ic]],ic, nconts, xlo, xhi, binwidthx, ylo, yhi, binwidthy
   axScatter.contour(np.arange(xlo,xhi,binwidthx), np.arange(ylo,yhi,binwidthy), r1, [s1[ic]],colors=col)
   #print 'OOOOAHHHHHHHHHH!!!!!!!', [sigconts[ic]],ic, nconts, xlo, xhi, binwidthx, ylo, yhi, binwidthy
## 
 
 #we need a sensible way to define the number of bins. The routine automatically tries to include everything
 #print xlo, xhi, ylo, yhi, binwidthx, binwidthy, nbins
 #raw_input()
 
 #print xlo, xhi, binwidthx, nbins
 #print ylo, yhi, binwidthy, nbins
 #raw_input()
 binsx = np.arange(xlo, xhi + binwidthx, binwidthx)
 binsy = np.arange(ylo, yhi + binwidthy, binwidthy)
 
 #print '10',binsx.shape[0],binsy.shape[0],xlo,xhi,ylo,yhi,binwidthx, binwidthy
 #print x
 #print translucency, binsx, col
 #raw_input()
 try:
  axHistx.hist(x, bins=binsx, color = col, label=l, alpha = translucency,histtype='step')
  #print 'Histogram good x', x.min(), x.max() 
 except: 
  axHistx.hist(x, bins=binsx, color = col, label=l, alpha = 0.3)
  print 'Histogram problem x', x.min(), x.max() 
 #axHistx.legend(l)
 
 
 #print '11'
 try:
  axHisty.hist(y, bins=binsy, orientation='horizontal', color = col, alpha = translucency,histtype ='step')
 except:
  axHisty.hist(y, bins=binsy, orientation='horizontal', color = col, alpha = 0.3)
  print 'Histogram problem y', y.min(), y.max()
  
 
 axHisty.set_xlabel(title)
 axHisty.xaxis.set_label_position('top')
 axHisty.xaxis.labelpad = 20
 
 axHistx.set_xlim( axScatter.get_xlim() )
 axHisty.set_ylim( axScatter.get_ylim() )

 yt = axHistx.get_ylim() 
 axHistx.vlines(truth[0],0,yt[1],color='k')
 yr = axHisty.get_xlim()
 axHisty.hlines(truth[1],0,yr[1],color='k')
 #print 'xlims', xlo,xhi, np.mean(x), np.std(x),x.shape[0]
 #print 'ylims', ylo,yhi, np.mean(y), np.std(y), y.shape[0]
 
 
 
 
 if (plotlp > 0):
  axline  = plt.axes(rect_line)
  #plot the true value
  axline.vlines(lineplot[0],lineplot[5], lineplot[6], color = 'k')
  axline.hlines(lineplot[0],lineplot[5], lineplot[6], color = 'k')
  
  #cram pt with errorbar then ccf with error bar
  axline.errorbar(lineplot[1],lineplot[3],xerr=lineplot[2],yerr=lineplot[4],color = col, linestyle = 'None')
  axline.plot([lineplot[5],lineplot[6]],[lineplot[5],lineplot[6]], color = 'k') #slope of 1
  axline.set_xlim(0,lineplot[5])
  axline.set_ylim(0,lineplot[6])
  axline.set_xticks([lineplot[5],lineplot[6]])
  axline.set_yticks([lineplot[5],lineplot[6]])
  
  #print lineplot

 return axHistx








#new function to plt the multiple correlation plots with input from a text file
#if lab style = 'bl', axis tick labels are placed only on xaxis of lower plots and y-axis of left plots

def myscatterhistfile(filename,colhist,nvert,nalong,parcols,sublabfile = '',xlim=[], ylim=[], labsin=['',''],labstyle='bl', spec=[[0],[0]], nsplit = 1, nperhist=1, sigconts=[], showconts=[],truth = [666,666], scatterskip=1,nbins=100):
 
 #number of datasets per correlatin plot determined by number of entries in colour list i.e ['r'] just 1 or ['r','b'] two etc 
 if (nsplit>1):
  nperhist = 1
 else:
  nperhist = len(colhist)

 if (nperhist > 1 or nsplit > 1):
  translucency = 0.5
 else:
  translucency = 1.0
 
 with file(filename) as f:
  fname = f.read().splitlines()

 nf  = len(fname)

 taumean_ccf=[]
 taumean_cream=[]
 taumean_real=[]
 tit = []
 #determine if we are making the special ccf plots if ccf and scream data are included
 
 lineplot_check = 0
 for i in range(nf):
  a = str.split(fname[i], ',')

  fname[i] = a[0]
  
  if (len(a) > 5):
   taumean_ccf.append((a[5], a[6]))
   taumean_cream.append((a[3], a[4]))
   taumean_real.append(a[1])
   lineplot_check = 1
   
  if (len(a) > 1):
   tit.append(a[1])
  else:
   tit.append('')
 
 if (lineplot_check > 0): 
  taumean_ccf   = np.double(np.array(taumean_ccf))
  taumean_cream = np.double(np.array(taumean_cream))
  taumean_real  = np.double(np.array(taumean_real))
 

  taumeanccf_min = np.min(taumean_ccf[:,0] - taumean_ccf[:,1] )
  taumeanccf_max = np.max(taumean_ccf[:,0] + taumean_ccf[:,1] )

  taumeancream_min = np.min(taumean_cream[:,0]- taumean_cream[:,1] )
  taumeancream_max = np.max(taumean_cream[:,0] + taumean_cream[:,1] )
 
  lineplotmin = 0#np.min((taumeancream_min,taumeanccf_min))

  
  lineplotmax = np.max((taumeancream_max,taumeanccf_max, taumean_real[i]))*1.2
  
   
 nsublab = len(sublabfile)
 if (nsublab > 0):
  with file(sublabfile) as f:
   sublab = f.read().splitlines()  
 


 
 dat=[]
 #!load the data
 for i in range(nf):
  print fname[i]
  try:
   a   = np.loadtxt(fname[i])
  except:
   a   = myreadlines.myreadlines(fname[i])
  ncol = a[0,:].shape[0]
  for i2 in range(ncol): #kde code for confidence curves doesn't like nans. Change these to median of data. Dodgy but no one will know!
   cmed = np.median(a[:,i2])
   nancheck = np.isnan(a[:,i2])
   a[nancheck,i2] = cmed
   
  dat.append(a)
  #print i, fname[i]
  #print a
  #print a.max(), a.min(), np.std(a)
  #raw_input()
  
 plt.figure(1, figsize=(8,8))
 
 idx = 0
 labs = ['','']
 for ia in range(nalong):
  if (ia == 0):
   ylab = labsin[1]
  else:
   ylab = ''
  
  for iv in range(nvert):
   if (iv == nvert-1):
    xlab = labsin[0]
   else:
    xlab = ''
   ivert  = iv 
   ialong = ia
   

   for i2 in range(nperhist):
    icol = i2
    ccfcomp = [0,0,0,0,0,0,0]
    if (lineplot_check > 0):
     ccfcomp[0] = taumean_real[idx]
     ccfcomp[1] = taumean_cream[idx][0]
     ccfcomp[2] = taumean_cream[idx][1]
     ccfcomp[3] = taumean_ccf[idx][0]
     ccfcomp[4] = taumean_ccf[idx][1]
     ccfcomp[5] = lineplotmin
     ccfcomp[6] = lineplotmax    


    for isplit in range(nsplit):
     if (nsplit > 1):
      icol = isplit
     
     #print idx, len(dat), i2, nperhist
     nd = len(dat[idx])
     npersamp = nd /nsplit
     samp = np.arange(isplit*npersamp,(isplit+1)*npersamp) 
     
     if (showconts[isplit] == 0):
      sigcontsnew = []
     else:
      sigcontsnew = list(sigconts)
     if (len(xlim) > 0 ):
      ax1 =  myscatterhist_2(dat[idx][samp,parcols[0]],dat[idx][samp,parcols[1]],ialong,ivert,nalong,nvert,xlo = xlim[0], xhi = xlim[1], ylo = ylim[0], yhi = ylim[1],col = colhist[icol], labstyle = labstyle, translucency = translucency, lineplot = ccfcomp,title = tit[idx], special = spec, sigconts=sigcontsnew, truth = truth,xlab=xlab,ylab=ylab,scatterskip=scatterskip,nbins=nbins)
     else:
      ax1 = myscatterhist_2(dat[idx][samp,parcols[0]],dat[idx][samp,parcols[1]],ialong,ivert,nalong,nvert,col=colhist[icol], labstyle = labstyle, translucency = translucency, lineplot = ccfcomp,title = tit[idx],special=spec,sigconts=sigcontsnew, truth = truth,xlab=xlab,ylab=ylab,scatterskip=scatterskip,nbins=nbins)
    #print dat[idx]
    #print dat[idx].min(), dat[idx].max(), np.mean(dat[idx]), idx
    
    idx = idx + 1
    

 return()
 



#n = 1000
#nf = 8
#nvert = 4
#nalong =2
#
#colhist=['r']#,'b']
#nperhist = 1
#
#dattest = np.zeros((n,2))
#
#
#fnam = 'histmultidat.dat'
#f = open(fnam,'w')
#for i in range(nf):
# 
# fname = 'test_myscatterhist_'+str(i)+'.dat' +' 1.1 1.0 0.5 1.2 0.5'
# f.write(fname+'\n')
# #ax1 = fig.add_subplot(2,2,i+1)
# dattest[:,0] = np.random.randn(n) + 0.1*i
# dattest[:,1] = np.random.randn(n) + 0.1*i
# np.savetxt(fname,dattest)
#
#f.close()
#
#
#
#ccfc = 'ccfcomp.dat' 
#
#a = myscatterhistfile(fnam,colhist,nvert,nalong,[0,1],xlim = [-5,5], ylim=[-5,5])
#
#
#plt.savefig('figcorplot.png')
#
#