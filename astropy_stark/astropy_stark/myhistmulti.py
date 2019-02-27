import numpy as np
import pylab as plt



#script to read an input file plot histograms of the data
#sublabfile if entered are the labels to label each of the plot windows
#xlim if entered are the limits to force all the histograms to lie within

def multihist(filename,nperhist,nvert,nalong,parcol,sublabfile = '',xlim=[], bin1=50):

 
 if (nperhist > 1):
  translucency = 0.6
 else:
  translucency = 0.0
 
 with file(filename) as f:
    fname = f.read().splitlines()
    
 nf  = len(fname)


 nsublab = len(sublabfile)
 if (nsublab > 0):
  with file(sublabfile) as f:
   sublab = f.read().splitlines()  
 


 
 dat=[]
 #!load the data
 for i in range(nf):
  print fname[i],i
  a   = np.genfromtxt(fname[i],skip_footer=1,usecols = np.arange(5))
  dat.append(a)
  #print a
  
#! plot the data
 fig = plt.figure()

#! alter python plotting to plot in column rather than row order 
 ic = 0
 isp = np.arange(nvert*nalong) + 1
 isp = np.reshape(isp,(-1,nalong))
 isp = np.reshape(isp,nvert*nalong,order='F')
 
 for i in range(nvert*nalong):
  ax1 = fig.add_subplot(nvert,nalong,isp[i])
    
  for i2 in range(nperhist):
   ax1.hist(dat[ic][:,parcol], edgecolor = 'None', alpha = translucency,bins = bin1)
   ic = ic + 1

 #label the subplots if file is provided
  if (nsublab > 1):
   
   a = np.array(ax1.axis()) #expand yaxis to make room for lab
   a[-1] = a[-1]*1.5
   ax1.axis(a)
   ax1.text(0.1,0.85,sublab[i],ha='left',va='top',transform=ax1.transAxes,fontsize = 12)
   
  #force the x axis limits to all be the same if entered
  if (len(xlim) > 0):
   a = np.array(ax1.axis())
   a[[0,1]] = xlim[0],xlim[1]
   ax1.axis(a)
   
  ax1.set_yticklabels('')
  if (np.mod(i+1,nvert) != 0):
   ax1.set_xticklabels('')



##script to test the above definition
#
#nf = 16
#nalong = 2
#nvert = 4
#parcol = 0
#nperhist = 2
#
#fnam = 'histmultidat.dat'
#f = open(fnam,'w')
#for i in range(nf):
# fname = 'test_multihis_'+str(i)+'.dat'
# f.write(fname+'\n')
# a = np.random.randn(1000)
# np.savetxt(fname,a)
#f.close()
#
#
#
#
#a = multihist(fnam,nperhist,nvert,nalong,parcol,sublabfile = fnam,xlim =[-5,5] )
#
#plt.savefig('myhistmultitest.png')
 

 
 
