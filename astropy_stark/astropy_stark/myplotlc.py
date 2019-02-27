#new sep4 different colours for each legens element unique colour python loop colour 
#plot all light curves with file names listed in fname
#if seperate == 1, plot on seperate graphs, else plot all on same
# if xminforce,xmaxforce entered, force the time limits of the plot

import numpy as np
#import pylab as plt
import os
import matplotlib.pylab as plt
import matplotlib.ticker as mtick

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.weight'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
def myplotlc(fname,seperate,xminforce = -1,xmaxforce = -1, wavfile = '',wav_abv_name='', savefigname='',axlab = ['',''],showticks=[1,1],plottype='e',multiplot=0,fig=0,forcecol=[],ann=[],hlin=[],vlin=[],conv2ergs=0,alpha=1.0,linewidth=1.0,ylimforce=[],yt=[-1],fiddle=[]):

 xmintemp = []
 xmaxtemp = []
 nann = len(ann)
 
 with file(fname) as f:
  lines = f.read().splitlines()


 datplot = []
 nf = len(lines)

 include_wavs = 0
 include_wav_abv = 0
 f.close()
 
# deal with wavelength files (if included) and abbreviations (if included) 
 if (len(wavfile) > 0):
  include_wavs = 1
  wav = np.array(np.loadtxt(wavfile))
  a = np.zeros(nf)
  a[:] = 1.*wav
  wav = 1.*a
  
  wav_scaled = (wav - wav.min()) / wav.ptp()
  wav_col    = plt.cm.nipy_spectral(wav_scaled)
 else:#plot according to wavelength or (if not entered) blue (=0.2 with nipy_spectral cmap)
  wav = [1000]*nf
  wav_col    = plt.cm.nipy_spectral(np.linspace(0,1,nf))
 
 if (len(wav_abv_name) > 0):
  f = open(wav_abv_name)
  wav_abv = f.readlines()
  f.close()
  include_wav_abv =1 
 else:
  wav_abv = ['']*nf
 
 
 if (len(forcecol) != 0):
  wav_col = forcecol
 
 
 
 plotlist=[] #list holding ax1.plot command for legend
 for i in range(nf):
  print 'Loading file ',i,' of', nf,'... ',lines[i].strip()
  
  
  print 'argh!!!', lines[i]
  print os.getcwd()
  print ''
  dat = np.loadtxt(lines[i].strip())
  
  if (len(fiddle) > 0):
   dat[:,1] = dat[:,1]*fiddle[i]
   dat[:,2] = dat[:,2]*fiddle[i]
  
  if ((conv2ergs == 1) and (wav[i] > 1.) ):
   #print i
   #print wav
   #print wav[i]
   #raw_input()
   #b = np.float(wav[i])
   #print b
   #raw_input()
   a = 2.9979e-8 *1.e16/ np.float(wav[i])**2
   dat[:,1] = dat[:,1]*a
   dat[:,2] = dat[:,2]*a
  
  
  print 'minimum spacing between points...',np.min(dat[1:,0]-dat[:-1,0])
  print ''
  #print 'diagnose ', len(dat)
  xmintemp.append(np.min(dat[:,0]))
  xmaxtemp.append(np.max(dat[:,0]))
  datplot.append(dat)


#sort the plot limits
 if (xminforce < 0):
  xmin = min(xmintemp)
 else:
  xmin = xminforce

 if (xmaxforce < 0):
  xmax = max(xmaxtemp)
 else:
  xmax = xmaxforce
#
  
 if (multiplot ==0):
  fig = plt.figure()
  
 if (seperate == 0):
  ax1 = fig.add_subplot(1,1,1) 
 
 nhlin = len(hlin) 
 nvlin = len(vlin) 
 l=[]
 for i in range(nf):
  if (seperate ==1):
   ax1 = fig.add_subplot(nf,1,i+1)

## add title
  if (nann > 0):
   ax1.text(ann[i][1],ann[i][2],ann[i][0],ha='left',transform=ax1.transAxes,color = wav_col[i])


## add horrizontal lines
  if (nhlin > 0):
   hlintemp = hlin[i]
   ax1.hlines(hlintemp,xmin = xmin,xmax=xmax,color = wav_col[i])
    
    
## add vertical lines
  if (nvlin > 0):
   vlintemp = vlin[i]
   for xline in vlintemp:
    ab=ax1.get_ylim()
    ax1.vlines(xline,ymin = ab[0],ymax=ab[1],color = wav_col[i])    
#label with wavelengths and abbreviations if present
  
  #if (include_wavs == 1):
  # ax1.text(0.9,0.92,str(int(round(wav[i],0)))+"$\AA$",ha='left',va='top',transform=ax1.transAxes, fontsize = 12)
  #if (include_wav_abv ==1):
  # ax1.text(0.9,0.42,wav_abv[i],ha='left',va='top',transform=ax1.transAxes,fontsize = 12)
  
  if (wav_abv[i] != ''):
   ltemp = wav_abv[i]+"$\AA$"
  else: 
   ltemp = str(int(round(wav[i],0)))+"$\AA$"
  
  
  if (plottype == 'e'):
   plotlist.append( ax1.errorbar(datplot[i][:,0], datplot[i][:,1], datplot[i][:,2],linestyle='None',color = wav_col[i],label = ltemp,alpha=alpha,linewidth=linewidth) )
  elif (plottype == 'l'):
   plotlist.append( ax1.plot(datplot[i][:,0], datplot[i][:,1],marker = 'None', color = wav_col[i],label = ltemp,alpha=alpha,linewidth=linewidth) )
  elif (plottype == 'p'):
   plotlist.append( ax1.scatter(datplot[i][:,0], datplot[i][:,1], color = wav_col[i],marker='o', edgecolors='None', label = ltemp,alpha=alpha,linewidth=linewidth) )
   #ax1.plot(datplot[i][:,0], datplot[i][:,1],linestyle = 'None', color = wav_col[i],marker='o', edgecolors='None')

  #only display xticks on bottom
  #l.append(ltemp)
  
  
  
  ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
  #if ((i < nf - 1) or (showticks[0] == 0) ):
  # plt.setp(ax1.get_xticklabels(), visible=False)
  
  if (showticks[1] == 0):
   plt.setp(ax1.get_yticklabels(), visible=False)

  if (i == nf-1):
   ax1.set_xlabel(axlab[0])
  if (i == nf/2):
   if (axlab[1] == ''):
    ax1.set_ylabel(r'$f_\nu ( \lambda )  / 10^{-14} \left( erg s^{-1} cm^{-2}  \AA^{-1} \right)$')
   else:
    ax1.set_ylabel(axlab[1]) 
  #if (i > nf -2):
  # ax1.tick_params(axis='x', pad = 30)
   
  
  ax1.set_xlim(xmin,xmax)
  #ax1.set_xticks(np.linspace(xmin,xmax,4))
  if (seperate ==1):
   ylims = ax1.get_ylim()
   yrange = ylims[1]-ylims[0]
   ax1.set_ylim(ylims[0]-yrange/20,ylims[1]+yrange/20)
   ax1.tick_params(axis='y',labelsize=16)
 
 leg = ax1.legend(numpoints = 1)
 leg.draw_frame(False)
 
 for i in range(nf):
  text = leg.get_texts()[i]
  text.set_color(plotlist[i][0].get_color())

  
 if (len(ylimforce) > 0 ):
  ax1.set_ylim((ylimforce[0],ylimforce[1]))
 
 if (len(yt) > 0):
  if (yt[0] != -1):
   ax1.set_yticks(yt)  
    
 
 #print xmintemp
 #print xmaxtemp
 if (len(savefigname) == 0):
  plt.show()
 elif (savefigname == 'wait'):
  print 'waiting...'
 else:
  plt.tight_layout()
  plt.savefig(savefigname,bbox_inches='tight',orientation='portrait')
  

 return(datplot)
 
 
 
 
 
def saveindir(fname,dir,wavfile='',wav_abv_name=''): 
 seperate = 1
 
 pwd = os.getcwd()
 
 
 with file(fname) as f:
  lines = f.read().splitlines()
 f.close()
 
 nf = len(lines)
 wavfil=['']*nf
 wavabv=['']*nf
 
 if (wavfile != ''):
  with file(wavfile) as f:
   wavfil = f.read().splitlines()
  f.close()
 

 if (wav_abv_name != ''):
  with file(wav_abv_name) as f:
   wavabv = f.read().splitlines()
  f.close()
 
 

 
 for idx in range(nf): 
  wfroutine = ''
  wabvroutine=''
  savefigname='fig_'+str(idx)+'.png'
 



  if (wavfile != ''):
   f = open('wavfiletemp.dat','w')
   f.write(wavfil[idx])
   f.close()
   wfroutine = 'wavfiletemp.dat'
   savefigname = 'fig_'+str(wavfil[idx])+'.png'
 
 
  if (wav_abv_name != ''):
   f = open('wavabvtemp.dat','w')
   f.write(wavabv[idx])
   f.close()
   wabvroutine = 'wavabvtemp.dat'
   savefigname = 'fig_'+str(wavabv[idx])+'.png'
   

  f = open('filetemp.dat','w')
  f.write(lines[idx])
  print 'file ... ',lines[idx]
  f.close()
   
   
  
  if not os.path.exists(dir):
   os.makedirs(dir)
     
  #move temporary files to output directory
  #str_sys  = 'mv '+wfroutine+' '+dir  
  #os.system(str_sys)
  #str_sys  = 'mv '+wabvroutine+' '+dir  
  #os.system(str_sys)
  #str_sys  = 'mv filetemp.dat '+dir  
  #os.system(str_sys)
  
  
  #save the light curve in the temporary directory
  #os.chdir(dir)
  savefigname = dir+'/'+savefigname
  a = myplotlc('filetemp.dat',seperate,xminforce = -1,xmaxforce = -1, wavfile = wfroutine,wav_abv_name=wabvroutine, savefigname=savefigname,axlab = ['',''],showticks=[1,1],plottype='e',multiplot=0,fig=0,forcecol=[],ann=[],hlin=[],vlin=[])
 
  #clean up temporary files
  str_sys  = 'rm '+wfroutine  
  os.system(str_sys)
  str_sys  = 'rm '+wabvroutine 
  os.system(str_sys)
  str_sys  = 'rm filetemp.dat'  
  os.system(str_sys)
  
  #os.chdir(pwd)
  

 return()
  
  
  
