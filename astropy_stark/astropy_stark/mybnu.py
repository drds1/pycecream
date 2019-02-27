import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
plt.rcParams['axes.linewidth'] = 2.0 #set the value globally
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['axes.labelsize'] = 2.0



h    = 6.626e-34
c    = 2.9979e8
kb   = 1.38e-23
Mpc  = 3.0857e22
Mpc2 = Mpc*Mpc
ang  = 1.e-10
tday = 3600.*24
ld = c*tday





#calculate flux observed from a disk annulus at radius r T(r) = T0 (r/r0)^(-3/4)
#dr, r0 in light days
def myfnu(T,wavang,r_r0,T0, D_Mpc,degi = 0.0, dr = 0.1, r0 = 0.1):

 bnu_p = bnu(T,wavang)
 radi = np.pi/180. * degi
 cosi = np.cos(radi)
 
 d = D_Mpc*Mpc 
 a = 2*np.pi/(d*d) * dr * r_r0 *r0 * ld*ld
 op = bnu_p * a/(d*d)
 return(op)




#define plank function fnu
def bnu(T,wavang):
 
 wav = wavang*ang
 wav2 = wav*wav
 wav3 = wav2*wav
 a = 2*h*c/wav3
 
 b = 0
 if (T == 0):
  bnu = 0
 else:
  b = h*c/(wav*kb*T)
  bnu = a/ (np.exp(b) - 1.)
 #if (T != 0):
 # print 'bnu'
 # print bnu
 # print 'b'
 # print a
 # print 'c'
 # print b
 # print 'wavang'
 # print wavang
 # print 'temp'
 # print T
  #raw_input()
 return(bnu)
 
 
 
 

#make plots of planck function for a list of T[t1,t2,t3,..]
def plotbnu(T,wav,colplot=[],fname='plotbnu.png',ylo=0.0,ymax=[],fudge=[]):
 
 pwd=os.getcwd()
 nT = len(T)
 dir = './plotbnu_dir'
 
 if (fudge == []):
  fudge = [1.]*nT
 
 if (colplot==[]):
  colplot = ['b']*nT

 
#make directory to save bnu's in (remove old ones if present)
 if not os.path.exists(dir):
  os.makedirs(dir)
 
 os.chdir(dir)
 #os.system('rm *.png')
 
 
 fig = plt.figure()
 ax1 = fig.add_subplot(111)
 ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
 ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
 ax1.set_xlabel(r'$ \lambda / \AA $')
 ax1.set_ylabel(r'$f_{\nu} ( \lambda ) $')
 ax1.get_yaxis().set_ticks([])

 for i in range(nT):
  Tnow = T[i]
  specint = bnu(Tnow,wav)
  specint = specint/specint.max()*fudge[i]
  
  ax1.plot(wav,specint,color=colplot[i])
  
  #print ymax,'klljlksdjsl'
  if (ymax != []):
   ax1.set_ylim(ylo,ymax[i])
   print 'new ylim',i,Tnow,ymax[i]
  
  #fname = 'plotbnu_'+str(Tnow)+'.png'
  plt.savefig(fname)
  
  os.chdir(pwd)