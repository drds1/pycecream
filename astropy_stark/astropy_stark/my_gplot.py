import numpy as np
import matplotlib.pylab as plt
from inspect import isfunction

#if priorfunc == -1 then just use random walk prior,
#else priorfunc should be an array of the prior with one element for each frequency
def gp(p,sd,freq,p0=1.0,dw=1.0,w0=1.0,plot_tit='G_plot.pdf',xline=[],labann=[],
priorfunc = -1,square=1):


 cfur = p[0:-1:2]
 sd_cfur = sd[0:-1:2]
 sfur = p[1::2]
 sd_sfur = sd[1::2]
 fnow = freq
 
 
 
 if (type(priorfunc) == int):
  sig0_2 = 0.5*p0*dw*((w0/2/np.pi)/fnow)**2
 else:
  sig0_2 = 1.*priorfunc
 
 sig_c2 = sd_cfur**2
 sig_s2 = sd_sfur**2
 
 
 Gsk = ( sig0_2 / ( sig0_2 + sig_s2 ) - 0.5)*2
 Gck = ( sig0_2 / (sig0_2 + sig_c2) - 0.5)*2
 
 if (square == 0):
  Gsk = np.sqrt(Gsk.clip(min=0))
  Gck = np.sqrt(Gck.clip(min=0))
  
  #Gsk = np.sqrt(( sig0_2 / ( sig0_2 + sig_s2 ) - 0.5)*2)
  #Gck = np.sqrt(( sig0_2 / (sig0_2 + sig_c2) - 0.5)*2)
  #Gsk = ( np.sqrt(sig0_2) / ( np.sqrt(sig0_2) + np.sqrt(sig_s2) ) - 0.5)*2
  #Gck = ( np.sqrt(sig0_2) / (np.sqrt(sig0_2) + np.sqrt(sig_c2)) - 0.5)*2
 
 nfurtot = np.shape(Gck)[0]
 fig = plt.figure()
 ax2 = fig.add_subplot(111)
 ax2.plot(fnow,Gsk,label='Sk',color='r')
 ax2.plot(fnow,Gck,label='Ck',color='b')
 sum_gs = np.sum(Gsk)
 sum_gc = np.sum(Gck)
 sum_ave = (sum_gc + sum_gs)/2
 ax2.set_xlabel('frequency cycles/day')
 ax2.set_ylabel(r'$G= \frac{\sigma_0^2 - \sigma^2} {  \sigma_0^2 + \sigma^2 } $')
 ax2.set_title('G factor')
 ax2.text(0.98,0.55,r'$N_{\mathrm{eff}}=\sum_k \frac{G_{Sk} + G_{Ck}}{2} = '+np.str(np.int(sum_ave))+'$ of '+np.str(np.int(nfurtot)),ha='right',transform=ax2.transAxes,fontsize=14)
 ax2.set_xscale('log')
 ax2.set_yscale('log')
 ylim = list(ax2.get_ylim())
 for ian in range(len(xline)):
  ax2.plot([xline[ian]]*2,ylim,ls='--',color='k',label=labann[ian])
 
 plt.legend()
 plt.tight_layout()
 if (plot_tit == ''):
  plt.show()
 else:
  plt.savefig(plot_tit)
 
 return(sum_ave,Gsk,Gck)