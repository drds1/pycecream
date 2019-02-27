# python code to produce an optimal average
#is spec = 1 then avoid points with zero error bars
import numpy as np


def myoptaverage(yin,sigyin,spec=1):
    
    ny = np.shape(yin)[0]
    if (spec == 1):
     idx = np.where(sigyin != 0)[0]
    else:
     idx = np.arange(ny)
    
    y = yin[idx]
    sigy = sigyin[idx]
    
    sigy2=sigy*sigy
    
    top=(y/sigy2).sum()
    bot=(1/sigy2).sum()
    optav=top/bot
    
    sig = np.sqrt(1/bot)
    return(optav,sig)
    
