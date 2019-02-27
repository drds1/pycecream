# python implementation of iterated optimal scaling
# we input a light curve (t,x,sig), a frequency f, and get out a sin and cosine amplitude and uncertainties for each
import numpy as np

#t=np.arange(1,100,0.5)
#n=t.shape[0]
#x=np.sin(0.5*t)
#sig=np.ones(n)*0.2
#f=0.5/2/np.pi

# function to fit the best s_k sin(w_k t) + c_k cos(w_k t) curve it can
# input f (frequecny NOT angular), t,x,sig times,fluxes,erros (1d arrays), 
# output s,sigs,c,sigc output sin and cosine amplitudes with uncertainties (all doubles)

def myoptscal(t,x,sig,f):
    #the bit above this is just a test to see if the thing below works. It will then be turned into a function
    sig2=sig*sig
    maxnits=100
    s=x.mean()
    sigs=0.1*s
    c=x.mean()
    sigc=0.1*c
    
    
    w=2*np.pi*f
    sw=np.sin(w*t)
    sw2=sw*sw
    cw=np.cos(w*t)
    cw2=cw*cw
    convfrac = 0.01 # 1 percent variation = convergence
    
    sconv = np.array([])
    cconv = np.array([])
    for idx in range(maxnits):
    #    print idx, s,sigs,c,sigc
    # deal with sin term    
        top = np.sum( (x-cw*c)*sw/sig2 )  
        bot = np.sum( sw2/sig2 )
    #    print top,bot,idx 
        s    = top/bot
        sigs = 1/bot
    
    # this is what happens if 0 freq is included    
        if ( (top == 0) and (bot == 0) ):
            s=0.0
            sigs=0.0
    # end of 0 freq problem 
        
    # Now deal with cosine term    
        top = np.sum( (x-sw*s)*cw/sig2 )
        bot = np.sum( cw2/sig2 )
        
        c    = top/bot
        sigc = 1/bot
        
        
    # now implement convergence check ########
        sconv = np.append(sconv,s)
        cconv = np.append(cconv,c)
        if (idx > 0):
            smean = sconv.mean()
            cmean = cconv.mean()
            sfrac = (s - smean)/smean
            cfrac = (c - cmean)/cmean
            
            if ( (np.abs(sfrac) < convfrac) and (np.abs(cfrac) < convfrac) ):
                break
                
            elif (idx == maxnits):
                print 'myoptscal.py warning: optimal scaling has not converged after',maxnits,' iterations!'
    # convergence check complete  ############          
    #    print s,sigs,c,sigc
    #    print ''
    #    print ''
    # print s,sigs,c,sigc
    # return(s,sigs,c,sigc)
    return(s,sigs,c,sigc)
        
        
        
    
    
    
    
    
    
    
