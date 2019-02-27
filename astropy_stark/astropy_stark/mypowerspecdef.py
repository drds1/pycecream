import numpy as np
#t,x,sig,f,np.array(ampspec),np.array(sigampspec),np.array(phasespec),np.array(sigphasespec)
# outputs t,x,sig data and errors
# ;       f, ampspec,sigampspec,phasespec,sigphasespec # frequencies, amplitude spectrum and errors and phase

proprtdtmin = 0.02

def mypowerspec(fname):    
    

    dat =np.loadtxt(fname)
    t   =dat[:,0]
    thi =t.max()
    tlo = t.min()
    tlen = thi - tlo
    x   =dat[:,1]
    sig =dat[:,2]
#    sig = sig*0 + 0.1
    nx  = x.shape[0]
    nxmod = nx*10
    dtmod = (t[-1]-t[0])/(nxmod-1)
    tmod=np.arange(t[0],t[-1]+dtmod, dtmod)
    nxmod=tmod.shape[0]
    
    xop = np.zeros(nxmod)
    sop=[]
    sigsop=[]
    cop=[]
    sigcop=[]
    ampspec=[]
    sigampspec=[]
    phasespec=[]
    sigphasespec=[]
    
#specify the lower and upper frequency limits
    dt=t[1:] - t[:-1]
    dt=dt[np.nonzero(dt)] #problem fixed if 2 datapoints taken at same time
    dtmin = dt.min() #this is the interval between points

# this fixes data points that are super close together making the nyquist frequency higher than can be computed in a reasonable time.    
    if (dtmin < properdtmin):
        print 'Time resolution too high, using average time resolution...'
        dtmin = (tlen/(nx-1))
        if (dtmin < properdtmin):
            print 'The time resolution on this lightcurve is better than sense! I dont like it and will use',properdtmin,'instead.'
            dtmin = properdtmin
    


    nf  = np.int( np.ceil(0.5*nx) )
    fhi = 0.5/dtmin              
    flo = 1/(t[-1]-t[0])
    df = 0.5*(fhi - flo)/(nf-1)
    f = np.arange(flo,fhi+df,df)
    nf = f.shape[0]
    
    w = 2*np.pi*f
# subtract the mean from the data (to avoid worry about the constant offset
    xmean = np.median(x)
    xinp = x - xmean
    xop = np.ones(nxmod)*xmean
    
    
    for i_f in range(nf):
        optscalop = myoptscal(t,xinp,sig,f[i_f])
        s    = optscalop[0]
        sigs = optscalop[1]
        c    = optscalop[2]
        sigc = optscalop[3]
        sop.append(s)
        sigsop.append(sigs)
        cop.append(c)
        sigcop.append(sigc)
        

             
#these bits are updated every i_f
        fit = s*np.sin(w[i_f]*t) + c*np.cos(w[i_f]*t)
        xinp = xinp - fit
        fitmod = s*np.sin(w[i_f]*tmod) + c*np.cos(w[i_f]*tmod)
        
#        print fitmod.mean(),i_f,'test',np.std(fitmod)
# if you fit to too high frequencies, the model goes to crazy high values, introduce a check for this below.
        if (np.abs( np.std(fitmod) ) > 100*xmean):
            print 'The fit has gone crazy! you have tried to fit to too high a frequency... stopping here.'
            nf = len(ampspec)
            break

# amplitude spectrum        
        sigssq = 2*s*sigs
        sigcsq = 2*c*sigc
        a =s*s+c*c
        siga = np.sqrt(sigssq*sigssq + sigcsq*sigcsq)
        ampspec.append(np.sqrt(a))
        sigampspec.append( 0.5/np.sqrt(a)* siga )

# phase spectrum
        soverc    = s/c
        sigsoverc = soverc * np.sqrt( (sigs/s)**2 + (sigc/c)**2 )
        phasespec.append( np.arctan(-1.*soverc) )
        sigphasespec.append( np.abs( sigsoverc/(1. + soverc*soverc) ) ) 

# update the plot model
        xop = xop + fitmod
        


# fit a broken power law to the amplitude spectrum
    ampspeclog = np.log10(np.array(ampspec))
    sigampspeclog = np.abs( np.array(sigampspec)/(np.array(ampspec)*loge10 ) )
   
    return(t,x,sig,f,np.array(ampspec),np.array(sigampspec),np.array(phasespec),np.array(sigphasespec)) 