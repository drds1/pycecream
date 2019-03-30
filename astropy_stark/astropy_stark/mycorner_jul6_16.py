import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

def mce(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip
    
    
    
    
def mc_b(cov, pos, nstd=2, ax=None, **kwargs):

 return mce(cov, pos, nstd=nstd, ax=ax)


import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
import matplotlib.cm as cm
import scipy
import scipy.stats
try:
 import myconfidence_ellipse as mc_b
except:
 pass

def corner_1(xs, weights=None, labels=None, extents=None, truths=None,
           truth_color="#4682b4", scale_hist=False, quantiles=[],
           verbose=True, plot_contours=False, plot_datapoints=False,
           fig_ax=None,header=[],sigconts=[100.-68,100.-95,100-99.7],figname='histplot.pdf',skip=1,xti=[-1],yti=[-1], annotate=[],
           medplot=1,sigplot=0.68,sigmedann=1,sigmedann_pre_in=[],title='',sigmedann_post_in=[],xtxtcoord=0.95,ytxtcoord=0.95,
           specax=[],specaxlab=[],plot_corr=0,
           lbm_in = 0.65,ltr_in = 0.4,xlab_coord=-0.3,ylab_coord=-0.3,**kwargs):
    
    npar = np.shape(xs)[1]
    if (len(specax) == 0):
     specaxin = [[] for i in range(npar)]
    else:
     specaxin = list(specax)
    if (len(specaxlab) == 0):
     specaxlabin = [[] for i in range(npar)]
    else:
     specaxlabin = list(specaxlab)
    
    
    if (len(header) ==0):
     header=['']*npar
    
    """
    Make a *sick* corner plot showing the projections of a data set in a
    multi-dimensional space. kwargs are passed to hist2d() or used for
    `matplotlib` styling.

    Parameters
    ----------
    xs : array_like (nsamples, ndim)
        The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space.

    weights : array_like (nsamples,)
        The weight of each sample. If `None` (default), samples are given
        equal weight.

    labels : iterable (ndim,) (optional)
        A list of names for the dimensions.

    extents : iterable (ndim,) (optional)
        A list where each element is either a length 2 tuple containing
        lower and upper bounds (extents) or a float in range (0., 1.)
        giving the fraction of samples to include in bounds, e.g.,
        [(0.,10.), (1.,5), 0.999, etc.].
        If a fraction, the bounds are chosen to be equal-tailed.
        OR (1,1) if same number is chosen then min and max values are used as limits

    truths : iterable (ndim,) (optional)
        A list of reference values to indicate on the plots.

    truth_color : str (optional)
        A ``matplotlib`` style color for the ``truths`` makers.

    scale_hist : bool (optional)
        Should the 1-D histograms be scaled in such a way that the zero line
        is visible?

    quantiles : iterable (optional)
        A list of fractional quantiles to show on the 1-D histograms as
        vertical dashed lines.

    verbose : bool (optional)
        If true, print the values of the computed quantiles.

    plot_contours : bool (optional)
        Draw contours for dense regions of the plot.

    plot_datapoints : bool (optional)
        Draw the individual data points.
    
    annotate : list[...ni] of strings to add to each plot (each histogram plot) NOTE annotate and sigmedann put things in the same place
    
    
    sigconts : plot contour levels by default at 1, 2, 3 sigma... e.g [100.-68,100.-95,100-99.7]. set sigconts[0] = 1 to fit n ellipse to you probability distribution and fit contours that way (must set sigconts[0] = 1 to do this where 1 is the 1 sigma elipse fitted).
    
    fig_ax : matplotlib.Figure (optional) and axes [fig,axes] format
        Overplot onto the provided figure object.
    
    plot_corr : 1 to annotate plots with correlation coefficient
    
    medplot : 1 to annotate histograms with median line
    
    sigplot : >0 to annotate histograms with sigplot confidence intervals (e.g 0.68) for 1 sig confidence upper lower limits

    sigmedann : 1 to label the histogram plots with uncertainties and median 
    
    sigmedann_pre_in [ndim] n-dimensional list of strings to put infront of the uncertainty and median parameter estmates if annotating
    
    sigmedann_post_in [ndim] n-dimensional list of strings to put at end of the uncertainty and median parameter estmates if annotating
    
    lbm_in: size of left and bottom margains
    
    ltr_in: size of top and right margains
    
    xlab_coord: distance of x-axis label from outside edge of axis -ve means plot axis label outside of frame
    
    ylab_coord: distance of y-axis label from outside edge of axis -ve means plot axis label outside of frame
    """
    nconts = len(sigconts)
    
    #deal with 1d parameter arrays
    ndtemp = np.shape(xs)[1]
    if (ndtemp == 1):
     xs = xs[:,0]
    
    
    #deal with annotations of uncertainties
    if (sigmedann_pre_in == []):
     sigmedann_pre = ['']*ndtemp
    else:
     sigmedann_pre = sigmedann_pre_in
    if (sigmedann_post_in == []):
     sigmedann_post = ['']*ndtemp
    else:
     sigmedann_post = sigmedann_post_in
 



    
    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
        
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    #assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
    #                                   "dimensions than samples!"

    if weights is not None:
        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError('weights must be 1-D')
        if xs.shape[1] != weights.shape[0]:
            raise ValueError('lengths of weights must match number of samples')

    # backwards-compatibility
    plot_contours = kwargs.get("smooth", plot_contours)
    
    
    
    K = len(xs)
    
    #print(K)
    factor = 1.0*6/K           # size of one side of one panel
    lbdim = lbm_in * factor   # size of left/bottom margin lbm_in = 0.65
    trdim = ltr_in * factor  # size of top/right margin #ltr_in = 0.4
    whspace = 0.05         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    if fig_ax is None:
        fig, axes = pl.subplots(K, K, figsize=(dim, dim))
    else:
        fig = fig_ax[0]
        axes = fig_ax[1]
    #    try:
    #        axes = np.array(fig.axes).reshape((K, K))
    #    except:
    #        raise ValueError("Provided figure has {0} axes, but data has "
    #                         "dimensions K={1}".format(len(fig.axes), K))
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    if extents is None:
        extents = [[x.min(), x.max()] for x in xs]

        # Check for parameters that never change.
        m = np.array([e[0] == e[1] for e in extents], dtype=bool)
        if np.any(m):
            raise ValueError(("It looks like the parameter(s) in column(s) "
                              "{0} have no dynamic range. Please provide an "
                              "`extent` argument.")
                             .format(", ".join(map("{0}".format,
                                                   np.arange(len(m))[m]))))
    else:
        # If any of the extents are percentiles, convert them to ranges.
        for i in range(len(extents)):
            try:
                emin, emax = extents[i]
                #if only want to manually set limits for one dimension
                if (emin == emax):
                 extens[i] =(x[i].min(),x[i].max())
            except TypeError:
                q = [0.5 - 0.5*extents[i], 0.5 + 0.5*extents[i]]
                extents[i] = quantile(xs[i], q, weights=weights)
    
    
    
      
    
    
    
    
    
    for i, x in enumerate(xs):
        if (ndtemp == 1):
            fig = pl.figure()
            ax  = fig.add_subplot(111)
        else:
            ax = axes[i, i]
        # Plot the histograms.
        n, b, p = ax.hist(x, weights=weights, bins=kwargs.get("bins", 50),
                          range=extents[i], histtype="step",normed=1,
                          color=kwargs.get("color", "k"))
        
        
        
        if (len(annotate) > 0):
         ax.text(xtxtcoord,ytxtcoord,annotate[i],ha='right',transform=ax.transAxes)
        
        ax.set_title(header[i],loc='left')
        xlimlo = extents[i][0]
        xlimhi = extents[i][1]
        binwidth_x = (xlimhi-xlimlo)/kwargs.get("bins", 50)
        ylhist = list(ax.get_ylim())
        ylhist[1] = ylhist[1]*1.1
        ax.set_ylim(ylhist)
        
        
        
        
        
        
        if (xti == []):
         ax.set_xticklabels([])
        if (yti == []):
         ax.set_yticklabels([])
        
        if truths is not None:
            if (truths[i] != ''):
             ax.axvline(truths[i], color=truth_color)
 
        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            qvalues = quantile(x, quantiles, weights=weights)
            for q in qvalues:
                ax.axvline(q, ls="dashed", color=kwargs.get("color", "k"))
            if verbose:
                print("Quantiles:")
                print(zip(quantiles, qvalues))

        # Set up the axes.
        ax.set_xlim(extents[i])
        if scale_hist:
            maxn = np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(5))

        
        
        #draw sigplot and medplot lines if optioned toggled on
        if ((medplot == 1 )or (sigplot > 0) or (sigmedann == 1)):
         ylnow   = list(ax.get_ylim())
         idxsort = np.argsort(x)
         nidx    = np.shape(idxsort)[0]
         xsort   = x[idxsort]
         
        if (medplot == 1 or sigmedann == 1):
         idmed = nidx/2
         xmed = xsort[idmed]
         if (medplot == 1):
          ax.plot([xmed,xmed],ylnow,color='k')
        if ((sigplot > 0) or (sigmedann == 1)):
         spp = np.abs(sigplot)
         if (np.abs(sigplot) >1):
          print('can only plot confidence region between 0 and 1. change sigplot lt 1')
          print('make it -1 gt sigplot lt 0 if dont want to plot uncertainties but still want to annotate them')
          spp = 0.68
         
         frac  = (1.-spp)/2
         idxlo = int(frac*nidx)
         idxhi = int((1.-frac)*nidx)
         xlo   = xsort[idxlo]
         siglo = xmed - xlo
         xhi   = xsort[idxhi]
         sighi = xhi  - xmed
         if (sigplot > 0):
          ax.plot([xlo,xlo],ylnow,ls='--',color='k')
          ax.plot([xhi,xhi],ylnow,ls='--',color='k')
        
        
        if ( sigmedann == 1):
         txt = r' '+sigmedann_pre[i]+' $='+str(np.round(xmed,2))+'^{+'+str(np.round(sighi,2))+'}_{-'+str(np.round(siglo,2))+'}$'+sigmedann_post[i]
         ax.text(xtxtcoord,ytxtcoord,txt,ha='right',transform=ax.transAxes)
        
        
         
        # Not so DRY.
        if i < K - 1:
            ax.set_xticklabels([])
        else:
            
            #change x axis labeling if custom option entered
            if (specaxlabin[i] != []):
             ax.set_xticklabels(specaxlabin[i])
             print('xtick labels set',specaxlabin[i])
             print(ax.get_xlim())
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[i])
                ax.xaxis.set_label_coords(0.5, xlab_coord)

        
        
        #print 'FUCK OFFFFF dont get it',i,specaxin[i],ndtemp
        if (specaxin[i] != []):
             ax.set_xticks(specaxin[i])
             #print 'xticks set', specaxin[i]
        #need go no further if only varying one parameter
        if (ndtemp == 1):
         [l.set_rotation(0) for l in ax.get_xticklabels()]
         print(i,labels,'herere')
         if labels is not None:
          ax.text(0.5,-0.07,labels[i],va='top',ha='center',transform=ax.transAxes)
         if (figname != ''):
          #pl.tight_layout()
          pl.savefig(figname)
         return fig, ax
        
        for j, y in enumerate(xs):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue
            elif j == i:
                continue

            hist2d(y[::skip], x[::skip], ax=ax, extent=[extents[j], extents[i]],
                   plot_contours=plot_contours,
                   plot_datapoints=plot_datapoints,
                   weights=weights, **kwargs)
            ylimlo = extents[j][0]
            ylimhi = extents[j][1]
            binwidth_y = (ylimhi-ylimlo)/kwargs.get("bins", 50)

 
            if (sigconts[0] == 1):
             cov = np.cov(y,x)
             mean = np.mean(y),np.mean(x)
             for scnow in sigconts:
              mc_b.mce(cov, mean, nstd=scnow,ax=ax,fill=False,edgecolor=kwargs.get("color", "k"),linewidth=2.0,color=kwargs.get("color", "k"))
            elif (sigconts[0] == 0):
             print('')  
            else:
             xtemp = y
             ytemp = x
             xlotemp = ylimlo
             xhitemp = ylimhi
             ylotemp = xlimlo
             yhitemp = xlimhi
             bwx     = binwidth_y
             bwy     = binwidth_x
             
             both = np.array((xtemp,ytemp)).T
             #both = np.array((x,y))
             #print('making pdf')
             pdf1=scipy.stats.kde.gaussian_kde(both.T)
             #print(ylimlo,ylimhi,binwidth_y,xlimlo,xlimhi,binwidth_x,extents)
             #q,w=np.meshgrid(np.arange(ylimlo,ylimhi,binwidth_y), np.arange(xlimlo,xlimhi,binwidth_x))
             q,w=np.meshgrid(np.arange(xlotemp,xhitemp,bwx), np.arange(ylotemp,yhitemp,bwy))
             r1=pdf1([q.flatten(),w.flatten()])
             #r1 = r1/np.max(r1)
             s1=[]
             #print('doing contours',nconts,sigconts,xlimlo,xlimhi,binwidth_x,ylimlo,ylimhi,binwidth_y,np.min(r1),np.max(r1))
             #print(pdf1)
             #print('pdf1')
             #raw_input()
             #print(q)
             #print('q')
             #raw_input()
             #print(w)
             #print('w')
             #raw_input()
             #print(r1)
             #print('r1')
             #raw_input()
             for ic in range(nconts):
              s1.append(scipy.stats.scoreatpercentile(pdf1(pdf1.resample(1000)), sigconts[ic]))
             r1.shape=(q.shape[0],q.shape[1])
             for ic in range(nconts):
              #ax.contour(np.arange(xlimlo,xlimhi,binwidth_x), np.arange(ylimlo,ylimhi,binwidth_y), r1, [s1[ic]])
              #print('contourminmax',xlimlo,xlimhi,ylimlo,ylimhi,np.min(r1),np.max(r1))
              #ax.contour(np.arange(xlimlo,xlimhi,binwidth_x), np.arange(ylimlo,ylimhi,binwidth_y), r1, [np.max(r1)/2],color='r',lwd=20)
              ax.contour(np.arange(xlotemp,xhitemp,bwx), np.arange(ylotemp,yhitemp,bwy), r1,[s1[ic]])
              #ax.text(0.5,0.5,'erewrw')
              #pl.show()
             
             
             
             
            if (xti == []):
             ax.set_xticklabels([])
            if (yti == []):
             ax.set_yticklabels([])
            
            if (plot_corr > 0):
             r_ann = np.corrcoef(np.array((y[::skip],x[::skip])))[0,1]
             nt = np.shape(y)[0]
             #for i in range(nt):
             # print(y[i],x[i])
             #print(r_ann)
             #print(np.corrcoef(np.array((y[1:],x[1:]))))
             #raw_input()
             ax.text(xtxtcoord,ytxtcoord-0.05,'R='+str(np.round(r_ann,2)),ha='right',transform=ax.transAxes)
             #ax.text(0.95,0.9,'R='+str(np.round(r_ann,2)),ha='right',transform=ax.transAxes)
            
            #if (len(annotate) > 0):
            # ax.text(0.95,0.9,annotate[i],ha='right',transform=ax.transAxes)
            
            if truths is not None:
                
                if (truths[i] != '' and truths[j] != '' ):
                 ax.plot(truths[j], truths[i], "s", color=truth_color)
                
                if (truths[j] != ''):
                 ax.axvline(truths[j], color=truth_color)
                
                if (truths[i] != ''):
                 ax.axhline(truths[i], color=truth_color)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
           
           
           
           
           
           
           
           
           
           
           
           
           
           
            #print 'FUCK OFFFFF dont get it',i,specaxin[i],ndtemp
            if (specaxin[j] != []):
             ax.set_xticks(specaxin[j])
            if (specaxin[i] != []):
             ax.set_yticks(specaxin[i])

             #print 'xticks set', specaxin[i]



            if i < K - 1:
                ax.set_xticklabels([])
            else:
                if (specaxlabin[j] != []):
                 ax.set_xticklabels(specaxlabin[j])
                 #print 'xtick labels set',specaxlabin[i]
                 #print ax.get_xlim()

                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j])
                    ax.xaxis.set_label_coords(0.5, xlab_coord)

            if j > 0:
                ax.set_yticklabels([])
            else:
                if (specaxlabin[i] != []):
                 ax.set_yticklabels(specaxlabin[i])
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i])
                    ax.yaxis.set_label_coords(ylab_coord, 0.5)

    
    if (title != ''):
     ax.set_title(title,loc='left')
     
    if (figname != ''):
     pl.savefig(figname)
    return fig,axes


def quantile(x, q, weights=None):
    """
    Like numpy.percentile, but:

    * Values of q are quantiles [0., 1.] rather than percentiles [0., 100.]
    * scalar q not supported (q must be iterable)
    * optional weights on x

    """
    if weights is None:
        return np.percentile(x, [100. * qi for qi in q])
    else:
        idx = np.argsort(x)
        xsorted = x[idx]
        cdf = np.add.accumulate(weights[idx])
        cdf /= cdf[-1]
        return np.interp(q, cdf, xsorted).tolist()


def error_ellipse(mu, cov, ax=None, factor=1.0, **kwargs):
    """
    Plot the error ellipse at a point given its covariance matrix.

    """
    # some sane defaults
    facecolor = kwargs.pop('facecolor', 'none')
    edgecolor = kwargs.pop('edgecolor', 'k')

    x, y = mu
    U, S, V = np.linalg.svd(cov)
    theta = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
    ellipsePlot = Ellipse(xy=[x, y],
                          width=2 * np.sqrt(S[0]) * factor,
                          height=2 * np.sqrt(S[1]) * factor,
                          angle=theta,
                          facecolor=facecolor, edgecolor=edgecolor, **kwargs)

    if ax is None:
        ax = pl.gca()
    ax.add_patch(ellipsePlot)

    return ellipsePlot


def hist2d(x, y, *args, **kwargs):
    """
    Plot a 2-D histogram of samples.

    """
    ax = kwargs.pop("ax", pl.gca())

    extent = kwargs.pop("extent", [[x.min(), x.max()], [y.min(), y.max()]])
    bins = kwargs.pop("bins", 50)
    color = kwargs.pop("color", "k")
    linewidths = kwargs.pop("linewidths", None)
    plot_datapoints = kwargs.get("plot_datapoints", True)
    plot_contours = kwargs.get("plot_contours", True)

    cmap = cm.get_cmap("gray")
    cmap._init()
    cmap._lut[:-3, :-1] = 0.
    cmap._lut[:-3, -1] = np.linspace(1, 0, cmap.N)

    X = np.linspace(extent[0][0], extent[0][1], bins + 1)
    Y = np.linspace(extent[1][0], extent[1][1], bins + 1)
    try:
        H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=(X, Y),
                                 weights=kwargs.get('weights', None))
    except ValueError:
        raise ValueError("It looks like at least one of your sample columns "
                         "have no dynamic range. You could try using the "
                         "`extent` argument.")

    V = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]

    for i, v0 in enumerate(V):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]

    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
    X, Y = X[:-1], Y[:-1]

    if plot_datapoints:
        ax.plot(x, y, "o", color=color, ms=1.5, zorder=-1, alpha=0.1,
                rasterized=True)
        if plot_contours:
            ax.contourf(X1, Y1, H.T, [V[-1], H.max()],
                        cmap=LinearSegmentedColormap.from_list("cmap",
                                                               ([1] * 3,
                                                                [1] * 3),
                        N=2), antialiased=False)

    if plot_contours:
        ax.pcolor(X, Y, H.max() - H.T, cmap=cmap)
        ax.contour(X1, Y1, H.T, V, colors=color, linewidths=linewidths)

    data = np.vstack([x, y])
    mu = np.mean(data, axis=1)
    cov = np.cov(data)
    if kwargs.pop("plot_ellipse", False):
        error_ellipse(mu, cov, ax=ax, edgecolor="r", ls="dashed")

    ax.set_xlim(extent[0])
    ax.set_ylim(extent[1])
    







def hist_1(xs, labels=None, lims=None, truths=None,nbins=[50],colin ='b',truth_color = 'r',annotate = None, figname = 'histplot.pdf'):
 
 fig = pl.figure()
 ax1 = fig.add_subplot(111)
 
 if (labels != None):
  ax1.set_xlabel(labels[0])
  
 
 if (lims != None):
  ax1.set_xlim(lims[0])
 
 
 if (len(nbins) > 0):
  bins = np.array(nbins)
 else:
  bins = nbins[0]
 
 if (truths != None):
  
  ylimnew = ax1.get_ylim
  ax1.vlines(truths[0],ylimnew[0],ylimnew[1],col = truth_color)
  
 ax1.hist(xs,bins = nbins[0], histtype='step',normed=1,color = colin)
 
 if (annotate != None):
  ax1.text(xtxtcoord,ytxtcoord-0.05,annotate,ha='right',transform=ax1.transAxes)
 
 pl.savefig(figname)
 
 return()
  
 
 
def corner(xs, weights=None, labels=None, extents=None, truths=None,truth_color="#4682b4", scale_hist=False, quantiles=[],verbose=True, plot_contours=False, plot_datapoints=False,fig=None,sigconts=[100.-68,100.-95,100-99.7],figname = 'histplot.pdf',annotate=[], **kwargs):


 
 ndim = np.shape(xs)
 
 if (ndim[1] == 1):
  if (len(annotate) > 0):
   text = annotate[0][0]
  hist_1(xs, labels=labels, truths=truths,colin = 'k',figname = figname, truth_color = truth_color, annotate = text)
  
 else:
  corner_1(xs, weights=weights, labels=labels, extents=extents, truths=truths,
           truth_color=truth_color, scale_hist=scale_hist, quantiles=quantiles,
           verbose=verbose, plot_contours=plot_contours, plot_datapoints=plot_datapoints,
           fig=fig,sigconts=sigconts,figname = figname, annotate = annotate, **kwargs)

 return()
  