import corner
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns; sns.set(style="ticks", color_codes=True)
import pandas as pd

#Set up the parameters of the problem.
ndim, nsamples = 3, 50000

# Generate some fake data.
np.random.seed(42)
data1 = np.random.randn(ndim * 4 * nsamples // 5).reshape([4 * nsamples // 5, ndim])
data2 = (4*np.random.rand(ndim)[None, :] + np.random.randn(ndim * nsamples // 5).reshape([nsamples // 5, ndim]))
data = np.vstack([data1, data2])

# Plot it.
figure = corner.corner(data, labels=[r"$x$", r"$y$", r"$\log \alpha$", r"$\Gamma \, [\mathrm{parsec}]$"],
                       quantiles=[0.45, 0.5, 0.55],
                       show_titles=True, title_kwargs={"fontsize": 12})
idx = 0
for ax in figure.axes:
    ax.annotate('axes fraction '+str(idx),xy = (0.9,0.95),
                xycoords='data',
            xytext=(0.9, 0.95), textcoords='axes fraction ',
            horizontalalignment='right', verticalalignment='top')
    idx += 1

#plt.show()


def corrfunc(x,y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    #r, _ = pearsonr(x, y)
    r = np.corrcoef(x,y)[0, 1]
    ax = ax or plt.gca()
    # Unicode for lowercase rho (ρ)
    rho = '\u03C1'
    ax.annotate(f'{rho} = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)



g = sns.pairplot(pd.DataFrame(data,columns = ['col 1','col 2','col 3']), corner=True)
g.map_lower(corrfunc)
