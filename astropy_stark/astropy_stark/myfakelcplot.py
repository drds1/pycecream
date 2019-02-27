import numpy as np
from pylab import *


def myfakelcplot(fname):
    from pylab import *
    dat=np.loadtxt(fname)
    errorbar(dat[:,0],dat[:,1],dat[:,2])

