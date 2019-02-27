import scipy.interpolate

def myextrapolate(xv,xarr,yarr):
    tck = scipy.interpolate.splrep(xarr, yarr, k=1, s=0)
    answer = scipy.interpolate.splev(xv, tck)
    return(answer)