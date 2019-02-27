# test_therm_plot scripr
import numpy as np
import pylab as plt
import scipy.optimize


function myplotxyfitline(fname,x,y,sig):




    dat=np.loadtxt(fname)

    x=dat[:,0]
    y=dat[:,1]
    n=y.shape[0]
    nmod=n*100


### analyse data
    def ymod(x,a,b):
        ym=a+b*x
        return(ym)
###


# analyse slope
    fit_coeff=scipy.optimize.curve_fit(ymod,x,y,sigma=sig)

    coeff_err=[fit_coeff[1][0,0],fit_coeff[1][1,1]]
    fit_coeff=fit_coeff[0]

    xmod=np.arange(x[0],x[-1],(x[-1]-x[0])/nmod)
    ymod=fit_coeff[0]+fit_coeff[1]*xmod
##



#plot log slope

    fig=plt.figure()
    ax1=fig.add_subplot(1,1,1)
    ax1.errorbar(x,y,sig)
    ax1.plot(xmod,ymod)
    plt.title('Something not very useful')

    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.text(0.5,0.5,'ic='+str(round(fit_coeff[0],2))+'+/-'+str(round(coeff_err[0]))+'  m='+str(round(fit_coeff[0],2))'+/-'+str(round(coeff_err[1],2)),ha='left',transform=ax1.transAxes)

    plt.show()


