#python code to take x,y,sigy return bet fit with uncertainty envelopes
import numpy as np

def myoptav(y,sigy):
    
    sigy2=sigy*sigy
    
    top=(y/sigy2).sum()
    bot=(1/sigy2).sum()
    optav=top/bot
    return(optav)
    



def mylinfitenv(x,y,sigy):

    xbar=myoptav(x,sigy)

# intercept b (optimal average)

    top=0.0
    bot=0.0
    si2=sigy*sigy
    temp=1/si2
    bot=temp.sum()
    top=(y*temp).sum()
    b=top/bot
    sigb=np.sqrt(1/bot)
    
    #print b,sigb,'b and sig b mylinfit.py'
    #bo=raw_input()
# slope (optimal scaling of patter)

    xsubbar=x-xbar
    temp=xsubbar/si2

    bot=(temp*xsubbar).sum()
    top=((y-b)*temp).sum()


    a=top/bot
    siga=np.sqrt(1/bot)

    #print a,siga,top, bot, xbar,'a and sig, top, bot, a mylinfit.py'

    
    #bo=raw_input()
#!!! this gets you the orthogonal slope(a, siga) but the intercept is just the optimal average of points
#! in the y direction. To get a reliable value for the intercept, need to find out where the curve
#! hits the axis intercept = ybar - (xbar * a)  where ybar = old intetcept b 
#! the above is in Keiths ada notes lecture 8, the below is written down on some long lost scrap of paper.
   
    ybar=b
    
    b = b - a*xbar
    
    return(b,sigb,a,siga,xbar,ybar)
    
    
    
    
# end of function





#
#x=np.arange(10)+1
#y=1.*x
#sigy=np.ones(10)*0.1
#y[3]=4.2
#y[6]=6.8
#
#dat=mylinfitenv(x,y,sigy)
#
#print dat



