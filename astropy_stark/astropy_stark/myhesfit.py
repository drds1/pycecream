#!! program to evaluate the uncertainties and parameter values of a model by minimising chi squared
#!! tested 3rd august with ax**2+bx +c and works (code below) 
#subroutine hesfit(x,y,sig,p,N,NP,pout,sigpout)

#inputs x(N),y(N),sig(N) x y and uncertainty values
# p(N,NP) the pattern to be fitted i.e fitting y=ax^2 + bsin(x) p(1:N,1) would be x**2 and p(1:N,2) would be an array of sin(x)

# outputs: pout(NP), sigpout(NP parameters and uncertainties

import numpy as np






#!!! test this using y=ax + b
#real y(10),x(10),sig(10),p(10,3),pout(3),sigpout(3),ymod(10)
#
#n=10
#n_p=3
#a=1.2
#b=2.8
#c=5.0

#x=np.arange(n)+1
#y=a*x+b*x**2 + c
#sig=np.ones(n)*0.2


#!!! set up patterns
#p=np.zeros((n,n_p))
#hes=np.zeros((n_p,n_p))
#p[:,0]=x
#p[:,1]=x**2
#p[:,2]=1



#below is what will be inside subroutine


#!! set up the Hessian matrix

def hesfit(x,y,sig,p):
 
 
    sig2=sig*sig
    n_p = p[0,:].shape[0]
 
#calculate the hessian matrix    
    for ip2 in range(n_p):    
        for ip1 in range(n_p):
            hes[ip1,ip2]=np.sum(p[:,ip1]*p[:,ip2]/sig2)

# calculate the covariance matrix    
    cov = np.linalg.inv(hes)
    
    
    
#!!! set up the c matrix in:  a * hes = c (where a is a(NP), hes is the hessian matrix and c = sum (y_i p_i / sig^2_i) see ada lecture 9)
    c=np.zeros(n_p)
    pout=np.zeros(n_p)
    sigpout=np.zeros(n_p)
    
    for ip in range(n_p):
        c[ip]=np.sum(y*p[:,ip]/sig2)
    
    #!stop
    #!! find the solution 
    for ip0 in range(n_p):
        pout[ip0]=np.sum(cov[:,ip0]*c)
        sigpout[ip0]=np.sqrt(cov[ip0,ip0])
    
    #!!
    return(pout,sigpout)

