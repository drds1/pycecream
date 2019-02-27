# monte carlo code to fit a brokem power law
import numpy as np

dat=np.loadtxt('testinput.dat')
x=dat[:,0]
y=dat[:,1]
sig=dat[:,2]


# define the model
def model(w,p0,w0,alpha):
    return( p0 / (1.+ (w/w0)**alpha ) )
model = np.vectorize(model)

nparm=3
nits=100
p=np.zeros((1,nparm))
pscale = np.zeros(nparm)
pnew=np.zeros((1,nparm))
pold=np.zeros(1,nparm))

prej=np.zeros(nparm)
ncheck =10
nrej = 4

# best guesses
p[0]=1.0
p[1]=0.01
p[2]=2
pscale[:]=0.1

pold=1.*p
pnew=1.*pold
ymod = model(x,p[0],p[1],p[2])
cisqold = ( (y - ymod) / sig )**2


for iteration in range(nits):

    for iparm in range(nparm):
        
        pold = p[iparm]
        pstep = pscale[iparm]
        
        p[iparm]=pold + myrandom.normdis(1,pold,pstep)
    
    
    
    