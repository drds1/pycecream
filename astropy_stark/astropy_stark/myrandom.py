import numpy as np

pi=np.pi
e=np.e
##This line generates n random numbers between a and b distributed uniformly

def unirand(no,a,b):
	dat=1.*(b-a)*np.random.rand(no)+1.*a
	return(dat)

	
###This definition is udentical to the previous but returns samples from a Gaussian distribution
def unirand_n(no,sigma,mu,a,b):
	dat_n=1./(sigma**2*2*pi)**0.5 * e**(-(unirand(no,a,b)-mu)**2/(2*sigma**2))
###	dat_n=np.random.randn(no)
	return(dat_n)
	

	
def unirand_e(no,tau,a,b):
	#dat_e=1./tau*e**(-np.random.rand(10**4)/tau)
	dat_e=1./tau*e**(-unirand(no,a,b)/tau)
	return(dat_e)
	
	
	

def normdis(n,mu,sigma):
#draw n random numbers from gaussian distribution with mean mu and sd sigma
 op = np.random.randn(n)*sigma + mu
 return(op)	
	
def normdisoldandbroken(n,mu,sigma):	
# for some reason, must set n to one otherwise get odd distribution of random numbers 24/10/2014

#n=1000000
#sigma=1
#mu=0
## draw n random numbers uniformly
	view_out = []
	count=0
	while count < n:
		randval=unirand(n, mu-10*sigma, mu+10*sigma)
		randval.sort()

		randprob=np.zeros((n,2))
		randprob[:,0]=1./(sigma**2*2*pi)**0.5 * e**(-(randval-mu)**2/(2*sigma**2))
		randprob[:,0]=randprob[:,0]/np.max(randprob[:,0])     ##the probablility that the chosen data point will remain
		randprob[:,1]=unirand(n,0,1)

#remove all data points where the random number just drawn between 0 a d 1 is less than the gausian prob for each f the original random numbers
		view=np.delete(randval,[np.where(randprob[:,0] < randprob[:,1])])
	
		view_out=np.append(view_out,view,axis=1)
		count = len(view_out)
		#print 'generating gaussian... ' + str(count) + ' of ' + str(n)
	## delete elements in the view_out array untill we have the correct number of points
	#diff= len(view) - n
	view_out=view_out[:n]
	return(view_out)
#normdis=np.vectorize(normdis)	
