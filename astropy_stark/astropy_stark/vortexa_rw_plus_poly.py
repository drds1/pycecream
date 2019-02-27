import numpy as np
import matplotlib.pylab as plt
import vortexa_polyfit as vpf
import vortexa_rw_simple as vrw




def rwp(t,y,tgrid=[],order=3,custom_freqs = []):

 nt = np.shape(t)[0]
 if (tgrid == []):
  tg = np.array(t)
 else:
  tg = tgrid
  
 si = np.ones(nt)  
  
  
 #fot the polynomial
 poly = vpf.fit(t,y,si,order=order,xgrid=tg,confidence=0.3173,nits=20000,figure_title='')
 #return(yg_med,yg_lo,yg_hi,cov,cisq,cisq_red,bic,aic,rmd)
 ygp = poly[0]
 sgp = (poly[2] - poly[1])/2
 #interpolate the poly model onto data
 yp_itp = np.interp(t,tg,ygp)
 
 
 #subtract the poly model from the data
 ysub = y - yp_itp
 
 #fit the random walk to the polynomial-subtracted data
 rw = vrw.rw(t,ysub,tgi = tg,fbreak=-1,custom_freqs = custom_freqs)
 #return(tgrid,ygridop,f,ckout,skout)
 a = rw[1]
 ygrwlo,ygrwm,ygrwhi = a[:,0],a[:,1],a[:,2]
 sgrw = (ygrwhi - ygrwlo)/2
 
 
 
 #combine the two signals
 yg_tot = ygrwm + ygp
 sg_tot = np.sqrt(sgrw**2 + sgp**2)
 
 
 #fit statistic - autocorrelation function
 yg_tot_itp = np.interp(t,tg,yg_tot)
 rmd = np.corrcoef(yg_tot_itp,y)[0,1]
 
 
 
 #return the output
 return(yg_tot,sg_tot,rmd)












##test model on fake data
#
#from mylcgen import *
#from myrandom import *
#import myresample as mrs
#timeall = []
#sigall  = []
#yall    = []
#
#nlc = 1
#
#tlo = 0.0
#thi = 100
#dtave = 1.0
#dtres = 0.1
#
#for i in range(nlc):
# #generate test light curve and add noise
# datpre      = mylcgen(tlo=tlo,thi=thi,dt=dtres,iseed=132423)
# npre     = np.shape(datpre[:,0])[0]
# datmean  = np.std(datpre[:,1])
# snow =  np.ones(npre)/10*datmean 
# dat = mrs.myresample(dir='',fname=[''],dtave=dtave,sampmin=0.5,sampcode=3,datin=np.array((datpre[:,0],datpre[:,1],snow)).T)
# ndat = np.shape(dat[:,0])[0]
# sig = dat[:,2]
# for i in range(ndat):
#  dat[i,1] = normdis(1,dat[i,1],sig[i])[0]
# sigall.append( sig )
# yall.append( dat[:,1] +50 )
# timeall.append( dat[:,0] )
# 
# 
#y = yall[0]
#t = timeall[0]
#s = sigall[0]
# 
# #dat =mylcgen(tlo=tlo,thi=thi,dt=dtave,iseed=132423)
# #nd = np.shape(dat[:,0])[0]
# #timeall = [dat[:,0]]
# #yall = [dat[:,1]]
# #sigall = [np.ones(nd)/10*datmean ]
# 
#y = y + 0.0*t**2
#print('going in')
#forecast = 50
#dtm = 0.4
#tmod = np.arange(tlo-forecast,thi+forecast,dtm)
#ntmod = np.shape(tmod)[0]
#xg = np.zeros((ntmod,3))
#ygrid,siggrid,rmid = rwp(t,y,tgrid=tmod,order=3)#rw(t,y,si=s,tgi = tmod,fbreak = 0.00005)
#xg[:,1] = ygrid
#xg[:,0] = ygrid - siggrid
#xg[:,2] = ygrid + siggrid
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.errorbar(t,y,s,ls='')
#ax1.fill_between(tmod,xg[:,0],xg[:,2],alpha=0.5)
#
#plt.show()
##from vortexa_mod_randomwalk import *
##a = fitrw([t],[y],[s],plot_tit='sh