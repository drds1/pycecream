import vortexa_makefake as vmf
import numpy as np

#dat ntimes x 2 matrix first column is the driver second is the response
def rli(dat):

 n = len(dat[:,0])
 npsi = n
 dmk = np.mean(dat[:,0])
 drive = dat[:,0]
 #smoothing matrix (see thesis page 123)
 s = np.zeros((n,npsi))
 for ik in range(n):
  if (ik == 0):
   s[ik,0] = 2
  elif (ik == 1):
   s[ik,1] = 2
  elif (ik == n-1):
   s[ik,-1] = 2
  elif (ik == n-2):
   s[ik,-2] = 2
  else:
   s[ik,ik-2:ik+3] = np.array([0.5,-2,3,-2,0.5])


 #driving matrix of look back times
 dmat = np.zeros((n,npsi))
 for ik in range(n):
  ik2 = np.arange(npsi)
  ilb = np.arange(ik,ik-npsi,-1)
  ilb[ilb < 0] = 0
  dmat[ik,ik2]=drive[ilb]
  #for ik2 in range(npsi):
  # idx = max(ik - ik2,0)
  # dmat[ik,ik2] = drive[idx]



#now construct hessian matrix from dmat and sigma 2 see also page 123 of thesis

 return()


#for ik in range(npsi):
# for ij in range(1,npsi):
#
##sum = 0.0
#  for it3 in range(1,nf):
#   sum = sum + (xitp(it3,ij) - dmj)*(xitp(it3,ik) - dmk)/sig2(it3)
#!write(*,*) it3,xitp(it3,ij),xitp(it3,ik),sig2(it3)
#enddo
#hesmat(ij,ik) = sum
#
#
#sum = 0.d0
#!now add constraint from smoothness function
#if ((ik .gt. 2) .and. (ik .lt. npsi-1)) then
# if ((ij .eq. ik - 2) .or. (ij .eq. ik + 2)) then
# sum = 0.25   !0.5
# else if ((ij .eq. ik -1) .or. (ij .eq. ik + 1)) then
# sum = -1.        !*2
# else if (ij .eq. ik) then
# sum = 1.5    !3
# endif
# 
# 
# !if (ij .eq. ik - 1) then 
# !sum = sum - alpha*0.5
# !else if (ij .eq. ik) then 
# !sum = sum + alpha*1.
# !else if (ij .eq. ik + 1) then
# !sum = sum - alpha*0.5
# !else
# !sum = sum - alpha*0.
# !endif
#
#else if ((ik .le. 2) .or. (ik .ge. npsi - 1)) then
# if (ij .eq. ik) then
# sum = 1. !+ alpha*1.!2
# else
# sum = 0.
# endif
#endif
#
#!if (ij .eq. 1) then  !! edge effects. This minimises psi^2 for the edges to pull tf down to 0 at ends
#!sum = sum - 0.0!1.e6*alpha*1.
#!else
#!sum = sum + alpha*0.
#!endif
#!else if (ik .eq. npsi) then
#!if (ij .eq. npsi) then
#!sum = sum - 0.0!1.e6*alpha*1.
#!else
#!sum = sum + alpha*0.
#!endif
#
#!endif
#!else if (ik .eq. 1) then !what to do about edge effects (ignore for now)
#!if (ij .eq. 1) then
#!sum = sum + alpha*1.
#!else if (ij .eq. 2) then
#!sum = sum + alpha*
#!endif
#
#smoothmat(ij,ik) = sum
#
#
#
#!write(*,*) ij,ik,sum
#!read(*,*)
#enddo














tlo = 0.0
thi = 100.0
dt = 0.01

a = vmf.makefake(tlo,thi,dt,shift = [2,3],noise=0.1)


din = np.array([a[0].values[:,1],a[1].values[:,1]]).T
rli(din)



