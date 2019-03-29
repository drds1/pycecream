import numpy as np
import os

#code to call fortran to calculate the response function


dirfort = '/Users/ds207/Documents/standrews/sta/fort/fortcode'
pwd = os.getcwd()

iexist = os.path.exists('tfb_pycall.exe')


def mytfb_pycall(tau,wav,uembh,uemdot,deginc,slopevisc=-0.75,slopeirad=-0.75,urin=3.0,r0=1.0,eta=0.1,uhx=3.0,forcetau0 = 1,tfx=0,T0v=1.e4,T0x=1.e4,sv=-0.75,sx=-0.75,newversion = 1):
 

 taulo = tau[0]
 tauhi = tau[-1]
 dtau = np.mean(tau[1:] - tau[:-1])

 if (tfx == 0):

  if (iexist == 0 or newversion == 1):
   print('compiling...')
   dirpycallf90 = '/Users/ds207/Documents/standrews/sta/fort/fortcode/delaydist/tfb_pycall.f90'
   os.system('cp '+dirpycallf90+' '+pwd)
   os.system('gfortran tfb_pycall.f90 '+dirfort+'/delaydist/tfb_redshifttest.f90 '+dirfort+'/delaydist/tfb_x.f90 '+dirfort+'/cgs.f90'+' '+dirfort+'/itp2.for '+dirfort+'/ran_new2.f90 -o tfb_pycall.exe')
   newversion = 0
   print('done compiling')

  f = open('tfb_pycall.par','w')
  f.write(np.str(wav)+'\n')
  f.write(np.str(taulo)+' '+np.str(tauhi)+' '+np.str(dtau)+'\n')
  f.write(np.str(urin) + '\n')
  f.write(np.str(uembh) + '\n')
  #f.write(np.str(uemdot) + '\n')
  f.write(np.str(deginc) + '\n')
  f.write(np.str(uhx) + '\n')
  f.write(np.str(eta) + '\n')
  f.write(np.str(-1*slopevisc)+'\n')
  f.write(np.str(-1*slopeirad)+'\n')
  f.close() 
  os.system('./tfb_pycall.exe')
  dat = np.loadtxt('tfb_pycall.op')
 
 else:

  
  if (iexist == 0 or newversion == 1):
   print('tfx mode is on!!!')
   print('compiling...')
   dirpycallf90 = '/Users/ds207/Documents/standrews/sta/fort/fortcode/delaydist/tfx_pycall.f90'
   os.system('cp '+dirpycallf90+' '+pwd)
   os.system('gfortran tfx_pycall.f90 '+dirfort+'/delaydist/tfb_x.f90 '+dirfort+'/myrs2ld.f90 '+dirfort+'/ran_new2.f90 '+dirfort+'/mytr.f90 -o tfx_pycall.exe')
   print('done compiling')
   
  f = open('tfx_pycall.par','w')
  f.write(np.str(wav)+'\n')
  f.write(np.str(taulo)+' '+np.str(tauhi)+' '+np.str(dtau)+'\n')
  f.write(np.str(urin) + '\n')
  f.write(np.str(uembh) + '\n')
  f.write(np.str(deginc) + '\n')
  f.write(np.str(T0v) + '\n')
  f.write(np.str(T0x) + '\n')
  f.write(np.str(-1*sv) + '\n')
  f.write(np.str(-1*sx) + '\n')
  f.close()
  os.system('./tfx_pycall.exe')
  dat = np.loadtxt('tfx_pycall.op')
  
 if (forcetau0 == 1):
  dat[0,1] = 0
 
 return(dat)
 
 