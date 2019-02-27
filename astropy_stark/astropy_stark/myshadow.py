import numpy as np

#inputs r(n,3) position vectors of all elements 
#a(n) area of elements
#norm(n, 3) normal vector for each surface element, 
#obs(3) unit vector to the oobserver

#returns
#shadow(n) integer vecor reads 1 if surface element can be seen or 0 if in shadow



#check all vectors are unit
def shadow(r,a,norm,obs):
 rmag    = np.sqrt(np.sum(r**2,axis=1))
 ru      = r/rmag[:,None]
 
 normmag = np.sqrt(np.sum(r**2,axis=1))
 normu   = norm/normmag[:,None]
 
 obsmag  = np.sqrt(np.sum(obs**2))
 obsu    = obs/obsmag
 
 nr = np.shape(r[:,0])[0]
 
 shadow = np.ones(nr)
 
 for i_i in range(nr):
  
  d_esubi = obs - r[i_i,:]
  mag_esubi = np.sqrt(np.sum(d_esubi**2))
  
  
  if (a[i_i] <= 0):#exclude lampposts
   shadow[i_i] = 0
   break
  for i_n in range(nr):
   if (i_i == i_n):#dont check point with itself
    continue
   if (a[i_n] <= 0):#exclude lampposts
    continue
   
   d_esubn = obs - r[i_n,:]
   mag_esubn = np.sqrt(np.sum(d_esubn**2))
   
   #if  closer to observer than n, dont bother with rest of loop because i cant be blocked by n 
   if (mag_esubn > mag_esubi):
    continue 
   
   d_nsubi = r[i_n,:] - r[i_i,:]
   mag_nsubi = np.sqrt(np.sum(d_nsubi**2))
   
   
   u_nsubi   = d_nsubi/mag_nsubi
   u_esubn   = d_esubn/mag_esubn
   
   #print 'predot left', u_nsubi*u_esubn
   dotleft = np.sum(u_nsubi*u_esubn)
   sinleft = np.sin(np.arccos(dotleft))
   
   left = mag_nsubi*sinleft#dotleft
   
   dotright = np.sum(u_esubn*normu[i_n])
   right = np.abs(dotright*np.sqrt(a[i_n]))
   
   
   if (left != left):
    print mag_nsubi,dotleft,u_nsubi,u_esubn,d_nsubi,i_i,i_n, 'problem with shadow code myshadow.py'
    print d_nsubi,mag_nsubi,r[i_n,:], r[i_i,:]
    raw_input()
   print i_n,r[i_n]
   print i_i,r[i_i]
   print i_i,left,right
   
   #right is the linear component of the area element projected into the direction of element n toward the observer, left is the distance between the line of site of n to point i, need right > left for point to be in shadow 
   if (right > left):
    shadow[i_i] = 0
    break
   
 return(shadow)
   
  
  
  
  
  
#code to test the shadow code
r = np.array(([0.,1.,0.],[0.,2.,0.]))
a = np.array([1.,1.])
norm = np.array(([0.,1.,0.],[0.,1.,0.]))
obs = np.array([0.,50.,0.])
idxshadow = shadow(r,a,norm,obs)







  