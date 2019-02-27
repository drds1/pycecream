r_out=100000

inc=0
D=1
M=1*10**8
M_dot=1
R_star=3
size=10000
from trprof_3 import *


def r_sch(M):
	rs=2*const.G*M/const.c**2
	return(rs)

r_in=3.*r_sch(M*const.Msun)
### PLot the spectra at severla redshifts
#redshifts=[0.1,0.15,0.2,0.3,0.5,1,2]
Mass=[1*10**8,2*10**8,3*10**8,4*10**8,5*10**8]

#rad=[100/const.ps*r_sch(const.Msun*M),1000/const.ps*r_sch(const.Msun*M),10000/const.ps*r_sch(const.Msun*M),100000/const.ps*r_sch(const.Msun*M)]
spectra=np.zeros((size,len(Mass)))
w_l=np.zeros((size,len(Mass)))
count=0
l=[]
for D in (Mass):
	l.append(D)
	temp=flat_disk(r_out/const.ps*r_sch(M*const.Msun),r_in,inc,1,D,M_dot,R_star)
	w_l[:,count]=temp[0]#np.log10(temp[0])
	spectra[:,count]=temp[1]#np.log10(temp[1])
	plot(np.log10(w_l[:,count]),np.log10(spectra[:,count]),linewidth=2,marker='+')
	count=count+1

text=temp[2]  ## plotting information


xlabel('log(Wavelength / m)')
ylabel('log(flux Wm^-2)')
xticks(fontsize=15)
yticks(fontsize=15)
title ('Z = '+str(text[0])+' .  r_in/r_out (pc)= 3/' +str(r_out)+'.  M='+str(text[4])+'.  M dot='+str(text[3])+'.  Inc='+str(text[5]))
#title('z=1 r_in (pc)= 3rs. M= '+str(text[4])+'/Msun.   M dot='+str(text[3])+'.  Inc='+str(text[5]))	

legend(l)
show()

