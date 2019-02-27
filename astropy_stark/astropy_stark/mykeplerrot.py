# script to calculate the keplarian rotation velocity of a flat disc at an input radius in ld

import numpy as np

G  = 6.673e-11
M0 = 2e30
ld = 2.59020684e13

print 'Welcome to mykeplerrot.py, The code that calculates rotation velocities and wavelength offsets, so you can check it later and find a bug in the code.'
print ''
urld=np.float( raw_input('Please enter disc radius (Ld):-') )
print ''
uembh=np.float( raw_input('PLease enter black hole mass (M0):-') )

print ''


# Temperature info ##################################################

rs=3000*uembh
r0rs=3.0
r0=r0rs*rs

urrs=urld*ld/rs

apower=0.75
bpower=0.75
pi=np.pi
emdot=1.0
year=365*24*3600
udmdt=M0/year
s=5.67e-8
c=2.9979e8
eta=0.1
fourpi=4.*pi
zx=3.0 * rs

T0v4=3/(8*pi)*10**(np.log10(G*uembh)+np.log10(emdot)-3*np.log10(r0)+np.log10(M0)+np.log10(udmdt)-np.log10(s))
T0x4=1/fourpi*10**(np.log10(zx)+np.log10(eta*c*c)+np.log10(emdot)+np.log10(udmdt)-3*np.log10(r0)-np.log10(s))

T=np.sqrt(np.sqrt(T0v4*(r0rs/urrs)**(bpower*4)+T0x4*(r0rs/urrs)**(apower*4)))
wav = 2.89711e-3/T *1e10 # peak wavelength at the input radius (in Angstroms)

# end of temperature calculation ##################################################



print ''
print ''

vkms = np.sqrt( G*uembh*M0/ (ld*urld) ) /1000
dwav = wav*vkms/(c/1000)

print 'Rotation velocity (kms-1):- ',vkms
print 'Emitted wavelength (Ang):- ',wav
print 'dwav:- ',dwav




