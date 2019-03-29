## module to return the black body plank function of an input temperature and wavelength
#output blackbody intensity (erg/cm2/s/Hz/ster)


#update 6/8/17 include feature to go from erg/cm^2/s/Hz to mJy
#update 6/8/17 include feature to go from erg/cm^2/s/Hz to erg/cm^2/s/A
#using conversions below http://www.stsci.edu/~strolger/docs/UNITS.txt
#[Y Jy]                        = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
#[Y Jy]                       = 1.0E+26 * [X W/m^2/Hz]
#[Y Jy]                        = 1.0E+23 * [X erg/cm^2/s/Hz]
#[Y erg/cm^2/s/Hz]             = 3.33564095e-19 * [X2 A]^2 [X1 erg/cm^2/s/A]


import numpy as np
def bnu(wave,temp,mjy=1,ergcmang = 0):
 
	c1 = 1.43883e8  ### c1 = hc/k (in cgs units) *10^8 as the wavelength is given in angstroms
	c2 = 1.95722e5  ### c2 = (c1 / (2hc))**1/3
	
	
	##convert to si
	
	#c1=1.43916e-2
	#c2=1.9577e6
	BNU = 0.
	
	if (temp > 0.0):
		X = c1/(wave*temp)
		if X > 85:
			bnuln=3.*np.log((c1/wave)/c2)-X
			bnu=np.e**(bnuln)
		if X < 1e-4:
			factor= 2. / ( X * ( X + 2. ) )
			X = X *temp /c2
			bnu = factor * X**3
		if X < 85 and X > 1e-4:
			factor = 1./ (np.e**(X) - 1.)
			X = X *temp /c2
			bnu = factor * X**3
	else:
		bnu = 0.0
	
	
	if (mjy == 1):
		bnu = bnu * 1.e26# * wave*wave * 3.33564095e7
	elif (ergcmang ==1):
		bnu = bnu/3.33564095e-19 / wave/wave
	return(bnu)




#if you have an array of EITHER  wav or temp then use the vectorized version
def bnuvec(wav,temp,mjy=1, ergcmang = 0):
	a = np.vectorize(bnu)
	return( a(wav,temp,mjy=mjy, ergcmang = 0))

