import pyfits
import numpy




def loadsdssspec(infile, redshift=False):
   """
   Load a 1D SDSS spectrum and extract wavelength solution.
   This does not work for DR9 spectra.
   @param infile: input file name
   @type infile: string
   @param redshift: redshift spectrum?
   @type redshift: bool
   """
   dat = pyfits.open(infile)
   if 'COEFF0' in dat[0].header is False or 'COEFF1' in dat[0].header is False:
       raise Exception('No proper header.')
   else:
       co0 = dat[0].header['COEFF0']
       co1 = dat[0].header['COEFF1']
   if dat[0].data is None:
       raise Exception('This is most likely a DR9 spectrum!!!!')
   dim = numpy.shape(dat[0].data)
   if numpy.shape(dim) == (1,):
       n = len(dat[0].data)
       flux = dat[0].data
       err = numpy.zeros_like(flux)
   else:
       n = len(dat[0].data[0])
       flux = dat[0].data[0]
       err = dat[0].data[2]
   wl = 10**(co0 + co1*numpy.arange(n))
   if redshift:
       wl /= 1 + dat[0].header['Z']
   dat.close()
   return {'wl': wl, 'flux': flux, 'error': err, 'z': dat[0].header['Z']}