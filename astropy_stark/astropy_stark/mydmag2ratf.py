# script to convert differences in magnitudes to ratio of fluxes

import numpy as np

## Input file
inputfile=raw_input('enter input file name...')
print ''
band = raw_input('Enter sloan band (u,g,r,i,z)...')
skiprows=0


def dmag2ratf(dmag):
    ratf=10**(-0.4*dmag)