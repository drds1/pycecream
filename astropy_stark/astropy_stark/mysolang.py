## a program to calcualte solid angles

import numpy as np
import const


def sol_ang(r,dr,i,D):
	sa=2*np.pi*r*dr*cos(i)/D**2
	return(sa)
	
