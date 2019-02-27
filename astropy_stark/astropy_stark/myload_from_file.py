import numpy as np
from pylab import *



## a function to load an n dimensional array and output a table of data



##open the text file and store the numbers in an array

##function to load n d array from file with name 'name', and number of rows n_row

def file_load(n,name,no_row):
	
	file=open(name, 'r')
	l=[]
	a=file.readlines()

##the file contents are now stored line by line in a list . Need to extract all the numbers from the values from the list and store them in an array


##len(a) is the number of columns, no_row is no of rows
## remember, rows and columns are opposite for python and fortran

##extract all the spaces values from file
	data=np.zeros((len(a), no_row))
	for i in range(len(a)):
		data[i,:]=a[i].split()

	return(data)
	
