# script to convert all files of type PS in directory to png

import numpy as np
import glob
import os


def myps2png(dir):


 pwd = os.getcwd()
 os.chdir(dir)
 fname = glob.glob('*.PS')
 nf = len(fname)


 for i in range(nf):
  file = fname[i]
  idxext  = file.find('.')
  com = 'convert -alpha off '+file+' '+file[:idxext]+'.png'
  os.system(com)
 os.chdir(pwd)