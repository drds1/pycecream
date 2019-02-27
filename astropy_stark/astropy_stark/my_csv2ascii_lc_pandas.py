import numpy as np
import sys
import pandas as pd
import os


def filt_id(fnow):
 #identify wavelength
 if (fnow == 'B'):
  wav = 4392.0
 elif (fnow == 'V'):
  wav = 5468.0
 elif (fnow == 'HX'):
  wav = 4.4
 elif (fnow == 'SX'):
  wav = 25.3
 elif (fnow == 'U'):
  wav = 3465.0
 elif (fnow == 'UVM2'):
  wav = 2246.0
 elif (fnow == 'UVW1'):
  wav = 2600.0
 elif (fnow == 'UVW2'):
  wav = 1928.0
 
 else:
  wav = 0.0
  print 'unknown wavelength filter:', fnow
  print 'aborting...'
  sys.exit()
 
 return(wav)









path = '/Users/ds207/Documents/standrews/sta/lightcurves/rick_newlc_feb2018'
filename = 'swift4agn.xlsx'#'nearfinal.xlsx'

dirsave = path


df = pd.read_excel(path + '/' + filename)
print df.head()





#identify object colun
o=df[u'Object']
f=df[u'Filter']

ounique = list(set(o))
funique = list(set(f))
wavunique = [filt_id(fn) for fn in funique]
idwavsort = np.argsort(wavunique)
funique = [funique[idw] for idw in idwavsort]

for onow in ounique:
 dsave = dirsave+'/'+onow
 os.system('rm -rf '+dsave)
 os.mkdir(dsave)
 
 
 
 #make creamnames file
 f = open(dsave+'/creamnames.dat','w')
 
 for fnow in funique:
  
  wav = filt_id(fnow)
  
  

  
  
  a = df.loc[(df[u'Object']==onow) & (df[u'Filter']==fnow),['MJD','Flux','Error']]
  asort = a.sort_values(['MJD'],ascending=True)
  asort = asort.round(9)
  fname = onow+'_'+fnow+'.dat'
  asort.to_csv(dsave+'/'+fname,header=False,index=False,sep=' ')
  
  f.write('"'+fname+'" '+np.str(wav)+'\n')
  
 f.close()
  
  
  
#iobject = 0
#ifilter = 1
#itime   = 3
#iflux   = 4
#isig    = 5
#
#
#
#
##df = pd.read_excel(path + '/' + filename,columns=['object','filter','cadence','mjd','flux','error','coad'])
##print df.head()
#
#
#with open(path+'/'+filename) as f:
#    content = f.readlines()
## you may also want to remove whitespace characters like `\n` at the end of each line
#content = [x.strip() for x in content] 
#nc = len(content)
#f.close()
#
#
#for i in range(nc):
# fnow = content[i].split(',')
# 
