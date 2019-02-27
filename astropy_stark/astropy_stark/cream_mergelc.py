import numpy as np
import glob

#read the creanames file and load wavelengths and files
def read_creamnames(dir,file='creamnames.dat'):

 f = open('../../creamnames.dat')
 filedat = f.readlines()
 filedat = [i.strip() for i in filedat]#remove \n
 wav     = [i.split(' ')[-1] for i in filedat]
 wav = [float(x) for x in wav]

 return(filedat,wav)
 
 
 
 

def cream_merge(dirin,filecontref=0,opfile_tf='',file_lc='./plots/modellc_merge.dat',file_lcsig='./plots/modellc_sig_merge.dat',opfile_tf='./plots/modeltf_merge.dat',file_tfsig='./plots/modeltf_sig_merge.dat',opfile_filewav='./outputpars3_merge.dat',file_osst='./outputpars2_merge.dat',file_fileth='./outputpars_th_merge.dat',file_ve='./outputpars_varexpand_merge.dat',file_creamnames='../creamneames_merge.dat'):
 
#cream_merge(dir,wav='',
#file='',
#opfile_tf=''
#,file_lc='./plots/modellc_merge.dat',
#file_lcsig='./plots/modellc_sig_merge.dat',
#opfile_tf='./plots/modeltf_merge.dat',
#file_tfsig='./plots/modeltf_sig_merge.dat',
#opfile_filewav='./outputpars3_merge.dat',
#file_osst='./outputpars2_merge.dat',
#file_fileth='./outputpars_th_merge.dat',
#file_ve='./outputpars_varexpand_merge.dat',
#file_creamnames='../creamneames_merge.dat'
  
 if 'output_20' in dirin:
  dir = dirin  
 else:
  print 'output file not in dir string...',dirin
  dir = glob.glob(dirin+'./output_20*')[-1]
  print 'looking here...',dir
  
   


 #input files are just outputfiles without the _merge.dat
 #load input files
 filelc,wav=read_creamnames(dir+'/../creamnames.dat')
 nlc = len(wav)
 wavu = np.unique(wav)
 nwavu = np.shape(wavu)[0]
 idwavu = [np.where(wav == wavu[iu])[0] for iu in range(nwavu)]
 idth = np.where(wavu == -1)[0]
 wavu = wavu[np.where(wavu != -1)[0]]
 nth  = np.shape(idth)[0]

 
 nop = nth + wavu
 
 datlc = np.loadtxt(file_lc)
 datlc_sig = np.loadtxt(file_lcsig)
 dattf = np.loadtxt(file_tf)
 dattf_sig = np.loadtxt(file_tfsig)
 
 ntgrid = np.shape(datlc[:,0])[0]
 ntaugrid = np.shape(dattf)[0]
 datlc_op = np.zeros((ntgrid,nop+1))
 dattf_op = np.zeros((ntaugrid,nop+1))
 
 tgrid = datlc[:,0]
 taugrid = dattf[:,0]
 
 
 #identify which continuum ight curve to use as the reference when merging
 idcont = []
 if (isinstance(filecontref,list)):
  for iu in range(nwavu):
   idcont.append([i for filelc[i] in range(nlc) if filecontref[iu] in filelc[i]][0])
 else:
  idcont.append(wavu[) 
 
 
 #model light curves: compute the median and standard deviation
 idxinc = [np.where((datlc[:,1+ilc] != 0))[0] for ilc in range(nlc)]
 median = [np.median(datlc[idxinc[ilc],1+ilc]) for ilc in range(nlc)]
 std    = [np.std(datlc[idxinc[ilc],1+ilc]) for ilc in range(nlc)]

 
 #populate the model arrays for the continuum wavelength
 for iw in range(nwavu):
  wavunow = wavu[iw]
  inow = np.where(wav == wavunow)[0]
  
  
  
  
  
  
  
  
  datlc_op[:,0] = tgrid
  datlc_op[:,iw] = inow[
  
  #light curves
  
 
 #if output_20xxx in dir string then look in there else look for subdir in dir string
 
