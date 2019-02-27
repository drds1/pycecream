import numpy as np
import os

#str_rep_in is a list [1:ncol][1:nrow] if you want to replace a column by strings, set str_rep_in[icol]='' if dnot want to alter
#code to make a latex table
#aw shift by xxcm
def mytextable(datin=[],filedatin='',opfile='mytable',caption='',colheadin='',ndpin=[],datsd=[],colunits='',str_rep_in=[],sideways=0,aw='0cm',tabref='',verbon = 1):

 
 

 if ((datin ==[]) & (filedatin == '')):
  raise Exception('mytextable.py: Both filedat ad dat arguments cannot be empty!')
 elif (datin == []):
  dat = np.loadtxt(fiiledatin)
 elif (filedatin==''):
  dat = 1.*datin
 else:
  raise Exception('mytextable.py: Something else is wrong!')
 
 
 
 
 ncol = np.shape(dat[0,:])[0]
 nrow = np.shape(dat[:,0])[0]
 
 
 
 if (str_rep_in == []):
  str_rep = ['']*ncol
 else:
  str_rep = str_rep_in
 
 if (datsd == []):
  datsdtab = np.ones((nrow,ncol))*-1.
 else:
  datsdtab = datsd

 
 if (ndpin ==[]):
  ndp = [2]*ncol
 else:
  ndp = list(ndpin)
 
 print ndp
 
 target = open(opfile+'.tex', 'w')
 
 cs = ''
 cs = [cs + 'c ' for i in np.arange(ncol)] 
 cs = ''.join(cs)+'}'
 cs = '{'+cs
 
 ch = ''
 if (colheadin == ''):
  ch = [ch + np.str(i) +' & ' for i in np.arange(ncol)]
 else:
  ch = [ch + i +' & ' for i in colheadin]
 ch = ''.join(ch)[:-2] +' \\\\ \n'

 cu = ''
 if (colunits == ''):
  cu = ''
 else:
  cu = [cu + i +' & ' for i in colunits]
 cu = ''.join(cu)[:-2] +' \\\\ \n'



 
 #write document ifo
 target.write('\\documentclass[10pt]{article}  \n')
 target.write('\\newcommand{\\HRule}{\\rule{\\linewidth}{0.5mm}} \n')
 target.write('\\parindent 0pt \n')
 target.write('\\parskip 10pt \n')
 target.write('\\usepackage{anysize} \n')
 target.write('\\usepackage{graphicx} \n')
 target.write('\\usepackage{rotating} \n')
 target.write('\\usepackage{epsfig} \n')
 target.write('\\usepackage{float} \n')
 target.write('\\usepackage{natbib} \n')
 target.write('\\usepackage{setspace} \n')
 target.write('\\usepackage{changepage} \n')
 target.write('\\marginsize{3.5cm}{3.5cm}{1cm}{1cm} \n')
 target.write('\\onehalfspacing \n')
 target.write('\\usepackage{caption} \n')
 target.write('\\usepackage{subcaption} \n')
 target.write('\\usepackage{amsmath} \n')
 target.write('\\begin{document} \n')
 
 
 #write preamble including table caption
 if (sideways == 1):
  target.write('\\begin{sidewaystable} \n')
 else:
  target.write('\\begin{table} \n')
  
 target.write('\\begin{adjustwidth}{'+aw+'}{} \n')
 target.write('\\caption{'+caption+'} \n')
 target.write('\\centering \n')
 target.write('\\begin{tabular}'+cs+' \n')
 target.write('\\hline\\hline \n')
 target.write(ch)
 if (colunits != ''):
  target.write(cu)
 target.write('\\hline \n')
 
 #now input table data
 #row by row
 
 for i in range(nrow):
  rownow = ''
  for i2 in range(ncol):
   
   
   if (str_rep[i2] == ''):
    
    if (datsdtab[i,i2] == -1):
     if (ndp[i2] != 0):
      rownow = rownow +  np.str(np.round(dat[i,i2],ndp[i2])) +' & '
     else:
      rownow = rownow +  np.str(np.int(dat[i,i2])) +' & '
    else:
     if (ndp[i2] != 0):
      rownow = rownow +  np.str(np.round(dat[i,i2],ndp[i2])) +' $\\pm$ '+ np.str(np.round(datsdtab[i,i2],ndp[i2])) + ' & '
     else:
      rownow = rownow +  np.str(np.int(dat[i,i2])) +' $\\pm$ '+ np.str(np.int(datsdtab[i,i2])) + ' & '
    
   else:
    if (verbon == 1):
     rownow = rownow + '\\verb|'+str_rep[i2][i] + '| & ' 
    else:
     rownow = rownow + str_rep[i2][i] + ' & ' 
  
  #rownow = [rownow + np.str(np.round(x,ndp[0])) +' & ' for x in dat[0,:]]
  #rownow = ''.join(rownow)[:-2] + '\\\\ \n'
  rownow = rownow[:-2] + '\\\\ \n'
  target.write(rownow)
 
 target.write('\\hline \n')
 target.write('\\end{tabular} \n')
 if (tabref != ''):
  target.write('\\label{'+tabref+'} \n')
 target.write('\\end{adjustwidth} \n')
 if (sideways==1):
  target.write('\\end{sidewaystable} \n')
 else:
  target.write('\\end{table} \n')
 
 target.write('\\end{document} \n')
 
 target.close()
 
 
 os.system('pdflatex '+opfile+'.tex '+opfile+'.pdf')
 #os.system('rm '+opfile+'.tex')
 os.system('rm '+opfile+'.aux')
 os.system('rm '+opfile+'.blg')
 os.system('rm '+opfile+'.log')
 os.system('rm '+opfile+'.synctex.gz')

 #print 'pdf2latex not available. I have just mde the tex file you will have to compile later'



#dat = np.ones((10,2))
#datsd = np.ones((10,2))
#datsd[:,-1] = -1
#datsd[3,0] = -1
#colheadin = ['col 1','col 2']
#units = ['','days']
#a = mytextable(datin = dat, datsd = datsd,colheadin = colheadin,colunits = units)
#
#