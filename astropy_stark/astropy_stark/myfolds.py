#take an input datasset with class variables and return nfold sets of nsub (<= Ndat - Ntest) subsamples of this data 
# for cross validation. Ntest will be extracted from the orriginal dataset to serve as the test comparison

import numpy as np
def my_folds(data,clas,nfolds,ntest,nsub):
 ndat,ndim = np.shape(data)
 
 #randomly select ntest samples to use as the benchmark accuracy test and remove these from the input
 #data prior to nfolds subselection
 idtest = np.random.choice(np.arange(ndat), size=ntest, replace=False, p=None)
 dattest   = data[idtest,:]
 clastest  = clas[idtest]
 datasub = np.delete(data,idtest,axis=0)
 classub = np.delete(clas,idtest)
 #select nfold new training data sets based on subsample of orriginal training data
 #for cross validation
 ndatnew = ndat - ntest
 datfold = []
 clasfold = []
 for i in range(nfolds):
  idtest = np.random.choice(np.arange(ndatnew), size=nsub, replace=True, p=None)
  dsnow = datasub[idtest,:]
  csub  = classub[idtest]
  datfold.append(dsnow)
  clasfold.append(csub)
 
 return(datfold,clasfold,dattest,clastest)