import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pylab as plt

#function to identify the contribution of constituent signals to a prent signal

def constituents_sklearn(y_main,y_parts):

 reg = linear_model.LinearRegression()
 reg.fit(y_parts,y_main)                                      
 linear_model.LinearRegression(copy_X=True, fit_intercept=True, n_jobs=None,
                 normalize=False)

 coefficients = reg.coef_
 
 
 #split into train and test
 ndat,ndim = np.shape(y_parts)
 frac = 0.85
 ntest = np.int(frac*ndat)
 yt_main = y_main[:ntest]
 yt_parts = y_parts[:ntest,:]
 # Create linear regression object
 reg_test = linear_model.LinearRegression()
 # Train the model using the training sets
 reg_test.fit(yt_parts, yt_main)
 # Make predictions using the testing set
 ytp_main = regr.predict(yt_parts[ntest:,:])
 mse = mean_squared_error(y_main[ntest:],ytp_main)
 # Explained variance score
 yp_main = reg.predict(y_parts)
 r2s = r2_score(y_main,yp_main)

 
 return(coefficients,mse,r2s)
 
 
 
 
 
 
#custom code to fit the glm to multivariate time series components
#y_main(npoints): numpy array of main time series (flows)
#y_parts(npoints,ndimension): 2d numpy array of points for each proposed constituent time series
def constituents_fit(y_main,y_parts,trainfrac = 0.8):
 
 ny = len(y_main)
 ny,npatterns = np.shape(y_parts)
 
 #normalise
 yp_mean = 0#np.mean(y_parts,axis=0)
 yp_std  = 1#np.std(y_parts,axis=0)
 ym_mean = 0#np.mean(y_main)
 ym_std  = 1#np.std(y_main)
 yn_main = (y_main - ym_mean)/ym_std
 yn_parts = (y_parts - yp_mean)/yp_std
  
 
 #compute hes and covariance
 def fit(x,y):
  hnow  = np.tensordot(x.T,x,axes=1)
  cnow  = np.dot(y,x)
  cov   = np.linalg.inv(hnow)
  parms = cov.dot(cnow)
  return(parms,cov)
 parms,cov = fit(yn_parts,yn_main)
 y_mod = np.sum(parms*yn_parts,axis=1)
 
 #compute the r2 score describes how much of the variance the model explains
 var = np.std(y_main)
 var_mod = np.std(y_main-y_mod)
 r2 = 1. - var_mod/var
 
 
 #train the model on a fraction of the input data to compute the model mse
 idtrain = np.int(frac*ny)
 parms_train, cov_train = fit(yn_parts[:idtrain,:],yn_main[:idtrain])
 y_mod_predict = np.sum(parms_train*yn_parts[idtrain:,:],axis=1)
 mse = np.sum((y_mod_predict - yn_main[idtrain:])**2)/(ny-idtrain)
 
 #importance will be the contribution of each amplitude to the model
 #(a measure of the contribution of each component to the explained variance)
 importance = parms**2/np.sum(parms**2)
 
 return(parms,cov,r2,mse,importance)
 
 
 
 
 
def constituents_predict(parms,x):
 predict = np.sum(parms*x,axis=1)
 return(predict) 
 
 
 
 



def constituents_diagplot(importance,labels=None,fig_title='show'):
 ncompare = len(importance)
 color = ['r','g','b','k','purple','cyan','orange','yellow','skyblue'] 
 fig = plt.figure()
 ax1 = fig.add_subplot(111)
 xbar = np.arange(ncompare)+1
 barwidth = 0.8
 
 for i in range(ncompare):
  npar = len(importance[i])
  c = color*npar
  if (labels == None): 
   lab_ann = ['parameter '+np.str(i2+1) for i2 in range(npar)]
   y = np.cumsum(importance[i])
   x = [xbar[i]]*npar
   ax1.bar(x,y[::-1],width=barwidth,color=c[:npar],label=None)
   for i2 in range(npar):
    if i2 == 0:
     yannotate = y[0]/2
    else:
     yannotate = (y[i2] - y[i2-1])/2 + y[i2-1]
    ax1.annotate(lab_ann[i2],(xbar[i]-barwidth/2,yannotate))
 ax1.set_ylabel('explained variance fraction')
 ax1.set_xlabel('model number')
 ax1.set_xticks(xbar)
 if fig_title == 'show':
  plt.show()
 else:
  plt.savefig(fig_title)
 return()
 
 
 

#test on fake data 
x = np.arange(0,100,0.1)
y = 2.3*x + 4.*x**2 + 0.3*x**3#np.sin(2*np.pi*x/30)
y = y
ypart = [x,x**2,x**3]
ypart = np.array(ypart).T
ny,ndim = np.shape(ypart)


#test sklearn glm code
reg = linear_model.LinearRegression()
reg.fit(ypart,y)                                      
a=linear_model.LinearRegression(copy_X=True, fit_intercept=True, n_jobs=None,
                normalize=True)



#test custom code
frac = 0.9
ntest = np.int(frac*ny)
ytrain = y[:ntest]
ytest = y[ntest:]
ypart_train = ypart[:ntest,:]
ypart_test  = ypart[ntest:,:]
parms,cov,r2,mse,importance = constituents_fit(ytrain,ypart_train,trainfrac = 0.8)
ypredict = constituents_predict(parms,ypart_test)





#diagnostic figure
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.scatter(x[:],y[:],color='b',label='Data')
ax1.scatter(x[ntest:],ypredict,color='r',label='prediction')

ax2=fig.add_subplot(212)
xbar = np.arange(ndim)
ax2.bar(xbar,importance*r2,label='explained variance')
ax2.set_xlabel('parameter')
ax2.set_ylabel('explained variance')

plt.show()


constituents_diagplot([importance]*3)






