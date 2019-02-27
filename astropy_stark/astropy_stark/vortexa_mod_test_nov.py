import numpy as np
import astropy.convolution as apc
import vortexa_mod_randomwalk as mf18
import mylcgen as mlc
import myresample as mrs
import itertools
import matplotlib.gridspec as gridspec
import matplotlib.pylab as plt
import pandas as pd
import vortexa_convert_df2lc as c2l
import vortexa_polyfit as vpf
import matplotlib.gridspec as gridspec
import vortexa_rw_plus_poly as vrwp
import vortexa_rw_simple as vrs
import scipy.signal as ss


#make syntax correct for data frame input
def mass_mod_test(meta_df_in,dateref=-1,forecast_period = 90,figure_title='',model='poly',
convolve=None,transform_model='same',labels=[],dtgrid=10.0,normalise=1,polyorder=1,
accuracy_check=10,custom_freqs = [],figure_label = ''):


 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #!!!!!!!!!!!!!!!!!! Data preparation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#4/12/2018 - incorporate accuracy check
#if accuracy check > 0then first fit model to ndat - 10 points and see how accurate the model is
#two metrics 
#1) how many points within uncertainty envelopes as percentage
#2) what is the model variance
 accuracy_envelope = []
 accuracy_variance = []
 accuracy_mape     = []


 nframes = len(meta_df_in)
 color   = ['k','r','b','cyan','purple','green']*nframes

 #identify earliest date from all the time axis and use this as the first time by default
 datemax_save = []
 meta_df = []
 if (dateref == -1):
  dates = pd.concat([meta_df_in[i].iloc[:,0] for i in range(nframes)])
  datemax_save.append([meta_df_in[i].iloc[-1,0] for i in range(nframes)])
  datemin = min(dates)
  datemax = max(dates) + pd.Timedelta(days=forecast_period)
  date_trunc = np.min(np.array(datemax_save))
  print ('Earliest dat in sample...',datemin)
 else:
  datemin = dateref[0]
  datemax = dateref[1] + pd.Timedelta(days=forecast_period)
  date_trunc = dateref[1]


 for i in range(nframes):
  idinclude = np.where((meta_df_in[i].values[:,0] > datemin) & (meta_df_in[i].values[:,0] < datemax))[0]
  df = meta_df_in[i].iloc[idinclude,:]
  #if normalise 1 then normalise
  if (normalise == 1):
   xn = df.values[:,1]
   xn0 = xn[np.nonzero(xn)]
   xnmean = np.mean(xn0)
   xnstd  = np.std(xn0)
   if (xnstd == 0):
    xnstd = 1
   else:
    xnstd = 1
    
   xnew = (xn - xnmean)/xnstd 
   df.values[:,1] = xnew  
   meta_df.append( df )
  
  

 #convert data frames to light curve format
 meta_lc = []
 meta_mod = []
 transform_mod = []
 fitstat_self = []
 fitstat_cross = []
 for i in range(nframes):
  x = c2l.convert_df_2_lightcurve(meta_df[i],reftime=datemin,noise = 0.01,normalise = 0)
  meta_lc.append(x)
 
  
  #model paramers
  if (i == 0):
   tlo = np.min(x[:,0])
   thi = np.max(x[:,0])
   dt  = x[1,0] - x[0,0]
   tgrid = np.arange(tlo,thi+dtgrid+forecast_period,dtgrid)
   flo = 1./(thi-tlo)
   fhi = 0.1/dtgrid
   ntgrid = np.shape(tgrid)[0]
   dategrid = [datemin + pd.Timedelta(days=i*dt) for i in range(ntgrid)]
   print('lower upper time grid limits',tlo,thi)
  
  
  
  
  
 #label each time series
 if (labels == []): 
  lab = ['time series '+np.str(i) for i in range(nframes)]
 else:
  lab = labels

  
  
  












  
  
  
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #!!!!!!!!!!!!!!!!!! fit model to each time series !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 for idx_check in range(2):
  if ((idx_check == 1) and (accuracy_check == 0)):
   accuracy_envelope.append( np.nan )
   accuracy_variance.append( np.nan )
   continue
   

  for i in range(nframes):
   ntimes = np.shape(meta_lc[i][:,0])[0]
   if (idx_check == 1):
    idlast = ntimes - accuracy_check
   else:
    idlast = ntimes
   
   if (model == 'rw'):
    #x = mf18.fitrw([meta_lc[i][:idlast,0]],[meta_lc[i][:idlast,1]],[np.ones(ntimes)*np.std(meta_lc[i][:idlast,1])*0.1],floin=flo,
    #fhiin=fhi,plot_tit='show',dtresin=tgrid,nits = 10000,extra_f=[])
    #(tgrid,ygridop,f,ckout,skout)
    x = vrs.rw(meta_lc[i][:idlast,0],meta_lc[i][:idlast,1],si=0,tgi = tgrid,fbreak=-1)
    xmod = x[1][:,1]
    sigmod = (x[1][:,2] - x[1][:,0])/2
    fitstat = x[5]
   elif (model == 'poly'):
    
    x = vpf.fit(meta_lc[i][:idlast,0],meta_lc[i][:idlast,1],np.ones(ntimes)[:idlast],order=polyorder,xgrid=tgrid
    ,confidence=0.3173,nits=10000) 
    xmod = x[0]
    sigmod = (x[2] - x[1])/2
    fitstat = x[8]#rcoef
   
   elif (model == 'poly-rw'):
    xmod,sigmod,fitstat=vrwp.rwp(meta_lc[i][:idlast,0],meta_lc[i][:idlast,1],tgrid=tgrid,order=polyorder)
   
   elif (model == 'poly-custom'):
    if (custom_freqs == []):
     raise Exception('Must supply a numpy array of fit frequencies when chosing "poly_custom" model')
    
    xmod,sigmod,fitstat=vrwp.rwp(meta_lc[i][:idlast,0],meta_lc[i][:idlast,1],tgrid=tgrid,
    custom_freqs=custom_freqs,order=polyorder)
   #here we asses the model accuracy using cross validation testing on the last 
   #'accuracy_check' number of points
   

   
   
   if (idx_check == 0):
    meta_mod.append(np.array([dategrid,xmod,sigmod]).T) #meta_mod.append(np.array([tgrid,xmod,sigmod]).T) 
    fitstat_self.append(fitstat)
    

    
   else:
    tdat        = meta_lc[i][idlast:,0]
    xdat        = meta_lc[i][idlast:,1]
    ndat        = np.shape(xdat)[0]
    xmod_itp    = np.interp(tdat,tgrid,xmod)
    smod_itp    = np.interp(tdat,tgrid,sigmod)
    xmod_itp_hi = xmod_itp + smod_itp
    xmod_itp_lo = xmod_itp - smod_itp
    var = np.sum( (xmod_itp - xdat)**2 )/ndat
    rms = np.sqrt(var)
    mape = rms/np.mean(xdat)
    idgood = np.where( (xdat > xmod_itp_lo) & (xdat < xmod_itp_hi) )[0]
    ngood  = np.shape(idgood)[0]
    goodfrac = np.float(ngood)/accuracy_check
    accuracy_envelope.append( goodfrac )
    accuracy_variance.append( var )
    accuracy_mape.append( mape ) 
 
 
 
 
 
 
 
 
 
 
 
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 #!!!!!!!!!perform operation on model to transform to flow time #!!!!!!!!#!!!!!!!!
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 for i in range(nframes):
  if (i == 0 or transform_model == 'same'):
   response = meta_mod[i]
   response0 = meta_mod[0]
   t0lo = response0[0,0]
   t0hi = response0[-1,0]
   idinc = np.where((response[:,0] < t0hi) & (response[:,0] > t0lo))[0]
   sd0 = np.std(response0[idinc,1])
   sd  = np.std(response[idinc,1])
   mean0 = np.mean(response0[idinc,1])
   mean   = np.mean(response[idinc,1])
   rnew = (response[:,1] - mean)*sd0/sd + mean0
   sdnew = response[:,2]*sd/sd0
   transform_mod.append(np.array([dategrid,rnew,sdnew]).T)
  
   #plt.clf()
   #plt.plot(dategrid,transform_mod[0][:,1],color='k')
   #plt.plot(dategrid,rnew)
   #print('mean,sd,driver',mean0,sd0)
   #print('mean,sd,response',mean,sd)
   #print('mean,sd,response new',np.mean(rnew),np.std(rnew))
   #print('i',lab[i])
   #plt.show()
   #plt.clf()
   
   
  
  elif (transform_model == 'convolve'):
   lagcent,lagwide = convolve[i]
   lg = np.arange(ntgrid)*dt
   shape_k = np.int(2*(ntgrid-1) + 1)
   kernel = np.zeros(shape_k)
   laggrid = np.zeros(shape_k)
   k0 = np.int(np.floor(shape_k/2))
   laggrid = np.zeros(shape_k)
   laggrid[k0:]=lg
   laggrid[:k0]=lg[-1:0:-1]
   kernel= np.exp(-0.5*((laggrid-lagcent)/lagwide)**2)/(2*np.pi*lagwide)**2
   modnow = meta_mod[i][:,1]
   modmean = np.mean(modnow)
   modstd = np.std(modnow)
   response = apc.convolve(modnow,kernel)
   response_std = np.std(response)
   response_mean = np.mean(response)
   
   response = (response - response_mean)/response_std*modstd + modmean
   transform_mod.append(np.array([dategrid,response,sigmod]).T)#transform_mod.append(np.array([tgrid,response,sigmod]).T)
   #print('asdasdfafafa')
   #input()
   #for it in range(ntgrid):
   # print(dategrid[it],response[it],transform_mod[0][it,1])
   #print('transform check',i)
   #input()
 
  #asses the correlation between transformed model and driver
  mm = np.array(meta_mod[0][:,1],dtype='float')
  rm = np.array(transform_mod[i][:,1],dtype='float')
  rcoef = np.corrcoef(mm,rm)[0,1]
  fitstat_cross.append(rcoef)
 
 
 
 
 
 
 
 
 
 
 
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 ##!!!!!!!!#!!!!!!!!#!!!!!!!! diagnostic plot #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 if (figure_title != ''):
  ndown = nframes + 2
  nalong = 4
  plt.clf()
  fig = plt.figure()
  gs1 = gridspec.GridSpec(ndown,nalong)
  gs1.update(left=0.2,right=0.85,wspace=0.1,hspace=0.0,bottom=0.08,top=0.90)
   

   
  uncertanties_old = []
  uncertainties_new = []
  improvement_potential = np.zeros((ntgrid,nframes))
  for i in range(nframes):
   ax1 = plt.subplot(gs1[i,:])
   #ax1.scatter(meta_df[i].values[:,0],meta_lc[i][:,1],s=4,color=color[i])
   ax1.bar(meta_df[i].values[:,0],meta_lc[i][:,1],color=color[i],width = 25.0,fill=False)
   tmod    = tgrid#meta_mod[i][:,0]
   xmod    = np.array(meta_mod[i][:,1],dtype='float')
   xmodsig = np.array(meta_mod[i][:,2],dtype='float')
   xmodlo  = xmod - xmodsig
   xmodhi  = xmod + xmodsig
   
   if (i == 0):
    for i2 in range(1,nframes):
     xmod_trans = np.array(transform_mod[i2][:,1],dtype='float')
     ax1.plot(dategrid,xmod_trans,color=color[i2],label=lab[i2])
   
   ax1.plot(dategrid,xmod,label=None,color=color[i])
   ax1.plot(dategrid,xmod,label=None,color=color[i])
   
   
   
   #ax1.plot(dategrid,xmodlo,label=None,color=color[i])
   #ax1.plot(dategrid,xmodhi,label=None,color=color[i])
   #ax1.fill_between(dategrid,xmodlo,xmodhi,alpha=0.25,label=None,color=color[i])
   dates_dat       = meta_df[i].values[:,0]
   #values_dat      = meta_df[i].values[:,1] 
   #dates_start     = dates_dat[0]
   dates_end       = dates_dat[-1]
   df_now          = meta_mod[i]
   #try:
   idfc = np.where(df_now[:,0] > dates_end)[0][0]
   #except:
   #idfx = -1
    
   ax1.fill_between(dategrid[idfc:],xmodlo[idfc:],xmodhi[idfc:],alpha=0.25,label=None,color=color[i])
   
   ax1.set_xlabel('dates')
   #ax1.set_ylabel(lab[i],rotation = 0)
   ax1.annotate(lab[i],(-0.1,0.5),xycoords='axes fraction',horizontalalignment='right')
    
   #if i = 0 add on the forecasts using all the other models
   if (i == 0):
    idlo = np.where(tgrid > tgrid[-1]-forecast_period)[0][0]
    ntgrid = np.shape(tgrid)[0]
    points_lo=[]
    points_med=[]
    points_hi = []
    
    
    sig_old = []
    sig_new = []
    for i2 in range(0,ntgrid,1):
     points = np.array([transform_mod[i3][i2,1] for i3 in range(nframes)],dtype='float')
     sigma  = np.array([transform_mod[i3][i2,2] for i3 in range(nframes)],dtype='float')
     r      = np.array([ft for ft in fitstat_cross])
     r[r<0] = 0#dont use -ve ccf info (technically should but fix later)
     mean = np.sum(points*r/sigma**2)/np.sum(r/sigma**2)#ccf weighted mean for forecast
     sd  = np.sqrt(1/np.sum(r/sigma**2))
     sig_old.append(sigma)
     sig_new.append(sd)       
     pl,pm,ph = mean-sd,mean,mean+sd#np.percentile(points,[15,50,85])
     points_lo.append(pl)
     points_med.append(pm)
     points_hi.append(ph)
     #print(i2,i,improvement_potential[i2,i])
     if (i2 == 0):
      for ilc in range(nframes):
       improvement_potential[:,ilc] = 1. - np.sqrt( 1. / ( 1. + r[ilc]    ) )#np.sqrt( 1./sigma[0]**2 / ( 1./sigma[0]**2 + r[i]/sigma[i]**2    ) )
    points_lo  = np.array(points_lo)
    points_med = np.array(points_med)
    points_hi  = np.array(points_hi)
    #ax1.fill_between(tgrid[idlo:],points_lo[idlo:],points_hi[idlo:],alpha=1.0,color=color[0],label='multivariate forecast')
    #ax1.plot(tgrid[idlo:],points_lo[idlo:],alpha=1.0,color=color[0],ls='--',label='multivariate forecast')
    #ax1.plot(tgrid[idlo:],points_hi[idlo:],alpha=1.0,color=color[0],ls='--',label=None)
    ax1.bar(dategrid[idlo:],points_med[idlo:],alpha=0.3)
    
    #for it in range(10):
    # print(dategrid[idlo+it],points_med[idlo+it],
    
    ax1.legend()
    sigma_mean = np.mean(np.array(sig_old))
    newsigma_mean = np.mean(np.array(sig_new))
    #replace the error bars on the driving series with the updated errors from the combined forecast
    points_lo = np.array(points_lo)
    points_hi = np.array(points_hi)
    points_sd = (points_hi - points_lo) / 2
    transform_mod[0][idlo:,2] = points_sd[idlo:]
    meta_mod[0][idlo:,2] = points_sd[idlo:]
    
  
   #set plot limits
   ymin,ymax = np.min(meta_lc[i][:,1]),np.max(meta_lc[i][:,1])
   ax1.yaxis.tick_right()
   if (i == 0):
    ymin = min(ymin,np.min(points_lo[idlo:]))
    ymax = max(ymax,np.max(points_lo[idlo:]))
    ax1.set_title(figure_label)
   yrange    = ymax - ymin
   ax1.set_ylim([ymin - yrange/10,ymax + yrange/10])
   if (i < nframes - 1):
    ax1.set_xticks([])
    
   ax1.set_xlim([datemin,datemax])
    
  ##bar plot of self-fit statistics
  #ax1 = plt.subplot(gs1[ndown-1,:2])
  #xpos = np.arange(nframes)
  #[ax1.barh([xpos[i]],[fitstat_self[i]],label=None,color='k') for i in range(nframes)]
  #ax1.set_title('correlation coefficient') 
  #ax1.set_yticks(xpos) 
  #ax1.set_yticklabels(lab)
  #ax1.tick_params(axis='y',rotation=0)
  
  #bar plot of cross-fit statistics
  ax1 = plt.subplot(gs1[ndown-1,3:])
  xpos = np.arange(nframes)
  [ax1.barh([xpos[i]],[accuracy_mape[i]*100],label=None,color='k') for i in range(nframes)]
  ax1.set_title('Mape accuracy metric') 
  ax1.set_yticks(xpos) 
  ax1.set_yticklabels(lab)
  ax1.set_xticks(np.arange(0,110,50))
  ax1.tick_params(axis='y',rotation=0)
  ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30)
  











 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 ##!!!!!output the model info into a dataframe to be converted to .csv spread sheet!!!!!#!!
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 
 writer = pd.ExcelWriter(figure_title[:-4]+'.xlsx',engine='xlsxwriter',datetime_format='dd mm yyyy')
 uncertainty_mod = []
 for i in range(nframes):
  title = lab[i]
  #add in point-to-point scatter to uncertainty estimate
  dates_dat       = meta_df[i].values[:,0]
  values_dat      = meta_df[i].values[:,1] 
  dates_start     = dates_dat[0]
  dates_end       = dates_dat[-1]
  df_now          = meta_mod[i]
  values_mod      = df_now[:,1]
  mape_now        = np.zeros(ntgrid)
  mean_now        = np.mean(values_mod[values_mod != 0])
  output_now      = np.zeros((ntgrid,4),dtype='object')
  for it in range(ntgrid):
   dates_mod = df_now[it,0]
   sd_now    = df_now[it,2]
   sd_tot    = sd_now
   if ((dates_mod > dates_start) & (dates_mod < dates_end)):
    idx          = np.where(dates_dat < dates_mod)[0]
    ilo          = idx[-1]
    ihi          = ilo + 1
    scatter      = np.std(values_dat)
    sd_extra     = np.std(scatter)
    sd_tot       = np.sqrt(sd_now**2 + sd_extra**2)
    df_now[it,2] = sd_tot
   mape_now[it] = sd_tot/mean_now
  
  output_now[:,:3] = df_now
  output_now[:,3]  = mape_now*100  
  df_output    = pd.DataFrame(output_now,columns=['Dates','model','model_range (uncertainty)','Median Absolute Percentage Error (MAPE)'])
  idfc = np.where(df_now[:,0] > dates_end)[0]
  
  if (len(idfc) > 0):
   df_output[:idfc[0]].iloc[:,:3].to_excel(writer,sheet_name = title,)
   df_output[idfc[0]:].to_excel(writer,sheet_name = title+' forecasts')
   ut = tgrid[idfc[0]:]
   x = df_output[idfc[0]:].values[:,1:3]
   ux,us = x[:,0],x[:,1]
   
   xhere  = np.array([ut,ux,us]).T
   uncertainty_mod.append(xhere)
   
   #print(np.shape(xhere))
   #print(len(uncertainty_mod))
   #input()
 writer.save()






 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 ##!!!!!output forecast figure!!!!!#!!
 #!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!
 nmonths = np.int(np.floor(forecast_period/30))
 nspace = 2

 xpos = np.arange(nmonths*(nframes+nspace))
 #fig = plt.figure()
 ysave = []
 for i in range(nframes):
  mm = np.array(uncertainty_mod[i])
  tm = mm[:,0]
  dtm = tm-tm[0]
  y = []
  values_dat      = meta_df[i].values[:,1] 
  mean = np.mean(values_dat[values_dat !=0])
  for im in range(nmonths):
   id = np.where(dtm > im*30)[0][0]
   pred = mm[id,1]
   sig = mm[id,2]
   pc = sig/mean
   y.append(pc*100)
  ysave.append(y)
 
 ax1 = plt.subplot(gs1[ndown-1,:2])
 [ax1.bar(xpos[i::nspace+nframes],ysave[i],label=lab[i],color=color[i]) for i in range(nframes)]
 lab_r = np.array(['']*len(xpos))
 idskip = np.int(np.floor(nframes/2))
 idr = np.arange(idskip,len(xpos),nframes+nspace)
 labs = [np.str(i) for i in range(nmonths)]
 for ir in range(nmonths):
  lab_r[idr[ir]] = labs[ir]
 

 ax1.set_xticks(xpos)
 ax1.set_xticklabels(lab_r)
 ax1.set_ylabel('MAPE')
 ax1.set_xlabel('forecast month')
 ax1.legend(fontsize='xx-small')
 ax1.tick_params(axis='x',which='both',length=0)
 
 if (figure_title !='show'):
  plt.savefig(figure_title) 
 else:
  plt.show() 
 
 


 return(fitstat_self,fitstat_cross,sigma_mean,newsigma_mean,
 accuracy_envelope,accuracy_variance,accuracy_mape,uncertainty_mod,improvement_potential) 
 
#return(parmout,covout,freq,tplot,xplot,xplotsd,xplotlo,xplothi,p0_out,w0,dw,sig2_prior,xplotsave)



##test make fake data
#import vortexa_makefake as mf
#tlo = 0
#thi = 100.0
#dt  = 1.0
#dfin = mf.makefake(tlo,thi,dt,shift = [10],noise=0.1,iseed=31234)
#
#
#a=mass_mod_test(dfin,dateref=-1,forecast_period = 50,figure_title='',model='poly',dtgrid=0.1,
#transform_model='convolve',convolve=[[0,0.1],[30,10.],[60,10]])






























#plot x and y using generalised poly fit
def reg_plot(xlist,ylist,lab = '',xlabel='',ylabel='',fig_title='reg_plot.pdf',nlag = 5):
 
 plt.clf()
 ny = len(ylist)
 rcoef = [0.0]*(ny)
 
 
 #try:
 ndown = 4
 nalong = ny 
 
 if (lab == ''):
  ax_title = ['']*ny
 else:
  ax_title = lab
 if (ylabel == ''):
  ylab = ['']*ny
 else:
  ylab = ylabel
 if (xlabel == ''):
  xlab = ['']*ny
 else:
  xlab = xlabel
  
 ngrid = 1000

 yglo = []
 yghi = []
 xgs = []
 ordersave = []
 lagsave = []
 ccfsave = []
 
 idx = 0
 for y in ylist:
  x = xlist[idx]
  idinclude = np.where((x > 0) & (y > 0))[0]
  
  
  
  
  #perform ccf 
  ccf = ss.correlate(x,y)
  nccf = np.shape(ccf)[0]
  neg = np.int(np.floor(nccf-1)/2)
  lag = np.zeros(nccf)
  idlag = np.arange(neg+1)
  lag[:neg] = idlag[-1:0:-1]
  lag[neg+1:] = idlag[1:]

  #identify lag that maximizes ccf only shift if abs(lag) < nlag maximum number of months
  lagmax = np.argmax(ccf)
  if (np.abs(lagmax) > nlag):
   lagmax = 0

  lagmax = 0
  kernel = np.zeros(nlag*2 + 1)
  kernel[nlag + lagmax] = 1

  ynew = apc.convolve(y, kernel)
  print(y)
  print(ynew)
  print()


  
  
  lagsave.append(lag)
  ccfsave.append(ccf)
  
  
  
  
  
  
  
  
  
  xmin,xmax = np.min(x),np.max(x)
  xrange = (xmax - xmin)/10
  xgrid = np.linspace(xmin-xrange,xmax+xrange,ngrid)
  xgs.append(xgrid)
  
  a = vpf.fit_search(x[idinclude],y[idinclude],maxorder=4,xgrid=xgrid)
  order_cisqred, order_aic, order_bic = a
  order = order_aic 
  ordersave.append(order)
  yg_med,yg_lo,yg_hi,cov,cisq,cisq_red,bic,aic,rmd = \
  vpf.fit(x[idinclude],y[idinclude],order=order,xgrid=xgrid,confidence=0.01,nits=20000,figure_title='')
  #correlation coefficient
  rcoef[idx] = np.corrcoef(x,y)[0,1] 
  yglo.append(yg_lo)
  yghi.append(yg_hi)
 #make figure
  idx = idx + 1
 
 gs1 = gridspec.GridSpec(ndown,nalong)
 gs1.update(left=0.1,right=0.98,wspace=0.4,hspace=0.0,bottom=0.2,top=0.9)
 for i in range(ny):
  xgrid = xgs[i]
  xmin,xmax = np.min(xgrid),np.max(xgrid)
  y = ylist[i]
  yg_lo = yglo[i]
  yg_hi = yghi[i]
  ymin,ymax = np.min(y),np.max(y)
  yrange = (ymax - ymin)/10
  ax1 = plt.subplot(gs1[0:3,i])
  ax1.scatter(x[idinclude],y[idinclude],color='k')
  ax1.plot(xgrid,yg_lo,color='b',label='fit')
  ax1.plot(xgrid,yg_hi,color='b',label=None)
  ax1.fill_between(xgrid,yg_lo,yg_hi,alpha=0.3,color='b',label=None)
  ax1.set_xlabel(xlab[i])
  ax1.set_ylabel(ylab[i])
  ax1.set_title(ax_title[i])
  ax1.set_xlim([xmin,xmax])
  ax1.set_ylim([ymin - yrange,ymax + yrange])
  ax1.annotate(r'$r_c = '+np.str(np.round(rcoef[i],2))+'$',(0.9,0.9),xycoords='axes fraction',horizontalalignment='right',color='b') 
  ax1.annotate('order = '+np.str(np.int(ordersave[i])),(0.9,0.8),xycoords='axes fraction',horizontalalignment='right',color='b') 
  
  
  
  #plot ccf figure only include first nlag time steps
  ax2 = plt.subplot(gs1[3,i])
  lag = lagsave[i]
  ccf = ccfsave[i]
  nccf = len(lag)
  nl0  = np.int(np.floor(nccf/2))
  ilo = nl0 - nlag
  ihi = nl0 + nlag
  ax2.bar(lag[ilo:ihi],ccf[ilo:ihi])
  ax2.set_xlabel('lag')
  ax2.set_ylabel(r'$r_c$')
  
  
  
  
  
  
  plt.savefig(fig_title)
   
 #except:
 # pass
 

 
 
 
 return(rcoef)



 
 
 





















