#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!
#configuration for the fake data plotting routine for paper in myplot_results.py#!!!!!!!!!
#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!

dirs = ['/Users/ds207/Documents/standrews/sta/fort/fortcode/lucky_results/sdss_fake_new3/output_20180125_001',
'/Users/ds207/Documents/standrews/sta/fort/fortcode/lucky_results/sdss_fake_new10/output_20180125_001',
'/Users/ds207/Documents/standrews/sta/fort/fortcode/lucky_results/sdss_fake_new19/output_20180125_001']
truths = [np.log10(0.766459),np.cos(30.*np.pi/180),0.75]
col = ['r','k','b']
binsin = 20
fparms = 'outputpars.dat'
emlog  = 8.4971866218073941
bfrac = 3./4#10./12#2./3
idmbh = 2
iddeg = 3
idslope = 4












#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!
#!!!!!!!!!!!configuration for the fake data plotting routine in py_plot_fakedata.py#!!!!!!!!!
#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!#!!!!!!!!!

idburnin = 3./4

#light curve ascii files from which to calculate snr (x-axis of plot)
filesnr=['','','']

#light curve ascii files from which to calculate the inferred parameters (mdot, inc, trslope)
fileparm=['','','']


#files saving the response function
fileresp_mod = ['','','']
fileresp_sig = ['','','']

wavlag = ['g','i']
lagtrue = [1.0,1.3]


#indexes of columns in filedat parameter correspondong to each disk parameter
idparm = [1,2,3]
labparm = ['mdot','inc','slope']
