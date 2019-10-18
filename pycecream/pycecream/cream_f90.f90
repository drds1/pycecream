!! NOTES AND BIG FIXES BELOW
!! update: 05/01/2018 include "shareos.par" file that tells cream to fix the offset difference and step together the offset parameters of pairs of light curves (new arrays shareos_diff(nlc),ilc_shareos_a(nlc),ilc_shareos_b(nlc),nlc_shareos)
!! update: 18/08/17 cosinc prior now gaussian prior on cosinc with mean inclination and sd inclination specified in creamin.par
!! NOTE:  30/04/2015 CREAM now runss a little slower at high redshift since the ookback time has increased (i.e the vertical line thrugh the top plot has moved to the right because a 10 day delay in the emitted frame corresponds to a 20 day lookback at redshift 1.
!! update 29/04/2015 main file. Ntaugrid=ceiling((maxtau*redshiftadd1)/dtgrid) redshiftadd1 added to increase number of evaluated points in transfer function because at high z, more points are needed to get to the same max tau value with points that are shrunk by 1+z
!! update 18/12/14 Problem with psisave not being properly updated after accept step fixed
!! .....2) 		   altertf introduced logical statement rather than if ip = npcosincidx .or. ... .or.   (easier to read)
!! update 12/12/14 floor on scaling of inclination, not allowed to drop below 1.e-3.
!! update 22/09/2014 changed 'f' in mcmcmulti_nobreak.f90 when deciding wlo. Now set to 1.0 from 0.5
!! update 14/08/14, if dtaugrid set to -ve in par file, will chose dtgrid=dtaugrid = 0.5 * 2 pi/whi
!! update 14/08/14, stop annotating with +/- if the inclination or mmdot is fixed
!! update 13/8/14, now if tf's go to zero before taugrid(ntaugrid)*taugrid_change_scale, scale the tau axis of the plots to better show the transfer function. Also no longer plot the y axis numbers
!! update 10/8/14 bug fix isubburnin = iteration - iburnin2 NOT iteration - iburnin2 +1
!! update 8/8/14 bug fixed  NPebmvidx = NPsigexpand+NPsigexpandidx+1 in mcmcmulti_nobreak.f90
!! update 8/8/14 in psres(iw)= ....... dw changed to df = dw*0.5/pi fix model being above points on power spectrum plot
!! update 17/7/14 program changed so that when adjusting cos(inc), a value not gt 0 and lt 1 is rejected rather than recast until it is
!! update 18/7/14 bug fix where echo lightcurve was being saved incorrectly.
!! - also introduced firstpsi (firstpsi is true whenever the routine is called until the first parameter is tested usefull for
!! -- forcing the code to accept the first step). Note not the same as firstrun, firstrun is the same but only true if iteration =1
!! firstcall is true whenever the subroutine is called until the subroutine ends, (i.e true for a whole run through of all parameters).
!! - also, set iteration in mcmcmulticont to iteration plus 1 (so that we are not repeating the same iteration).

!! code is now set to accept first move when mcmcmulticont.f90 is called
!! bug fix 16/09/2014 now removes old backups properly i.e keeps only nbackups of backups

!! update 29/09/2014 P(NPfurconstidx) (the offset term in the driving light curve) no longer scales. This removes the degeneracy between
! this term and the offset term in the echo light curve.

!update 6/10/14
!! added a whi term into the par file abov NW to allow better customisation of start and end frequencies
!! changed f to 0.5 when calculating number of frequencies

! 19/10/14 plots now normalise so that the yaxis  always goes from 0 to max of each tf for each lc.

! 20/10/14 New 'pubplot' feature. Chops off the pspec and bof plot, makes the others bigger, adds time axis to bottom echo plot and makes plot edges bolder.
!! program to read in a set of light curves and run monte carlo fits including the time delay parameter

! 31/10/14 New introduced iburnin(Nburnin) array that has multiple stages of burnin and resets the error envelopes each time a new stage is triggered
!!!! optional apply tf scalings rather than re-ealuate (ONLY FOR FLAT DISK WITH T(R)^3/4)




module creammod
integer NLC,NT,NTgrid,Ntaugrid,NPF,NPpspec,NPpspecidx,NPscaleidx,NPscale,NPdelay,NPumbhidx,NPcosincidx,NW,NP,NPdel,NPdisk,NPtr,&
NPtridx,NPxray,NPxrayoffsetidx,Nxray,interpidxmin,interpidxmax,NPfurconstidx,NPsigexpand,NPsigexpandidx,&
NPebmvidx,NPebmv,Nstore, NPvarexpand, NPvarexpandidx, NPpoly,NPpolytot,idxpoly,NPpolyidx,&
idaffine2,pAT_affine, pRT_affine

integer NPgalidx,NPdlidx,NPgal,NPoffset,NPoffsetidx,nr,version,AT,RT,iburnin1,iburnin2,iplotsave,&
NPth, NPthcentidx, NPthfwhmidx, nofur, pricream_idxnow, Nkeep_skipcor, idref_skipcor,&
iaffine_count,nlc_shareos

integer,allocatable:: lo(:),hi(:),interpidx(:),xrayinterpidx(:),taugrididx(:),iburnin(:),&
ilc_shareos_a(:),ilc_shareos_b(:),pRT_save(:)

integer nbackuptimesiplotsave,nbackup,NWres,Nburnin,iburninnow,iburninprev,ibnow,&
it_BIG_atstart, idxtaulo, idxtau0, npricream, ipnow, nparfix, iosfix, idaff_sigexp,&
idaff_varexp,n_affine,i_affinestart,iaffinenow,idaffine,itmin_sigrej,&
iat_ann, i_tfx,idlim


real,parameter:: pi=3.14159265,rttwopi=sqrt(2*pi),twopi=2*pi
real dtaugridtemp,urin,urout,alphamean,w0mean,p0mean,f0mean,diskk,diskka,alb,eff,&
alphasig,meanlogw0,stretchscale,ur0,uhx,umdot,sigalphasquare,dw,dtgrid,xraymed,sigp0square,&
sigw0square,bofmed,furconst,siglogp0,cisqxray,ebmvmw,ebmvagn,ebmvagnsteplog,redshift,&
tgridloextra,tgridhiextra, embhref, wavref, ddeginc, betamean, betasig, sigbetasquare,&
dtgrid_em,pri_mmdotmean,sigmmdotsquare,alpha_cosinc, cosinc_cut,redshiftadd1, omega_m,&
omega_l, taugridlo, sigcosinc2, valcream_now, parfixnow,umdotref,stepchange_affine,&
frac_stepchange_affine,T_ann, frac_ann,T1vnorm, T1vnorm_scale,T1inorm, T1inorm_scale,&
plonow,phinow, varerln, varerlnxray, varerlnold,varerlnold_affine,&
varerlnxrayold,varerlnxrayold_affine, sdprior


real,allocatable::xxray(:),txray(:),erxray(:),cisqplot(:),xgridold(:),cisqnew(:),fracalongxray(:),&
erxrayvar(:),ervar(:),erxrayvarold(:),ervarold(:),&
ervarold_affine(:),erxrayvarold_affine(:), bg_now(:),p_affine(:)
real,allocatable:: fdisk(:),xgridplot(:,:),xgridplotsave(:,:),&
xgridplotsave_affine(:,:),bofsave(:),bofrejectnew(:),bofreject(:),&
shareos_diff(:)
real,allocatable:: swgrid(:,:),cwgrid(:,:),bofrejectbef(:),psisave(:,:),&
psisave_affine(:,:),xgridop(:),w(:),&
fres(:),psres(:),wavem(:),fdisksave(:),fdisksave_affine(:),&
greytr(:),bofme(:), psistore(:,:),xgridop_em(:), ave_overflow(:),sd_overflow(:)
real,allocatable:: wavobs(:),fvf(:),fvb(:),fursum(:),fursumold(:),sigexpandparm(:),&
sigexpandsteplog(:),respcent(:),respsig(:), varexpandparm(:), varexpandsteplog(:),&
thtempcent(:), thtempfwhm(:),thtempcent_scale(:), thtempfwhm_scale(:),&
thtempfwhm_pri_width(:),thtempcent_pri_width(:),thtempfwhm_pri(:),&
thtempcent_pri(:), sigexpandpri_cent(:), sigexpandpri_fwhm(:),&
polyin(:),polyin_step(:,:),polyin_pri(:,:),senv1_bg(:,:),senv2_bg(:,:),&
bg_save(:,:),polyin_pri_mean(:,:),tref_poly(:), sigrej_k(:), covmat(:,:),&
parsave_affine_T(:,:),psinormold_affine(:),xgridold_affine(:),fursumold_affine(:),&
covmat_affine(:,:),echosave(:,:),sigechosave(:,:)

double precision,allocatable:: senv1(:,:),senv2(:,:),senv1drive(:),senv2drive(:),&
senv1tf(:,:),senv2tf(:,:),senv1fur(:),senv2fur(:),sumpar(:), covsum(:,:),sum2par(:),&
senv1p(:),senv2p(:)
real, allocatable:: parmean(:), sigpar(:),parsave(:,:), sdparm(:), &
avparm(:),psimeansave(:),pricream_mean(:), pricream_sd(:),&
pricream_bofold(:), pricream_bofnew(:), plo(:), phi(:)

double precision, allocatable:: senv1offset(:),senv2offset(:)
double precision senv1furconst,senv2furconst,senv1mmdot,senv1cosinc,b1,b2
double precision, save:: sumcosinc=0.d0, summmdot=0.d0,sum2mmdot=0.d0,sum2cosinc=0.d0,&
sumtrirad=0.d0,sum2trirad=0.d0,sumtrvisc=0.d0,sum2trvisc=0.d0
integer,allocatable:: isharelc(:), pricream_idx(:), sigrej_mode(:), sigrej_mark(:),&
sigrej_mark_save(:), ip_affine(:)

logical,allocatable:: pversion_affine(:)


character(len=256)filename2,filenametemp2,cwd,fname_savelcl,tit_save
character(len=256),allocatable:: file_lc_in(:)

character(1000) dirworking, dirmain
logical noconvolve,yesxray,defaults/.True./,plotsave,plotscreen,positivity,break,sigexpand,&
savelightcurves/.True./,greyscale/.true./, lstore/.True./, pubplot, BIGsave,dtauskip/.True./,&
tfcheat/.False./, varexpand, samescale/.false./,bgvary /.false./, backupon /.false./,&
pythonplot /.true./, sleepon /.false./, pricream /.false./, sigrej /.false./,&
skipcor /.false./, affine/.false./, chaos/.false./, autoacc/.false./, stepdiag /.false./,&
quick_conv_itp /.true./

integer, save:: idxskip = 0, Ntau_hr_in = 10


contains

!!!!!!!!!!! SUPERFITDATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine to read light curve filenames from an input text file, open up the files and read the data into t,x,and er arrays

!! inputs::: Nignore the number of lines of text at the start of the ascii file which should be ignored

!! outputs, t,x,er the allocated and filled arrays containing the light curve info
!Npoints:: and iteger array indicating the start and end indicees of the x,t, and er arrays for each light curve

subroutine superfitdata(t,x,er,Npoints,Nignore)

implicit none

integer Nres,Nex,i,nits,Nignore,iignore
character*10 reply
character(200) filename
real val1,dw
integer, allocatable:: Npoints(:)
real, allocatable:: x(:),er(:),t(:)
!integer Nt,NLC,NPF,NP,NW,Ntgrid,Ntaugrid
integer count,it,numline,a

!common /Nums/ Nt, NLC,NPF,NP,NW,Ntgrid,Ntaugrid


!! read the light curves from the files. The file names must be listed in the fitlcnames.dat file
!       NLC = numline('fitlcnames.dat')
!call chdir('./input') !! change to the directory with all the input information
open(2,file='creamnames.dat')
allocate(Npoints(NLC+1))
Npoints(1)=0
!! get the total number of data points added over all the light curves and the number of points in each file
Nt=0
write(*,*)'reading',NLC,'light curves from fitlcnames.dat'

if (yesxray) read(2,*)					!! skip the xray data if present
do i =1,NLC

read(2,*) filename


Npoints(i+1)=numline(trim(filename))+Npoints(i)-Nignore-1
!       write(*,*) trim(filename),numline(trim(filename))
write(*,*) 'file:  ',trim(filename),Npoints(i+1),'points'

!       Ntot=Ntot+Npoints(i+1)
enddo
close(2)
Nt=Npoints(NLC+1)
allocate(t(Nt))
allocate(x(Nt))
allocate(er(Nt))



open(2,file='creamnames.dat')
if (yesxray) read(2,*)						!!!!  skip the xray data if present
count=1 ! store the array element
do i = 1,NLC
read(2,*) filename
open(3, file=trim(filename))
iignore=0
do while (iignore .lt. Nignore)
read(3,'(A)')
iignore=iignore+1
end do
do it=1,Npoints(i+1)-Npoints(i)
!     write(*,*) count,NLC
read(3,*) t(count), x(count), er(count)
count=count+1
enddo
close(3)
enddo
close(2)
!call chdir('../') !! change back to the main directory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine





!!!!! New 20th march
!!!!INPUTS
!! subroutine to interpolate the input t(N), and x(N) grid with error sig(N)
!takes the median error and sets an extended array ergrid to haev the median of er
! inputs are the t(N),x(N),er(n),w(NW),p(2*NW),tgridlo,dtgrid
! outputs are each of the p(2*NW) terms see /fortcode/furoptscal.f90













!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superfitplot(tgrid,tgridminidx,tgridmaxidx,xgrid,&
taugrid,y,Npoints,t,x,err,p,pscale,psafe,w,iteration)

integer iteration,tgridmaxidx,tgridminidx!NLC, Ntgrid, Nt, NP,ic
real tgrid(NTgrid),x(Nt),t(Nt),err(Nt),p(NP),taugrid(Ntaugrid)
logical psafe(NP)
real xgrid(NTgrid),itarray(iteration),pscale(NP)
real miny,maxy,minx,maxx,gap,upy,loy,xgridmax,xgridmin
real frac,fracsub1,dy,moreplotspace,flatlinex(2),flatliney(2),med
integer county,Npoints(NLC+1),nplots
character*1000 title,plotinfo,tit,tit2,fname_savelc
character*10 txraymed,twav,ttaucent,titN,ttaucentlo,ttaucenthi,tse,bofmin_text_a
character*100 nam,nam2,fname_tit

integer it,iw,i_done_this_before
real miny0,maxy0,maxdelay,w(NW),deltabofsave
real,save:: maxybof,minyps=-3.0,maxyps=10,xraymean
real vertlinex(2),vertliney(2),y(Ntaugrid,NLC),xw0(2),yw0(2),pspecpt(NW)
integer ier,pgopen
character(len=100) temp1,temp2,temp3,temp4,tdel1,tdel2,tcisq,temp5,temp5step,&
tsigexpand,temp4sd,temp3sd,temp6,temp6sd, temp5sd, bofmin_text
real bl0,br0,bl1,taumax,tgridmax,tgridplotmin,tgridplotmax
logical, save:: firstrun =.True.
real xgridplotcent(Ntgrid),xgridsd,xgridenvlo(Ntgrid),xgridenvhi(Ntgrid),&
xgridenvlodrive(Ntgrid),xgridenvhidrive(Ntgrid),xgridmean(NTgrid)
real tfplotcent(Ntaugrid), tfenvlo(Ntaugrid), tfenvhi(Ntaugrid),taumean(NLC),&
pspecsd,pspeccent(NWres),pspeclo(NWres),pspechi(NWres),rt
real tfsave(Ntaugrid,NLC),&
xdrivesave(NTgrid), bofsaveplot(iteration),&
sigxdrivesave(NTgrid),sigtfsave(Ntaugrid,NLC),&
offset_bg(Ntgrid), offset_bg_sd(Ntgrid),bgop_save(Ntgrid,NLC), bgopsig_save(Ntgrid,NLC)
logical posterplot /.False./, plotergs /.False./, aas_lc /.True./, showchisq /.True./,&
shade /.False./, reddotslast /.False./, scale_y2fit_mod /.False./, showsig /.False./,&
showsigbars /.True./, showres /.False./, baron /.False./

integer,allocatable:: cusnumlo(:), ilc_1page(:), cusnumhi(:),ikey_temp(:),&
ilc_temp(:),ilc_tempold(:)
integer cuscol(NLC), NLC_1page
logical CREAM_plot_ex, wav_ann_input /.False./, cuscol_input /.False./, cusnum_input /.False./,&
psi_echo_lab_input /.False./, one_page_plot /.False./
character(len=100) t_xaxlab, dlclab, psilab(NLC), echolclab(NLC)
character(len=10) wav_ann(NLC)
real,allocatable,dimension(:):: tf_cheatlo,tf_cheathi,tf_cheatmid,taucheat,taumax_man,&
xres,xmod_itp,x_temp,t_temp,errtemp




!initialize required arrays to 0
do ilc = 1,NLC
sigechosave(1:Ntgrid,ilc) = 0
echosave(1:ntgrid,ilc) = 0
tfsave(1:Ntaugrid,ilc) = 0
sigtfsave(1:Ntaugrid,ilc) = 0
enddo

!reddots last replots datpoints on top with big red dots for visual paper aid.
! scale_y2fit_mod changes plotlimits of driving light curve in case model has higher amplitude variations than data (arbitrary i.e *10 not optimal)

!! options which specify (taugrid_change_scale) when to change from plotting delay function on same scale as rest of plot to enhanced scale and where to extend expanded scaled plots to (what fraction of psimax do we cut off the plots0
taugrid_change_scale=5.0 !! if taugrid_change_scale taugrid / taugrid_change_scale
!if taugrid(idxtau0) .lt. taugrid(Ntaugrid)*taugrid_change_scale) then  !! taugr when to change transfer function plotting scale

if (sigexpand) then
icsigexpand=5
endif







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! axis labels, showchisq, showsig, units 3/9/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
t_xaxlab='Time / days'
dlclab =  '\gDF\dx\u  (Arbitrary Units)'
psilab(1:NLC) = ''
!psilab(NLC/2) = '\gq(\gt|\gl)  (Arbitrary Units)'
echolclab(1:NLC) = ''
!echolclab(NLC/2) = 'F\d\gn\u(\gl)  (Arbitrary Units)'
inquire(file='../../CREAM_plot.dat',exist=CREAM_plot_ex)
if (CREAM_plot_ex) then
open(unit=1,file = '../../CREAM_plot.dat')
read(1,*)
read(1,*) t_xaxlab
read(1,*) dlclab
read(1,*) !description bit
read(1,*) psi_echo_lab_input,Nplot_tot
if (psi_echo_lab_input) then
do ilc = 1,Nplot_tot
read(1,*) psilab(ilc)
read(1,*) echolclab(ilc)
write(*,*) psilab(ilc), echolclab(ilc), ilc
enddo

endif
read(1,*)
read(1,*) showchisq, showsig, plotergs, showsigbars
read(1,*)
read(1,*) cuscol_input
if (cuscol_input) then
read(1,*) cuscol(1:NLC)
endif
read(1,*)
read(1,*) cusnum_input, icusnum_nplot
if (cusnum_input) then
allocate(cusnumlo(icusnum_nplot),cusnumhi(icusnum_nplot),taumax_man(icusnum_nplot))
read(1,*) cusnumlo(1:icusnum_nplot)
read(1,*) cusnumhi(1:icusnum_nplot)
read(1,*) taumax_man(1:icusnum_nplot)
endif
read(1,*)
read(1,*) one_page_plot, NLC_1page


if (one_page_plot) then !!new 14th october plot a single plot with a user specific set of light curves displayed
allocate(ilc_1page(NLC_1page))
read(1,*) ilc_1page(1:NLC_1page)
else
read(1,*)
endif					  !!

read(1,*)
read(1,*) wav_ann_input
if (wav_ann_input) then
do ilc = 1,Nplot_tot
read(1,*) wav_ann(ilc)
enddo
endif
close(1)
endif
if (wav_ann_input .eqv. .false.) then
do ilc = 1,NLC
write(wav_ann(ilc),'(I10)') int(wavobs(ilc))
enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


tmax=tgrid(Ntgrid)
tmin=tgrid(interpidxmin)!tmin-tgridloextra
tstart=tgrid(1)
taumax=taugrid(Ntaugrid)
tgridplotmax=tgrid(NTgrid)
tgridplotmin=tmin-taumax!-tgridloextra!tgrid(1)

!!!!!!!!! plotting information
bltf0 = 0.1 ! new 19/10/2014 where the transfer function plots start (sometimes these are squashed and too small)
bl0=0.2    ! left and right x axis limits
if (plotergs) then
br0 = 0.88
else
br0 = 0.88
endif

!tgridplotmax=tmax
!tgridplotmin=tmin-taumax
bl1=bl0 + (tmin-tstart)/(tmax-tstart)*(br0-bl0)!*(taumax/redshiftadd1)/(tgridplotmax-tgridplotmin)+bl0!bl0+taumax/(taumax+(tgridplotmax-tgridplotmin))*(br0-bl0)
!bl1 = tmin -

!write(*,*) tmin
!write(*,*) tstart + bl1*(tmax-tstart)
!stop

br1=br0!bl1+(tmax-tmin)/(tgridplotmax-tmin) * (br0-bl1)  !! help help
!! just to properly partition the x and y plots above
!!!!!!!!!!
!write(*,*) bl0,bl1,br0
!write(*,*) taumax*redshiftadd1, tmin,tmax
!stop



upy=0.87    !!! set as 0.90 to leave room for plotting information bar at the top of the screen

!! information bar
bofmin=minval(bofsave(1:iteration))
write(temp1,"(I10.2)") NP!/ deltabofsave
deltabofsave=bofmed-bofmin
write(temp2,"(F10.2)") deltabofsave

!if (sigexpand) then
!write(temp5,'(F10.2)') f
!write(temp5step,'(F10.2)') pscale(NPsigexpandidx)
!temp5=adjustl(temp5)
!temp5step=adjustl(temp5step)

!endif
itsubburn=iteration-iburninprev
iburnin2 = iburninprev





if (iburnin2 .eq. 0) then  !! if burnin finished write the mean 16th may 14
write(temp3,"(F10.2)") acos(P(NPcosincidx))*180/pi
write(temp4,'(F10.2)') alog10(p(NPumbhidx))

write(temp5,'(F10.2)') p(NPtridx)
write(temp6,'(F10.2)') p(NPtridx+1)

else



rln10 = alog(10.0)
radincplot=avparm(NPcosincidx)

avcosinc = avparm(NPcosincidx)
sdcosinc = sdparm(NPcosincidx)
sdinc    = sdcosinc/sin(radincplot)!abs( 1/sin(radincplot) ) * sdcosinc *180/pi
deginc   = acos( avcosinc / (1.+0.5*sdinc*sdinc) ) * 180/pi

apmmdot = avparm(NPumbhidx)
apmmdotlog = alog10(apmmdot)

sumtrvisc = avparm(NPtridx)
sumtrirad = avparm(NPtridx+1)

sdmmdot=sdparm(NPumbhidx)!sqrt( abs((sum2mmdot - summmdot*summmdot/itsubburn)/(itsubburn) ) )
summmdot = apmmdotlog + 0.5/rln10*sdmmdot*sdmmdot/(apmmdot*apmmdot)
sdtrvisc = sdparm(NPtridx)!sqrt( abs((sum2trvisc - sumtrvisc*sumtrvisc/itsubburn)/(itsubburn) ) )
sdtrirad = sdparm(NPtridx+1)!sqrt( abs((sum2trirad - sumtrirad*sumtrirad/itsubburn)/(itsubburn) ) )

!write(*,*) 'fuck up checking...', sum2cosinc, sumcosinc, itsubburn, degincplot, iburninprev, iburnin2
!read(*,*)
! heelp
write(temp3,"(I0)") int(deginc)
if (pscale(NPcosincidx) .gt. 0) then
temp3=trim(temp3)//'\(2233)'

if (sdinc .gt. 1.) then
write(temp3sd,'(I0)') int(sdinc)
else
write(temp3sd,'(F10.1)') sdinc
endif


temp3sd=adjustl(temp3sd)
temp3=trim(temp3)//trim(temp3sd)
endif



!a = sumalog10(real(summmdot)/itsubburn)

!write(*,*) sdmmdot, summmdot,apmmdot, 'kjffjsdkjfdskjfs'
!read(*,*)
write(temp4,'(F10.2)') summmdot
if (pscale(NPumbhidx) .gt. 0) then
temp4=trim(temp4)//'\(2233)'
write(temp4sd,'(F10.2)') sdmmdot/(apmmdot*rln10)
temp4sd=adjustl(temp4sd)
temp4=trim(temp4)//trim(temp4sd)
!write(*,*) real(summmdot)/itsubburn,sdmmdot,trim(temp4),trim(temp4sd), 'see line 345'
!stop
endif

write(temp5,'(F10.2)') real(sumtrvisc)/itsubburn
if (pscale(NPtridx) .gt. 0) then
temp5=trim(temp5)//'\(2233)'
write(temp5sd,'(F10.2)') sdtrvisc
temp5sd=adjustl(temp5sd)
temp5=trim(temp5)//trim(temp5sd)
endif
write(temp6,'(F10.2)') real(sumtrirad)/itsubburn
if (pscale(NPtridx+1) .gt. 0) then
temp6=trim(temp6)//'\(2233)'
write(temp6sd,'(F10.2)') sdtrirad
temp6sd=adjustl(temp6sd)
temp6=trim(temp6)//trim(temp6sd)
endif
!write(*,*) real(summmdot)/itsubburn,sdmmdot,trim(temp4),trim(temp4sd)
endif


temp1=adjustl(temp1)
temp2=adjustl(temp2)
temp3=adjustl(temp3)
temp4=adjustl(temp4)
temp5=adjustl(temp5)
temp6=adjustl(temp6)


!'   bofdif='//trim(temp2)&//
!if (sigexpand) then
!plotinfo='NP='//trim(temp1)//'   f='//trim(temp5)//' step='//trim(temp5step)//&
!'    inc='//trim(temp3)//'   log\d10\u(M\dBH\u/M\d\2281\u): '//trim(temp4)//'   Iteration: '//trim(tit)
!else
plotinfo='N\dP\u='//trim(temp1)//'  i\uo\d='//trim(temp3)//&
'   MM\u\b\.\d (10\u7\dM\d\(2281)\u\u2\dyr\u-1\d): 10\u'//trim(temp4)//&
'   \d \ga='//trim(temp5)//'   \ga (irad)='//trim(temp6)
!endif
!!
extrafracbof=1.0 !! go a different amount above the range for the bof plot
extrafrac=0.3   !!! what fraction of the range to extend plot above and below.
loy=0.05
dy=upy-loy
county=1
gap=0.02
nrows=3+NLC






if (iburnin2 .eq. 0) then !! color of delay and tf's (changes after convergence)
ic=2
else
ic=4
endif

icpt=2
icnorm=1
icenv=4




ilc_lo = 1
ilc_hi = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Work out how many plots to put on each page. New sep 7 !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (cusnum_input) then
Nlcplot = icusnum_nplot
else
Nplots = 0
idxcount = 0
ilc_lo = 1
ilc_hi_temp = NLC
do ilctemp = 1,NLC
if (isharelc(ilctemp) .eq. 0) idxcount = idxcount + 1
if (idxcount .eq. 7) ilc_hi_temp = ilctemp
enddo
ilc_hi = ilc_hi_temp

NLCmax=min(7,idxcount)
Nlcplot = int(max(1.0,real(ceiling(real(idxcount)/NLCmax)))) !!! How many plots to make
endif



!! pub plot feature is new 20/10/14 doe not plot the pspec and bof plot
if (pubplot) then
nplotsreal =  1
else
nplotsreal=2
endif
!!
!nplotsreal = nplots

nplots = nplotsreal + NLCmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (one_page_plot) then
NLCplot = 1 !only make one plot if we have the user specific 'one_age_plot' option turned on new 14 october 2015
nplots = 1 + NLC_1page
endif


idxplot_cumulate = 1

do ilcplot = 1,NLCplot !! new 17th july if more than 7 light curves plot on separate plots
!nplots = nplotsreal
!if (ilcplot .gt. 1) nplots = nplots - 1

!!!!! Work out how many plots to make

if (cusnum_input) then
if (ilcplot .lt. NLCplot) then
ilc_lo_next = cusnumlo(ilcplot+1)
else
ilc_lo_next = cusnumlo(ilcplot)
endif
ilc_lo = cusnumlo(ilcplot)
ilc_hi = cusnumhi(ilcplot)
else

ilc_lo_next = min(ilc_hi+1,NLC)
endif




idxcount = 0
do ilctemp = ilc_lo_next,NLC
if (isharelc(ilctemp) .eq. 0) idxcount = idxcount + 1
enddo





if ((ilcplot == NLCplot) .and. (ilcplot .gt. 1) .and. (cusnum_input .eqv. .false.)) then
ilc_hi = NLC
endif
ilc_hi_temp = NLC
idxcount = 0
do ilctemp = ilc_lo,ilc_hi
if (isharelc(ilctemp) .eq. 0) idxcount = idxcount + 1
if (idxcount .eq. 7) ilc_hi_temp = ilctemp
enddo
ilc_hi_next = ilc_hi_temp
if (one_page_plot .eqv. .false.) nplots = idxcount + nplotsreal





write(tit2,'(I10.5)') ilcplot
tit2=adjustl(tit2)//'_'
write(tit,'(I10.5)') iteration
tit=adjustl(tit)
!!! if this is the first plot ,then call pgbeg
if (plotsave) then     !! save to file if subroutine input is on
title=trim(adjustl(tit2))//'_iteration_'//trim(adjustl(tit))
title=adjustl(title)
ier = 1
!ier=pgopen(trim(title)//'.PS/VCPS')

!if (pubplot) call pgslw(2) !! 20/10/14


else if (plotscreen) then
!ier=pgopen('/Xserve')
!call pgbgw
endif
!write(*,*) tgridplotmin,tgridplotmax,t(1)-taumax

if (ilcplot == 1) then  !! only plot the driving light curve on 1st plot if too many lc's for 1 page 17th july 15

if ((posterplot .eqv. .false.) .and. (one_page_plot .eqv. .false.)) then
!call pgsvp(bltf0,br0,upy+0.07,upy+0.08)
!if (baron) call pgtext(0.0,upy+0.08,trim(plotinfo))     !! add information bar
!call pgsch(0.65)
endif

if (ier .ne. 1) then
stop
write(*,*) 'plotting error'
endif

!!!!
if (posterplot) then
!call pgslw(6)
!!call pgsch(1.0)
endif
!call pgsch(0.7)

if (yesxray) then             !!!! Alter plot range if using xray data
xgridmax=maxval(xxray)
xgridmin=minval(xxray)

else
xgridmax=maxval(xgrid)
xgridmin=minval(xgrid)
endif
extra= extrafrac*(xgridmax-xgridmin)
yminx=xgridmin - extra
ymaxx=xgridmax + extra


!!call pgsch(0.6)
!tmin-taumax

if (one_page_plot) then
!call pgsvp(bl0,br0,upy-(upy-loy)/nplots+gap+0.02,upy+gap+0.02)
else
!call pgsvp(bl0,br0,upy-(upy-loy)/nplots+gap,upy+gap)
endif

!write(*,*) 'this swin xray'

!call pgswin(tstart,tgridplotmax,yminx,ymaxx)
if (scale_y2fit_mod) then
!call pgswin(tstart,tgridplotmax,yminx/1.1,1.1*ymaxx)
endif
write(txraymed,'(F10.2)') xraymed
txraymed=adjustl(txraymed)
!call pglabel('',trim(adjustl(dlclab)), '' )



!call pgslw(2)

!! x ray specific plotting instructions




if (yesxray) then

texty=0.8*(ymaxx-yminx)+yminx



if (sigexpand) then  !! if sigexpand turned on, plot the expanded error basrs
f=p(NPsigexpandidx)
f2 = f*f
!!call pgsci(icsigexpand)
!call pgsci(15)
if (showsigbars) then
!if (f .gt. 1) call pgerrb (6,Nxray,txray,xxray,f*erxray,2.0) !! only plot expanded error bars if expansion factor > 1
endif
!call pgsci(1)
!call pgsls(1)
tsigexpand='\gs \(2261) '
call append_n( tsigexpand, f, 3 )
tsigexpand=trim(tsigexpand)//' \gs'
tsigexpand=adjustl(tsigexpand)
textx=0.80*(t(NT)-t(1))+t(1)
!if(showsig) call pgmtxt('T',-1.3-isharedplot*1.2,0.02,0.0,trim(adjustl(tsigexpand)))
endif

!!!!! Plot the data points
yw0(1:2)=p(NPfurconstidx)
xw0(1)=tstart
xw0(2)=tmax
!call pgline(2,xw0,yw0)
!call pgpt(Nxray,txray,xxray,6)
!call pgerrb (6,Nxray,txray,xxray,erxray,1)

!! Add reduced chi square
tcisq='\gx\u2\d/'
write(titN,'(I7.3)') Nxray
titN=adjustl(titN)

tcisq=trim(tcisq)//trim(titN)//'= '

!if (sigexpand) then
!call append_n( tcisq, cisqxray/Nxray/f2, 3 )
!else
call append_n( tcisq, cisqxray/Nxray, 3 )
!endif

textx=0.03*(t(NT)-t(1))+t(1)
!!call pgtext(textx,texty,trim(tcisq))
!if (showchisq) call pgmtxt('T',-1.3,0.97,1.0,trim(tcisq))

!
endif
!call pgbox( 'bcmst', 0., 10.0, 'bcnst', 0., 10. )
if (showres) then
!call pgmtxt('T',2.3,0.5,0.5,trim(adjustl(t_xaxlab))) !HJD - 245000
else
!call pgmtxt('T',2.3,0.5,0.5,trim(adjustl(t_xaxlab)))
endif




!! Plot 0 level
if (reddotslast .eqv. .false.) then
xw0(1)=tgrid(1)
xw0(2)=tgrid(NTgrid)
yw0(1:2)=0.0
!call pgline(2,xw0,yw0)
endif


!!call pgsci(ic) !! ic = red if burnin2 not yet reached, blue otherwise


!!! Plot constant driving ligthcurve offset (furconst) (p(NPfurconstidx)
if (itsubburn .ge. 1 .and. iburnin2 .gt. 0) then  !!! Decide whether to plot current or average with envelopes
av=avparm(NPfurconstidx)!senv1furconst/itsubburn
yw0(1:2)=av
!call pgsls(1)
!!call pgsci(icnorm)
!call pgline(2,xw0,yw0)

sd= sdparm(NPfurconstidx)!sqrt( abs((senv2furconst - senv1furconst*senv1furconst/itsubburn)/(itsubburn) ) )
!call pgsls(1)
!!call pgsci(icenv)
yw0(1:2)=av + sd
!call pgline(2,xw0,yw0)
yw0(1:2)=av - sd
!call pgline(2,xw0,yw0)
!call pgsls(1)
else
!call pgsls(1)
yw0(1:2)=p(NPfurconstidx)
!!call pgsci(icnorm)
!call pgline(2,xw0,yw0)
!call pgsls(1)


if (reddotslast .eqv. .false.) then
!! plot a vertical line where the first time point starts
xw0(1:2)=tmin
yw0(1)=yminx
yw0(2)=ymaxx
!write(*,*) t(1), 'green line'
!call pgline(2,xw0,yw0)
!!call pgsci(1)
endif


endif



if (reddotslast .and. yesxray) then
!call pgsch(1.0)
!call pgslw(3)
!call pgsci(2)
!call pgpt(Nxray,txray,xxray,6)
!call pgsci(1)
!call pgslw(2)
!call pgslw(0.7)
endif



!! add error envelopes to driving lc if burnin is complete
if (iburnin2.gt. 0 .and. itsubburn .ge. 1 ) then

do it=1,Ntgrid
xgridmean(it) = senv1drive(it)/itsubburn
xgridsd= sqrt( abs((senv2drive(it) - senv1drive(it)*senv1drive(it)/itsubburn)/(itsubburn) ) )
xgridenvlodrive(it)=xgridmean(it) - xgridsd
xgridenvhidrive(it)=xgridmean(it) + xgridsd
enddo

else

do it=1,Ntgrid
xgridmean(it) = xgrid(it)
xgridsd= 0.!sqrt( abs((senv2drive(it) - senv1drive(it)*senv1drive(it)/itsubburn)/(itsubburn) ) )
xgridenvlodrive(it)=xgridmean(it) - xgridsd
xgridenvhidrive(it)=xgridmean(it) + xgridsd
enddo

endif



!!!! save driving light curve and envelopes if post burnin !!!!
call chdir('../')
open(unit = 1,file = 'xray_plot_mod.dat')
do it = 1,Ntgrid
write(1,*) tgrid(it),xgridenvlodrive(it), xgridmean(it), xgridenvhidrive(it)
enddo
close(1)
call chdir('./plots')
!!!!!

if (iburnin2.gt. 0 .and. itsubburn .ge. 1 ) then





!!call pgsci(12)
!!call pgsci(ic)
!call pgsls(1)
!!call pgsci(icenv)
!call pgline(Ntgrid,tgrid,xgridenvlodrive)
!call pgline(Ntgrid,tgrid,xgridenvhidrive)
!!call pgsci(icnorm)
!call pgsls(1)
!call pgline(Ntgrid,tgrid,xgridmean)
!!call pgsci(1)


if (savelightcurves) then
xdrivesave=xgridmean
do itt = 1,Ntgrid
sigxdrivesave(itt) = (xgridenvhidrive(itt) - xgridenvlodrive(itt))/2
enddo
endif
else
!!call pgsci(icnorm)
!call pgline(Ntgrid,tgrid,xgrid(1:Ntgrid))
if (savelightcurves) then
xdrivesave=xgrid
sigxdrivesave(1:Ntgrid) = 0
endif
endif
!!call pgsci(1)
!!




endif !! end ponly plot driving light curve 1nce if too many lc 17th july 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!Plot the delay light curves!!!!!!!!!!!!!!!!
Nlab=0 !! number of coordinate labels per plot (stops the sides getting too crowded with numbers)
if (Nlab == 0) Nlf = 0

idx_plotlc = 1
isharedplot = 0

if (ilcplot .gt. 1) nplotsreal = 0
nplots = ilc_hi - ilc_lo + 1 + nplotsreal



do i = ilc_lo,ilc_hi
!call pgsch(0.75)
if (isharelc(i) == 1) nplots = nplots -1
enddo



if (one_page_plot) then
nplots = NLC_1page+1 !new october 14 plot only 1 plot with user specified light curves plus driver
ilc_lo = 1
ilc_hi = NLC_1page
endif



nplots_ppage = ilc_hi - ilc_lo + 1
idxplot_now = 1
do i_14oct=ilc_lo,ilc_hi
if (one_page_plot) then
i = ilc_1page(i_14oct)
else
i = i_14oct
endif


!! determine how many light curves will be sharing the plot
idxcount_share = 0
do i2 = ilc_lo, ilc_hi
if (isharelc(i2) == 1) idxcount_share = idxcount_share + 1
enddo

iptmarker = 2 + isharedplot
isymbol = 0841 + isharedplot


call wavcol(wavobs(i),iclcnow)
if (cuscol_input) then
if (cuscol(i) .gt. 0) iclcnow = cuscol(i)
endif

!if (i .eq. 1) iclcnow = 1
!if (i .eq. 2) iclcnow = 6
!if (i .eq. 3) iclcnow = 5

! if plot ergs true convert the echo light curves to erg/s/cm2/A from mJy
if (plotergs) then
wavconv = wavobs(i)
wavconv2 = wavconv*wavconv
convlc = 2.9979e-8 *1e14/ wavconv2
endif

if (isharedplot == 0) then
yminlc=minval(x(lo(i):hi(i)))              !! delay lc plot limits
ymaxlc=maxval(x(lo(i):hi(i)))
else
yminlc=minval(x(lo(i-isharedplot):hi(i-isharedplot)))              !! if plotting on same plot just use same limits as previous light curve
ymaxlc=maxval(x(lo(i-isharedplot):hi(i-isharedplot)))
endif


extra=(ymaxlc-yminlc)*extrafrac
yminlc=yminlc-extra

if ((showchisq) .or. (showsig)) then
ymaxlc=ymaxlc+1.5*extra !help
else
ymaxlc = ymaxlc + 1.0*extra
endif


textlctop=0.8*(ymaxlc-yminlc)+yminlc


county=county+1

if (ilcplot .eq. 1) then
frac=(idx_plotlc+1)*dy/nplots
fracsub1=(idx_plotlc)*dy/nplots
else
frac=idx_plotlc*dy/nplots
fracsub1=(idx_plotlc-1.)*dy/nplots
endif


!!call pgsci(1)

top_old_plot = upy-fracsub1
bot_old_plot = upy-frac

if (showres) then
frac_up = 0.3
else
frac_up = 0.0
endif

diff = frac_up*(top_old_plot - bot_old_plot)  !diff is plot for residuals

!!call pgsvp(bl1,br1,bot_old_plot,bot_old_plot + diff )
!call pgsvp(bl1,br1,bot_old_plot + diff,upy-fracsub1 )

if (plotergs) then
yminlc = yminlc*convlc
ymaxlc=ymaxlc*convlc
endif

!write(*,*) 'this swin echo',i,ilc_lo,ilc_hi
!call pgswin(tmin,tmax,yminlc,ymaxlc)



!!! axis labels on middle echo plot


if (psi_echo_lab_input) then
if (isharelc(i) == 0) then
!call pgmtxt('R',5.0,0.5,0.5,trim(adjustl(echolclab(idxplot_cumulate))))
endif
else if (idxplot_now .eq. floor(nplots_ppage/2.))  then
!!if (aas_lc .eqv. .False.) then
if (plotergs) then
!call pgmtxt('R',5.0,0.5,0.5,'F\d\gl\u(\gl) / 10\u-14\d (erg s\u-1\d cm\u-2\d \A\u-1\d)')
else
!call pgmtxt('R',5.0,0.5,0.5,'F\d\gn\u(\gl) (mJy)')
endif

endif
!!!

!!!! add time axis if pubplot
if (one_page_plot) then
!if (i .eq. ilc_1page(NLC_1page)) call pgmtxt('B',2.3,0.5,0.5,trim(adjustl(t_xaxlab)))
else
if (i .eq. ilc_hi) then
if (showres) then
!!call pgmtxt('B',6.0,0.5,0.5,trim(adjustl(t_xaxlab))) !HJD -245000
else
!call pgmtxt('B',2.3,0.5,0.5,trim(adjustl(t_xaxlab)))
endif
endif
endif


!if ( (idx_plotlc .eq. NLCmax) ) then
Nysub = 0!2!0
ytick = 0.0!2.0!Nlf*(ymaxlc-yminlc)/(Nlab-1)
if (i .eq. ilc_hi) then

!!call pgsch(1.16)


if (isharedplot == 0) then
if (showres) then
!call pgbox( 'bcst', 0., 0, 'bcmstv',ytick,Nysub )
else
!call pgbox( 'bcsnt', 0., 0, 'bcmstv',ytick,Nysub )
endif
!!call pgsch(0.65)
!if (aas_lc .eqv. .false.) call pgbox( '', 0., 0, 'm', ytick,Nysub )
!!call pgsch(1.16)
!call PGMTXT ('B', 2.0, 0.5, 0.5, '')
!!call pgsch(0.65)
else
!call pgbox( '', 0., 0, '', ytick,Nysub )
endif


else

if (isharedplot == 0) then
!call pgbox( 'bcst', 0., 0, 'bcmstv', ytick,Nysub )
!if (aas_lc .eqv. .false.) call pgbox( '', 0., 10.0, 'm',ytick,Nysub )
else
!call pgbox( '', 0., 0, '',ytick, Nysub )
endif

endif


!write(*,*) 'lcplot',i,ytick,Nysub,'dsfsfd',Nlf,ymaxlc,yminlc,Nlab


!!!!! plot data points first
!!call pgsci(icpt)
!i_done_this_before = 0
!911 if (done_this_before .gt. 1) goto 999
!call pgsci(2)
do it = lo(i),hi(i),1
fnowpt = x(it)
if (plotergs) fnowpt = fnowpt*convlc
if (sigrej) then
!if (sigrej_mark_save(it) == 1) call pgpt(1,t(it),fnowpt,isymbol)
endif
enddo
!call pgsci(1)


if (sigexpand) then  !! if sigexpand turned on, plot the expanded error basrs
!!call pgsci(icsigexpand)


if (yesxray) then
f=p(NPsigexpandidx+i)
else
f=p(NPsigexpandidx+i-1)
endif
f2 = f*f

!write(*,*) sigexpand,showsigbars,f,i,'what the hell is going on!!',lo(i), hi(i),plotergs

!call pgsci(15)
if (f .gt. 0) then !! only plot the expanded error bars if > 0
!!call pgsls(3)

if (plotergs) then

if (showsigbars) then
!call pgsci(15)
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),convlc*x(lo(i):hi(i)),&
!convlc*f*err(lo(i):hi(i)),2.0)
endif
!call pgsci(1)
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),convlc*x(lo(i):hi(i)),convlc*err(lo(i):hi(i)),1.0)

else

!call pgsci(15)
if (showsigbars) then
!call pgsci(15)
do i1 = lo(i), hi(i)
!call pgerrb(6,1,t(i1),x(i1),f*err(i1),2.0)
!if(i .ge. 11) then
!write(*,*) t(i1),x(i1), err(i1), f*err(i1),'herepleaseeeee',i,f,p(NPsigexpandidx+i-1)
!read(*,*)
!endif
enddo


endif
!call pgsci(1)
!!call pgsci(icpt)
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i)),err(lo(i):hi(i)),1.0)

endif




else
!!call pgsci(ipt)

if (plotergs) then
!call pgsci(1)
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),convlc*x(lo(i):hi(i)),convlc*err(lo(i):hi(i)),1.0)
!!call pgsci(icsigexpand)
!call pgsci(15)
if (showsigbars) then
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),convlc*x(lo(i):hi(i)),convlc*f*err(lo(i):hi(i)),2.0)
endif
!call pgsci(1)
else
!call pgsci(1)

!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i)),err(lo(i):hi(i)),1.0)
!call pgsci(15)
!!call pgsci(icsigexpand)
!call pgsls(3)
if (showsigbars) then
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i)),f*err(lo(i):hi(i)),2.0)
endif
!call pgsls(1)
endif


endif
!call pgsci(1)
!!call pgsls(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! add the expansion factor to the plot
if (yesxray) then !annotate with error expansion info

if (iburnin2 .gt. 0) then
call turnp2str(p(NPsigexpandidx+i),avparm(NPsigexpandidx+i),sdparm(NPsigexpandidx+i),2.0,nam)

else
call turnp2str(p(NPsigexpandidx+i),avparm(NPsigexpandidx+i),sdparm(NPpspecidx+i),-2.0,nam)
endif

else
if (iburnin2 .gt. 0) then
call turnp2str(p(NPsigexpandidx+i-1),avparm(NPsigexpandidx+i-1),sdparm(NPsigexpandidx+i-1),2.0,nam)
else
call turnp2str(p(NPsigexpandidx+i-1),avparm(NPsigexpandidx+i-1),sdparm(NPsigexpandidx+i-1),-2.0,nam)
endif

endif



textx=0.02*(t(NT)-t(1))+t(1)


plotysch = -1.3-isharedplot*1.2 -0.3

tsigexpand='\gs \(2261) '//trim(adjustl(nam))
write(tse,'(I10)') isymbol
if (idxcount_share .gt. 0) then
tsigexpand = trim(adjustl(tsigexpand))//'  \('//trim(adjustl(tse))//')'
endif
!if (showsig) call pgmtxt('T',plotysch,0.02,0.0,trim(adjustl(tsigexpand)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


else

if (plotergs) then
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),convlc*x(lo(i):hi(i)),convlc*err(lo(i):hi(i)),1)
else
!call pgerrb (6,hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i)),err(lo(i):hi(i)),1)
endif


endif

!call pgsci (1)

!!!!!!!!!!!!!!!!!!





!!!!! ADD OFFSET LEVEL
!!call pgsci(icnorm)
flatlinex(1)=tmin
flatlinex(2)=tmax
if (version .eq. 1) then

flatliney(1:2)=p(NPgal+(i-1))+fdisk(i)				!!!!!
else if (version .eq. 3) then

if (itsubburn .ge. 1 .and. iburnin2 .gt. 0) then  !! plot offset envelopes
sd=sdparm(NPoffsetidx+i-1)!  sqrt( abs( (senv2offset(i) - senv1offset(i)*senv1offset(i)/itsubburn)/(itsubburn) ) )
offset=avparm(NPoffsetidx+i-1)!senv1offset(i)/itsubburn
flatliney(1:2)=offset

!!assign error envelopes to background component
do itg = 1,Ntgrid
bgmean = senv1(itg,i)/itsubburn
bgsd2  = abs((senv2_bg(itg,i) - senv1_bg(itg,i)*senv1_bg(itg,i)/itsubburn)/(itsubburn) )
a = offset + bgmean
b = sqrt(sd*sd + bgsd2)
offset_bg(itg)    = a
offset_bg_sd(itg) = b
enddo

!call pgsls(1)
if (plotergs) then
!call pgline(Ntgrid,tgrid,convlc*offset_bg)
else
!call pgline(Ntgrid,tgrid,offset_bg)
endif




if (plotergs) then
!call pgline(2,flatlinex,(offset_bg + offset_bg_sd)*convlc)
else
!call pgline(2,flatlinex,offset_bg + offset_bg_sd)
endif

!call pgsls(1)

else

offset = p(NPoffsetidx+i-1)
do itg = 1,Ntgrid
offset_bg(itg) = bg_save(itg,i) + offset
enddo


!call pgsls(1)

if (plotergs) then
!call pgline(Ntgrid,tgrid,convlc*offset_bg)
else
!call pgline(Ntgrid,tgrid,offset_bg)



endif

!call pgsls(1)
endif

endif
!call pgsci(1)
!!!!!!!





!if (i_done_this_before == 1) goto 999

!! Either plot updated echo lc or, if after iburnin2, error envelopes and average.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	!! Error envelopes (11th April 2014) See keith folder (avg rms routine)
!call pgsci(iclcnow)

!write(*,*) 'burnin',iburnin2,itsubburn,iteration
!read(*,*)
if (iburnin2 .gt. 0 .and. itsubburn .ge. 1 ) then  !! if past the burnin
do it=1,Ntgrid
!itpick =NTgrid-20




xgridplotcent(it) = (senv1(it,i)+senv1_bg(it,i)) /itsubburn






xgridsd_a = abs( (senv2(it,i) - senv1(it,i)*senv1(it,i)/itsubburn)/(itsubburn) )
xgridsd_b = abs( (senv2_bg(it,i) - senv1_bg(it,i)*senv1_bg(it,i)/itsubburn)/(itsubburn) )
xgridsd=sqrt( xgridsd_a + xgridsd_b )
xgridenvlo(it)=xgridplotcent(it) - xgridsd
xgridenvhi(it)=xgridplotcent(it) + xgridsd


bgopsig_save(it,i) = sqrt(xgridsd_b)
bgop_save(it,i)    = senv1_bg(it,i)/itsubburn

enddo
!call pgsls(1)
!!call pgsci(icenv)

if (plotergs) then
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridenvlo(1:interpidxmax)*convlc)
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridenvhi(1:interpidxmax)*convlc)
if (shade) call mypgshade(tgrid(1:interpidxmax),&
xgridenvlo(1:interpidxmax)*convlc,xgridenvhi(1:interpidxmax)*convlc,interpidxmax)
else
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridenvlo(1:interpidxmax))
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridenvhi(1:interpidxmax))
if (shade) call mypgshade(tgrid(1:interpidxmax),&
xgridenvlo(1:interpidxmax)*convlc,xgridenvhi(1:interpidxmax)*convlc,interpidxmax)

endif

!call pgsls(1)
else !! else just keep plotting the current echo lightcurve
xgridplotcent(1:Ntgrid) = xgridplot(1:Ntgrid,i)
xgridenvhi(1:Ntgrid) = xgridplot(1:Ntgrid,i)
xgridenvlo(1:Ntgrid) = xgridplot(1:Ntgrid,i)

bgopsig_save(1:Ntgrid,i) = 0
bgop_save(1:Ntgrid,i)    = bg_save(1:Ntgrid,i)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!call pgsci(icnorm)

if (plotergs) then
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridplotcent(1:interpidxmax)*convlc) !ome back here
else
!call pgline(interpidxmax,tgrid(1:interpidxmax),xgridplotcent(1:interpidxmax))
endif
if (savelightcurves) then
echosave(1:Ntgrid,i)=xgridplotcent(1:Ntgrid)*sd_overflow(i)+ave_overflow(i)
do itt = 1,Ntgrid
sigechosave(itt,i) = (xgridenvhi(itt) - xgridenvlo(itt))/2*sd_overflow(i)
enddo
endif
!!call pgsci(icpt)




!call pgsci(1)
!! Add reduced chi square
tcisq='\gx\u2\d/'
write(titN,'(I7.3)') hi(i)-lo(i)+1
titN=adjustl(titN)
tcisq=trim(tcisq)//trim(titN)//'  ='


call append_n( tcisq, cisqplot(i), 3 )



!read(*,*)

textx=0.45*(t(NT)-t(1))+t(1)

!if (showchisq) call pgmtxt('T',-1.3-isharedplot*1.2,0.97,1.0,trim(tcisq))
!!call pgtext(textx,textlctop,trim(tcisq))





!if (iburnin2 .gt. 0) then
!i_done_this_before = i_done_this_before + 1
!goto 911
!endif


if (reddotslast) then
!call pgsci(2)
!call pgslw(5)
!call pgsch(1.0)
if (plotergs) then
!call pgpt(hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i))*convlc,isymbol)
else
!call pgpt(hi(i)-lo(i)+1,t(lo(i):hi(i)),x(lo(i):hi(i)),isymbol)
endif
!call pgsch(0.7)
!call pgslw(2)
!call pgsci(iclcnow)
endif




!if (ilc .ge. 10) then
!Nnow = hi(i) - lo(i) + 1
!do it = 1,Nnow
!write(*,*) 'oldpt', wavem(i),i, t(lo(i)+it),x(lo(i)+it),err(lo(i)+it), plotergs
!enddo
!endif


!do ilc2 = 1,15
!do it = lo(ilc2),hi(ilc2)
!if (t(it) == 6703.0166) then
!write(*,*) t(it), x(it), err(it), p(NPsigexpandidx+10),'afafgwrew',ilc2
!!stop
!endif
!
!enddo
!enddo

if (showres) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Residual Plots !!!!!! only if post burnin



!!! inputs 1D array: x1(N1),x2(N2),y1(N1)
!!! integer N1 or N2 elements in each array
!!! output 1D array y(N2)

idxfinal =0
do i2 = 1,ilc_hi
if (isharelc(i2) .eq. 0) idxfinal = i2
enddo

if ((isharelc(i) .eq. 0) .and. (i .lt. NLC)) then! try this if plots crash
!if ((isharelc(i) .eq. 0) .and. (i .le. NLC)) then !nonly plot this if this is the first plot of a group of rttelescopes
irhi = NLC
do  i2 = i+1,NLC
if (isharelc(i2) .eq. 0) then
irhi = i2-1
exit
endif
enddo

irlo = i

itemp  = hi(irhi)
lotemp = lo(irlo)
Ntemp  = itemp-lotemp + 1

allocate(xres(Ntemp),xmod_itp(Ntemp),t_temp(Ntemp), x_temp(Ntemp),ikey_temp(Ntemp),&
ilc_temp(Ntemp),errtemp(Ntemp),ilc_tempold(Ntemp))

!irlo = min(irlo,NLC-1)
idx = irlo
!write(*,*) irlo
do i2 = 1,Ntemp
!write(*,*) irlo,idx,NLC,'what the fuck is ging on@!'
if (i2 .gt. lo(idx+1)-lo(irlo)) idx = min(idx + 1,NLC-1)
ilc_tempold(i2) = idx
enddo

call sort1d(t(lotemp:itemp),Ntemp,t_temp,ikey_temp)
do it = 1,Ntemp
x_temp(it) = x(ikey_temp(it)+lotemp-1)
errtemp(it) = err(ikey_temp(it)+lotemp-1)
ilc_temp(it) = ilc_tempold(ikey_temp(it))
enddo

call itp1d(xmod_itp,ichoplo,ichophi,tgrid,t_temp,xgridplotcent,Ntgrid,Ntemp,0)


!idx = lotemp
do it = 1,Ntemp
xres(it) = x_temp(it) - xmod_itp(it)
!idx = idx + 1
enddo

!write(*,*) 'resplot fix',i, idxfinal,i2,isharelc
!read(*,*)
!call pgsvp(bl1,br1,bot_old_plot,bot_old_plot + diff )

!xmax = max(maxval(xres),maxval(xgridplotcent))
!xmin = min(minval(xres),minval(xgridplotcent))
call avgrms( Ntemp, xres, xres_bar, xres_rms )
xmrange = 5*xres_rms
xmin = xres_bar - xmrange
xmax = xres_bar + xmrange

!call pgswin(tmin,tmax,xmin,xmax)
if (i .eq. idxfinal) then
!call pgbox('bcnst',0.0,0,'bcst',ytick,Nysub)
!call pgmtxt('B',2.5,0.5,0.5,trim(adjustl(t_xaxlab))) !HJD -245000
else
!call pgbox('bcst',0.0,0,'bcst',ytick,Nysub)
endif

!write(*,*) 'resplot',i,ytick,Nysub
!write(*,*)
inow = lotemp
ilcnow = irlo


do idx = 1,Ntemp
if (sigexpand) then
if (yesxray) then
sigexpnow = p(NPsigexpandidx + ilcnow)
else
sigexpnow = p(NPsigexpandidx + ilc_temp(idx)-1)
endif
else
sigexpnow = 1
endif
signow = errtemp(idx)


xresnow = xres(idx)
tresnow = t_temp(idx)

!if (i .eq. 10) then
!write(*,*) ilc_temp(idx), irlo, tresnow,xresnow,signow,sigexpnow
!read(*,*)
!endif

if (sigexpnow .lt. 1) then
!call pgerrb(6,1,tresnow,xresnow,signow,0.0)
!call pgsci(15)
!if (showsigbars) call pgerrb(6,1,tresnow,xresnow,signow*sigexpnow,1.0)
!call pgsci(1)
else if (sigexpnow .lt. 5) then
!call pgsci(15)
!if (showsigbars) call pgerrb(6,1,tresnow,xresnow,signow*sigexpnow,1.0)
!call pgsci(1)
!call pgerrb(6,1,tresnow,xresnow,signow,0.0)
else
!call pgsci(15)
!if (showsigbars) call pgerrb(6,1,tresnow,xresnow,signow*5,1.0)
!call pgsci(1)
!call pgerrb(6,1,tresnow,xresnow,signow,0.0)
endif
!write(*,*) tresnow, xresnow, signow,sigexpnow
!if (ilcnow .eq. 11) then
!if (tresnow == 6703.0166)then
!
!do it2 = 1,Ntemp
!ttnownow = t_temp(it2)
!xxnownow = x_temp(it2)
!write(*,*) ttnownow,xxnownow,errtemp(it2)
!enddo
!stop
!endif
!endif




if (inow .gt. hi(ilcnow)) ilcnow = ilcnow + 1
inow = inow + 1
enddo


!idx = lotemp
!do it = 1,Ntemp
!write(*,*) isharelc(i),i,it,Ntemp, t_temp(it),x_temp(it), t(idx),x(idx)
!idx = idx + 1
!enddo
!read(*,*)

!!call pgline(Ntemp,t_temp,xres)
deallocate(xres,xmod_itp,ikey_temp,x_temp,t_temp,errtemp,ilc_temp,ilc_tempold)

xw0(1) = tmin
xw0(2) = tmax
yw0(1:2) = 0
!call pgsci(iclcnow)
!call pgline(2,xw0,yw0)
!call pgsci(1)

!!call pgsvp(bl1,br1-0.2,bot_old_plot,bot_old_plot + diff )
!!call pgswin(tmin,tmax,xmin,xmax)
!!call pgbox('',0.0,0,'m',0.0,0)
endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TF PLOTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


999 ymind=0.0!minval(y(1:Ntaugrid,1))  !! Delay Plot limits
ymaxd=maxval(y(1:Ntaugrid,1))
ymaxdnlc=maxval(y(1:Ntaugrid,ilc_hi))
range=ymaxd-ymind
extra=extrafrac*range
ymaxd=ymaxd+extra
ymind=ymind
textdeltop=0.8*(ymaxd-ymind)+ymind
textdelsep=0.1*(textdeltop-ymind)



!! !! upgrade 13th August, this works out where the plots goes to zero
top = 0.d0
bot = 0.d0

if (one_page_plot) then !new october 14  make sure the plotted tf's all have same scale and is appropriate so tf is not chopped of edge of graph
ilc_hi_current = ilc_1page(NLC_1page)
else
ilc_hi_current = ilc_hi
endif
ilc_label = ilc_hi_current

do itau=1,Ntaugrid  !! upgrade 13th August, this works out where the plots goes to zero
taunow = taugrid(itau)
tfnow = y(itau,ilc_hi_current)


if (tfnow .eq. ymaxdnlc) then
idxxmaxdnlc=itau
xmaxnlc=taunow
endif
top = top + taunow*tfnow
bot = bot + tfnow
enddo
taumeannlc = top/bot

!do i3 = 1,NLC
!write(*,*)i3, y(1:10,i3)
!enddo
!write(*,*) xmaxnlc,ilc_hi_current,ilc_hi,taumeannlc,'bug here start'
!read(*,*)

if (cusnum_input) then
xmaxnlc = taumax_man(ilcplot)
taumeannlc = xmaxnlc/3
endif

!!!!!! upgrade 19 october, plots now normalise so that the yaxis  always goes from 0 to max of each tf for each lc.
tfmaxplot = 1.4*maxval(y(1:Ntaugrid,i))




!write(*,*) ilc_lo, ilc_hi
!do i2 = 1,NLC
!write(*,*) i, i2, x(lo(i2)), x(hi(i2)), t(lo(i2)), t(hi(i2))
!enddo
!if (i .eq. NLC) stop
!read(*,*)






!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! !DECIDE WHETHER TO PLOT MEAN OR CURRENT DELAY AND MEAN TRANSFER FUNCTION
!!!!!!!!!!!!!!!!

!!call pgsci(ic)    !! plot mean delay and delay function
if (iburnin2 .gt. 0) then!!!!!! decide whether to plot the mean or current delay
!!! plot the error envelopes and mean delays
top=0.d0
topb=0.d0
topt=0.d0
bot=0.d0
botb=0.d0
bott=0.d0
do itau=1,Ntaugrid
tfplotcent(itau) = senv1tf(itau,i)/itsubburn
if (savelightcurves) tfsave(itau,i)=tfplotcent(itau)

tfsd=sqrt(abs(senv2tf(itau,i) - senv1tf(itau,i)*senv1tf(itau,i)/itsubburn)/(1.*itsubburn-1.))
tfenvlo(itau)=tfplotcent(itau) - tfsd
tfenvhi(itau)=tfplotcent(itau) + tfsd
top=top+tfplotcent(itau)*taugrid(itau)
bot=bot+tfplotcent(itau)
topb=topb+tfenvlo(itau)*taugrid(itau)
botb=botb+tfenvlo(itau)
topt=topt+tfenvhi(itau)*taugrid(itau)
bott=bott+tfenvhi(itau)

if (savelightcurves) sigtfsave(itau,i) = tfsd

enddo





tfmaxplot = 1.2*maxval(tfenvhi)
!call pgsci(1)
!!!!!!!!! set up the plot window
!write(*,*) xmaxnlc,taugrid_change_scale,i,iclcnow,'bug here'
!read(*,*)
!write(*,*) cusnum_input
!stop
if ((xmaxnlc .lt. taugrid_change_scale) .or. (cusnum_input)) then  !! taugrid_change_scale at start of subroutine
!call pgsvp(bltf0,bl1-(bl1-bl0)/5,upy-frac,upy-fracsub1)

if (tfcheat) then
!call pgswin(taugrid(1),3*(taumeannlc-taugrid(1)),0.0,1.1)
else
!call pgswin(taugrid(1),3*(taumeannlc-taugrid(1)),ymind,tfmaxplot)
endif

else
!call pgsvp(bltf0,bl1,upy-frac,upy-fracsub1)  !!! update 19/10/2014 this makes the tf plots bigger

if (tfcheat) then
!call pgswin(taugrid(1),taugrid(Ntaugrid),0.0,1.1)
else
!call pgswin(taugrid(1),taugrid(Ntaugrid),ymind,tfmaxplot)
endif
!!      call pglabel( 'time delay', 'psi(tau)', '' )
endif
!!!
if (idxplot_now .eq. nplots_ppage) then
!call pgbox( 'bcnst', 0., 10, 'bcstv', 0., 0)
!call PGMTXT ('B', 2.5, 0.5, 0.5, '\gt (days)')
!!if (i .eq. ceiling(ilc_lo + 0.5*(ilc_hi-ilc_lo))) call PGMTXT ('L', 2.0, 0.5, 0.5, '\gq (\gt | \gl)')
else
if (isharedplot == 0) then
!call pgbox( 'bcst', 0., 10, 'bcstv', 0., 0)!ceiling((ymaxd-ymind)/Nlab))
else
!call pgbox( '', 0., 10, '', 0., 0)
endif
endif
!else
!!if (pubplot) call pgsch(1.16)
!hell
!if (isharedplot == 0) then

!else
!!call pgbox( '', 0., 10, '', 0., 0)
!endif

!!if (pubplot) call pgsch(0.65)
!endif
!!!!!!!!!




taucent=top/bot
taucentlo=topb/botb
taucenthi=topt/bott
taucent=top/bot !! record the position of the centroid of the time delay
xw0(1:2)=taucent
yw0(1)=0
yw0(2)=tfmaxplot







if (isharelc(i) .eq. 0) then

if (tfcheat) then

Ntaucheat = 4*Ntaugrid

allocate(tf_cheatlo(Ntaucheat),tf_cheathi(Ntaucheat),tf_cheatmid(Ntaucheat),taucheat(Ntaucheat))
dtaucheat = (taugrid(Ntaugrid) - taugrid(1))/(Ntaucheat-1)
taucheatlo = taugrid(1)
do itau = 1,Ntaucheat
taucheat(itau) = taucheatlo + (itau-1)*dtaucheat

enddo


if (pscale(NPtridx) .eq. 0) then
call cheat_tf(Ntaucheat,&
avparm(NPumbhidx), sdparm(NPumbhidx), avparm(NPcosincidx), sdparm(NPcosincidx),wavem(i),&
taucheat,tf_cheatlo,tf_cheatmid,tf_cheathi,p)
else
!write(*,*)  avparm(NPumbhidx), sdparm(NPumbhidx), avparm(NPcosincidx), sdparm(NPcosincidx),&
!avparm(NPtridx), sdparm(NPtridx), wavem(i)

call cheat_tf_alpha(Ntaucheat,&
avparm(NPumbhidx), sdparm(NPumbhidx), avparm(NPcosincidx), sdparm(NPcosincidx),&
avparm(NPtridx), sdparm(NPtridx), wavem(i),&
taucheat,tf_cheatlo,tf_cheatmid,tf_cheathi,p)

!do it = 1,Ntaucheat
!write(*,*) taucheat(it), tf_cheatlo(it), tf_cheatmid(it), tf_cheathi(it),taugrid(1),3*taumeannlc,0.0,1.1
!enddo
!stop
endif
endif


!call pgsci(iclcnow)


if (tfcheat) then
!call pgline(Ntaucheat,taucheat,tf_cheatlo)
!call pgline(Ntaucheat,taucheat,tf_cheathi)
!call pgline(Ntaucheat,taucheat,tf_cheatmid)
yw0(1) = 0.0
yw0(2) = 1.1
xw0(1:2)=taucentlo
!call pgline(2,xw0,yw0)
xw0(1:2)=taucenthi
!call pgline(2,xw0,yw0)
!call pgsls(1)


deallocate(tf_cheatlo,tf_cheathi,tf_cheatmid,taucheat)

else
!call pgsls(1)
!!call pgsci(icnorm)
!call pgline(Ntaugrid,taugrid,tfplotcent)
!!call pgsci(1)
!call pgline(Ntaugrid,taugrid,tfenvlo)
!call pgline(Ntaugrid,taugrid,tfenvhi)
!!call pgline(2,xw0,yw0)
!!call pgsls(1)
!!call pgsci(icenv)
xw0(1:2)=taucentlo
!call pgline(2,xw0,yw0)
xw0(1:2)=taucenthi
!call pgline(2,xw0,yw0)
!!call pgsci(icenv)
!call pgsls(1)

endif


if (shade) call mypgshade(taugrid,tfenvlo,tfenvhi,Ntaugrid)

!call pgsci(1)
write(ttaucent,'(F10.2)') taucent
ttaucent=adjustl(ttaucent)
write(ttaucentlo,'(F10.2)') taucentlo
ttaucentlo=adjustl(ttaucentlo)

respcent(i) = taucent
respsig(i) = (taucenthi-taucentlo)/2

!write(*,*) wavem(i), respcent(i), respsig(i)
!read(*,*)
write(ttaucenthi,'(F10.2)') abs(taucenthi-taucentlo)/2
ttaucenthi=adjustl(ttaucenthi)
tdel1='<\gt>='//trim(ttaucent)//'\(2233)'//trim(ttaucenthi) !tauhelp
tdel1=adjustl(tdel1)
!!call pgsci(1)
!!call pgsch(1.16)
!call PGMTXT ('T', -1.8, 0.97, 1.0, trim(tdel1))
!!call pgsch(0.65)


endif
!!!!!!!!!!!!!!!!!!
else
!!!!!!!!!!!!!!!!!!! if not burned in plot current tf



!!call pgsci(1)
!!!!!!!!! set up the plot window
if (xmaxnlc .lt. taugrid_change_scale) then  !! taugrid_change_scale at start of subroutine
!call pgsvp(bltf0,bl1-(bl1-bl0)/5,upy-frac,upy-fracsub1)
!call pgswin(taugrid(1),2*(taumeannlc-taugrid(1)),ymind,tfmaxplot)
!write(*,*) taugrid(1), taumeannlc, ymind, tfmaxplot, 'after tauswin'

else
!call pgsvp(bltf0,bl1,upy-frac,upy-fracsub1)  !!! update 19/10/2014 this makes the tf plots bigger
!call pgswin(taugrid(1),taugrid(Ntaugrid),ymind,tfmaxplot)
!!      call pglabel( 'time delay', 'psi(tau)', '' )
endif
!!!
if (idxplot_now .eq. floor(nplots_ppage/2.)) then
!call pgbox( 'bcnst', 0., 10, 'bcstv', 0., 0)
!!call PGMTXT ('B', 2.0, 0.5, 0.5, '\gt (days)')
!!call PGMTXT ('L', 2.0, 0.5, 0.5, '\gq (\gt | \gl)')
else
if (isharedplot == 0) then
!call pgbox( 'bcst', 0., 10, 'bcstv', 0., 0)!ceiling((ymaxd-ymind)/Nlab))
else
!call pgbox( '', 0., 10, '', 0., 0)
endif
endif


if (savelightcurves) tfsave(1:Ntaugrid,i)=y(1:Ntaugrid,i)
if (isharelc(i) .eq. 0) then
!!! add vertical lines to indicate the mean delay
top=0.d0
bot=0.d0
do itau=1,Ntaugrid
top=top+y(itau,i)*taugrid(itau)
bot=bot+y(itau,i)
enddo
taucent=top/bot !! record the position of the centroid of the time delay
xw0(1:2)=taucent
respcent(i) = taucent
respsig(i) = 0.1 !! fake out the time delay and errorbars for fvg plot if not converged
yw0(1)=0
yw0(2)=tfmaxplot
write(ttaucent,'(F10.2)') abs(taucent)
ttaucent=adjustl(ttaucent)
tdel1='<\gt> = '//trim(ttaucent)
!write(*,*) i, ttaucent, 'mean tau nlc plot',taugrid(1), taumeannlc, ymind, tfmaxplot
!write(*,*) taugrid(1), taugrid(Ntaugrid), maxval(y(1:Ntaugrid,i))
!read(*,*)
!!call pgsci(1)
!call PGMTXT ('T', -1.8, 0.97, 1.0, trim(tdel1)) !!tauhelp


!!call pgsci(icnorm)
!call pgsci(iclcnow)
!call pgline(Ntaugrid,taugrid,y(1:Ntaugrid,i))
!!call pgsci(1)
!call pgline(2,xw0,yw0)
!do itau = 1,Ntaugrid
!write(*,*) taugrid(itau), y(itau,i)
!enddo
!write(*,*) 'IM hereerere!!!', i

!read(*,*)
endif
!call pgsci(1)
!write(twav,'(I7.3)') int(wavobs(i))
!twav=adjustl(twav)



if (isharelc(i) == 0) then
!call PGMTXT ('L', 2.0, 0.5, 0.5, trim(adjustl(psilab(idxplot_cumulate))))

if (wav_ann_input) then
tdel2=trim(adjustl(wav_ann(idxplot_cumulate)))
else
tdel2=trim(adjustl(wav_ann(idxplot_cumulate)))//'\A'
endif

!call pgmtxt('T',-3.2,0.97,1.0,trim(tdel2))
endif
!!if (i .eq. 1) call pgmtxt('T', -1.5, 0.9, 1.0,'u')
!!if (i .eq. 2) call pgmtxt('T', -1.5, 0.9, 1.0,'g')
!!if (i .eq. 3) call pgmtxt('T', -1.5, 0.9, 1.0,'r')
!!if (i .eq. 4) call pgmtxt('T', -1.5, 0.9, 1.0,'i')
!!if (i .eq. 5) call pgmtxt('T', -1.5, 0.9, 1.0,'z')

!!call pgsch(0.65)

!! new 2august 15 plot light curves with equal transfer functions on same plot (this whole things is so we can expand error bars on individual points within a light curve)

endif

if ((idx_plotlc .lt. ilc_hi) .and. (i .lt. NLC)) then
if (isharelc(i+1) == 0) then
idx_plotlc = idx_plotlc + 1
isharedplot = 0
else
isharedplot = isharedplot + 1
endif
else
idx_plotlc = idx_plotlc + 1
isharedplot = 0

endif


if (isharelc(i) == 0) then
!call PGMTXT ('L', 2.0, 0.5, 0.5, trim(adjustl(psilab(idxplot_cumulate))))

if (wav_ann_input) then
tdel2=trim(adjustl(wav_ann(idxplot_cumulate)))
else
tdel2=trim(adjustl(wav_ann(idxplot_cumulate)))//'\A'
endif

!call pgmtxt('T',-3.2,0.97,1.0,trim(tdel2))
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF TIME DELAY PLOT

idxplot_now = idxplot_now + 1

if (isharelc(i) == 0) idxplot_cumulate = idxplot_cumulate + 1

!write(*,*) i,idxplot_cumulate,

!write(*,*) idxplot_now,idx_plotlc
enddo !! end the ilc loop over all lightcurves

!write(*,*) idxplot_now,idx_plotlc
!read(*,*)


ilc_lo = ilc_lo_next
ilc_hi = ilc_hi_next


!!!!! save the transfer function echo lightcurves
if (savelightcurves) then


!!call pgsci(ic)    !! plot mean delay and delay function
do i = 1,NLC
if (iburnin2 .gt. 0) then!!!!!! decide whether to plot the mean or current delay
!!! plot the error envelopes and mean delays
top=0.d0
topb=0.d0
topt=0.d0
bot=0.d0
botb=0.d0
bott=0.d0
do itau=1,Ntaugrid
tfplotcent(itau) = senv1tf(itau,i)/itsubburn
tfsave(itau,i)=tfplotcent(itau)

tfsd=sqrt(abs(senv2tf(itau,i) - senv1tf(itau,i)*senv1tf(itau,i)/itsubburn)/(1.*itsubburn-1.))
tfenvlo(itau)=tfplotcent(itau) - tfsd
tfenvhi(itau)=tfplotcent(itau) + tfsd
top=top+tfplotcent(itau)*taugrid(itau)
bot=bot+tfplotcent(itau)
topb=topb+tfenvlo(itau)*taugrid(itau)
botb=botb+tfenvlo(itau)
topt=topt+tfenvhi(itau)*taugrid(itau)
bott=bott+tfenvhi(itau)

if (savelightcurves) sigtfsave(itau,i) = tfsd

enddo

else
tfsave(1:Ntaugrid,i)=y(1:Ntaugrid,i)
endif
enddo




open(unit=1,file='modeldrive.dat')
do it=1,Ntgrid
write(1,*) tgrid(it),xdrivesave(it),sigxdrivesave(it)
enddo
close(1)

open(unit=1,file='modeltf.dat')
do itau=1,Ntaugrid
write(1,*) taugrid(itau),tfsave(itau,1:NLC)
enddo
close(1)

open(unit=1,file='modeltf_sig.dat')
do itau=1,Ntaugrid
write(1,*) taugrid(itau),sigtfsave(itau,1:NLC)
enddo
close(1)



open(unit=1,file='modellc_bg.dat')
do it=1,Ntgrid
write(1,*) tgrid(it),bgop_save(it,1:NLC)
enddo
close(1)

open(unit=1,file='modellc_bgsig.dat')
do it=1,Ntgrid
write(1,*) tgrid(it),bgopsig_save(it,1:NLC)
enddo
close(1)

!if quick_conv_itp is on then the code needs to calculate the model from parameters
!as it hasn't been evaluated at all points on time grid

!ilback = 1000!min(2,iburnin2)
!lblo = max(1,iteration - ilback)
!lbhi = max(2,iteration-1)
!nlb = lbhi - lblo + 1
!if (quick_conv_itp .and. iteration .gt. 2) then
! call cream_modop_ps(parsave(1:NP,lblo:lbhi),wavem,embhref,z,&
! Nlc,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,NPscaleidx,NPthcentidx,NPtridx,&
! nlb,Ntgrid,Ntaugrid,w,taugrid,tgrid,echosave,sigechosave)!,&
! !psiop,sigpsi)
!endif

open(unit=1,file='modellc.dat')
do it=1,Ntgrid
write(1,*) tgrid(it),echosave(it,1:NLC)
enddo
close(1)

open(unit=1,file='modellc_sig.dat')
do it=1,Ntgrid
write(1,*) tgrid(it),sigechosave(it,1:NLC)
enddo
close(1)



open(unit=1,file='cream_miscpar.dat')
write(1,*) embhref, eff
close(1)

open(unit = 1,file = 'savelc_list.txt')
do ilcs = 1,NLC
fname_tit= file_lc_in(ilcs)
fname_savelc = 'data_echo_dat_'//trim(adjustl(fname_tit))
write(1,*) trim(adjustl(fname_savelc))
open(unit = 2,file = trim(adjustl(fname_savelc)))
do its = lo(ilcs),hi(ilcs)
if (sigrej) then
write(2,*) t(its), x(its)*sd_overflow(ilcs)+ave_overflow(ilcs),&
err(its)*sd_overflow(ilcs),ervar(its)*sd_overflow(ilcs),&
sigrej_mark_save(its)
else
write(2,*) t(its), x(its)*sd_overflow(ilcs)+ave_overflow(ilcs),&
err(its)*sd_overflow(ilcs),ervar(its)*sd_overflow(ilcs), 0
endif
enddo
close(2)
enddo
close(1)






!!save merged light curves for light curves with combined wavelengths
do ilc = 1,NLC
!!start a new reference offset and stretch
!do it = 1,Ntgrid
!xnow = echosave(it,ilc)
!idlo = it
!if (xnow .ne. 0) then
!exit
!endif
!enddo

!do it = idlo,Ntgrid
!xnow = echosave(it,ilc)
!idhi = it
!if (xnow .eq. 0) then
!exit
!endif
!
!enddo
!write(*,*)'here',Ntgrid, interpidxmin,interpidxmax,Ntgrid

call avgrms(interpidxmax - interpidxmin + 1,echosave(interpidxmin:interpidxmax,ilc),&
xres_bar, xres_std)

if (isharelc(ilc) .eq. 0) then !if the first of a new wavelength of continuum light curves,
xref_bar = xres_bar
xref_std = xres_std
endif

!light curves have all been transformed to by sub mean and divide sd
!undo transformation first
ave_trans = ave_overflow(ilc)
sd_trans = sd_overflow(ilc)



!make merged data light curves
fname_tit= file_lc_in(ilc)
fname_savelc = 'merged_dat_'//trim(adjustl(fname_tit))//'.dat'

open(unit = 2,file = trim(adjustl(fname_savelc)))
do its = lo(ilc),hi(ilc)
    xt = x(its)*sd_trans + ave_trans
xnew = (xt - xres_bar)*xref_std/xres_std + xref_bar
signew = err(its)*sd_trans*xref_std/xres_std
sigmodnew = ervar(its)*sd_trans*xref_std/xres_std
if (sigrej) then
write(2,*) t(its), xnew, signew,sigmodnew, sigrej_mark_save(its)
else
    write(2,*) t(its), xnew, signew,sigmodnew, 0
!write(2,*) t(its), x(its), err(its),ervar(its), 0
endif
enddo
close(2)

!write(*,*) trim(adjustl(fname_savelc)),xc/Nc
!write(*,*) isharelc(ilc),ilc,xref_bar,xres_bar,xref_std,xres_std
!read(*,*)

!make merged model light curves
fname_tit= file_lc_in(ilc)
fname_savelc = 'merged_mod_'//trim(adjustl(fname_tit))//'.dat'
open(unit = 2,file = trim(adjustl(fname_savelc)))
do its = interpidxmin,interpidxmax
xnew = (echosave(its,ilc) - xres_bar)*xref_std/xres_std + xref_bar
signew = sigechosave(its,ilc)*xref_std/xres_std
write(2,*) tgrid(its), xnew, signew
enddo
close(2)

enddo







!prepare data to be used to get flux flux plot
!open(unit = 1, file = 'data_times.dat')

!close(1)

endif





!!call pgsch(1.0)



county=county+1
frac=county*dy/nplots
fracsub1=(county-1.)*dy/nplots


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! POWER SPECTRUM PLOT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! plot the model power spectrum based on the w0,lpha,p0 parms  !! in cycles per day


!return

minx=max(alog10(fres(1)),-4.0)
!fres(1)=10**minx
maxx=alog10(fres(NWres))
!write(*,*) fres(1), w(1)/twopi
!stop
if (mod(iteration,10) .eq. 0 .or. firstrun) then  !! update the plotting scale every 10 iterations

minyps=alog10(psres(1))-7.0
maxyps=alog10(psres(1))+3.0


endif





!! define the plotting points don't allow - inf points
!      write(*,*) 'before pspect plotting'
iw=1
do it =1,NPF,2
pspecpt(iw)=p(it)**2+p(it+1)**2
if (pspecpt(iw) .lt. 10**minyps) then
pspecpt(iw)=10**minyps
endif
if (pspecpt(iw) .gt. 10**maxyps) then
pspecpt(iw)=10**maxyps
endif
!      psres(it)=pspecpt(it)
!heelp
iw=iw+1
enddo



!xgridplotcent(it) = senv1(it,i)/itsubburn
!xgridsd=sqrt( abs( (senv2(it,i) - senv1(it,i)*senv1(it,i)/itsubburn)/(itsubburn) ) )
!xgridenvlo(it)=xgridplotcent(it) - xgridsd
!xgridenvhi(it)=xgridplotcent(it) + xgridsd

!!      call pgsch(2.0)
!      write(*,*) 'made it here pspec plot'



!call pgsci(1)
!!call pgsch(0.6)





!!!!! NEW 20/10/2014 ONly PLOT THE PS AND BOF PLOTS IF 'pubplot' is turned off.
!if (pubplot .eqv. .False.) then
if (pubplot) then
!call pgend !! end the previous plot
bl0 = 0.3
br0 = 0.9
bb0 = 0.3
br0 = 0.9

!ier=pgopen('pspec.PS/CPS')
!call pgsvp(0.2,0.9,0.2,0.9)
!call pgslw(4)
!!call pgsch(1.36)

else
!! temporary for Keith's telescope prop
!call pgsvp(bl0,0.50,upy-frac,upy-fracsub1-2*gap )
endif

!call pgswin(minx,maxx,minyps,maxyps)


!call pgbox('G',1.0,1,'G',2.0,1) !! must set pox after pgswin


!call pgmtxt('B',3.0,0.5,0.5,'(cycles / day)')
!call pgmtxt('L',3.2,0.5,0.5,'log\d10\uP(f)')


p1temp = p(NPpspecidx+2)
p2temp = p(NPpspecidx+3)
if ((iburnin2 .gt. 0) .and. (pscale(NPpspecidx+2) .gt. 0)&
.and. (pscale(NPpspecidx+3) .gt. 0) ) then !! annotate with slope info
call turnp2str(p2temp,avparm(NPpspecidx+3),sdparm(NPpspecidx+3),2.0,nam)
call turnp2str(p1temp,avparm(NPpspecidx+2),sdparm(NPpspecidx+2),2.0,nam2)
else
call turnp2str(p2temp,avparm(NPpspecidx+3),sdparm(NPpspecidx+3),-2.0,nam)
call turnp2str(p1temp,avparm(NPpspecidx+2),sdparm(NPpspecidx+2),-2.0,nam2)
endif

tit=trim(adjustl('S\d1\u='))//trim((adjustl(nam2)))
!if (showchisq) call pgmtxt('RV',-3.0,0.85,1.0,trim(adjustl(tit)))

if (p1temp .ne. p2temp) then
tit=trim(adjustl('S\d2\u='))//trim((adjustl(nam)))
!if (showchisq) call pgmtxt('RV',-3.0,0.7,1.0,trim(adjustl(tit)))
endif



!call pgbox( 'bcnsl', 0., 10, 'bcnsv', 0., 1e6)!ceiling((maxyps-minyps)/Nlab))
!      write(*,*) 'after pbbox'
!      stop


!call pgsci(1)
!call pgpt(NW,alog10(w(1: NW)/twopi),alog10(pspecpt),6)
if (w(1) .eq. 0.0) then
!call pgsci(3)
!call pgpt(1,minx,alog10(pspecpt(1)),6)

!call pgsci(1)
endif

if (1/psres(1) .eq. 0.0) psres(1) = pspecpt(1)
!! to stop the model line going dow to - inf at w=0


!call pgsci(2)



!!! decide whether to plot mean or error envelopes

if (iburnin2 .gt. 0 .and. itsubburn .ge. 1) then
psprob = 0
do iw=1,NWres
pspecsd=sqrt( abs( (senv2fur(iw) - senv1fur(iw)*senv1fur(iw)/itsubburn)/(itsubburn) ) )

rt = senv1fur(iw)/itsubburn
pspeccent(iw)= alog10(rt) + 1./rln10/(rt*rt) * pspecsd*pspecsd


!!! hell
psc = pspeccent(iw)
if ((1./psc .eq. 0) .or. (psc .ne. psc)) then
write(*,*) iw,psc, rt, itsubburn, 'pspec model problem cream_mod_new..'
psprob = 1
endif
!!!


pspeclo(iw)=pspeccent(iw) - pspecsd/rt/rln10
pspechi(iw)=pspeccent(iw) + pspecsd/rt/rln10




enddo
!call pgsls(1)
!call pgsci(icenv)

if (psprob == 0) then
!call pgline(NWres,alog10(fres),pspeclo)
!call pgline(NWres,alog10(fres),pspechi)
endif

!call pgsls(1)
!call pgsci(icnorm)
!if (psprob == 0) call pgline(NWres,alog10(fres),pspeccent)
!call pgsci(1)

else
!!!!!!!hell
psprob = 0
do iw = 1,NW
psc = psres(iw)
if ((1./psc .eq. 0) .or. (psc .ne. psc)) then
write(*,*) iw,psc, rt, itsubburn, 'pspec model problem cream_mod_new..'
psprob = 1
exit
endif
enddo
!!!!!!!

!call pgsci(icnorm)
!if (psprob  ==0) call pgline(NWres,alog10(fres),alog10(psres)) ! plot the model and ...
!call pgsci(1)
endif





yw0(1)=minyps!minval(psres)
yw0(2)=maxyps!maxval(psres)




!!!!!!!!!!!!!!!!!! Add lines to the power spectrum plot to indicate the break frequency (if present) or reference frequency
if (break) then
xw0(1:2)=alog10(p(NPpspecidx+1)*0.5/pi)
else  !! if no break
xw0(1:2)=alog10(w0mean*0.5/pi)
endif

yw0(1)=minyps
yw0(2)=maxyps
!call pgsls(3)
!call pgline(2,xw0,yw0)
!call pgsls(1)

!
!
!call pgsci(1)


!if (pubplot) call pgend
!endif !!! end if pubplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Include the step sizes of all the Fourier coefficients
!if (posterplot .eqv. .false.) then
!ip1=1
!do iw1=1,NW
!!call pgsci(3)
!ssteplog=alog10(pscale(ip1))
!csteplog=alog10(pscale(ip1+1))
!
!ypt1=max(minyps,ssteplog) !! sin step
!ypt1=min(ypt1,maxyps)
!!if ((ypt1 .eq. minyps) .or. (ypt1 .eq. maxyps)) call pgsci(1)
!
!
!ypt2=max(minyps,csteplog) !! cosine step
!ypt2=min(ypt2,maxyps)
!!if ((ypt2 .eq. minyps) .or. (ypt2 .eq. maxyps)) call pgsci(1)
!
!!call pgpt(1,alog10(0.5*w(iw1)/pi),ypt1,8)
!!call pgpt(1,alog10(0.5*w(iw1)/pi),ypt2,8)
!ip1=ip1+2
!enddo
!!call pgsci(1)
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BOF PLOT

!!!!!! include a badness of fit plot in the lower right section of the plot window


do i =1,iteration
itarray(i)=i
enddo

!!! rescale plot in order to be able to see the smaller bofs at later times
miny=minval(bofsave(1:iteration))
bofsaveplot(1:iteration) = bofsave(1:iteration) - miny
if (iteration .lt. 20) then                     !!! periodically switch between different
maxybof=bofsaveplot(1)!maxval(bofsave(1:iteration))  !! zooms on the bof plot
else
maxybof=maxval(bofsaveplot(iteration - 19:iteration))
extra=extrafracbof*(maxybof)

endif

!write(*,*) 'PLotting BOF'
if (pubplot) then


!ier=pgopen('bofplot.PS/CPS')
!call pgsvp(0.05,0.78,0.1,0.9)
!call pgslw(4)
!!call pgsch(1.36)
else
!call pgsvp(0.51,br0-0.05,upy-frac,upy-fracsub1-2*gap )
endif

!call pgswin(1.0,1.*iteration,0,maxybof+3*extra)

!call pgbox( 'bcnst', 0.0, 0, 'bcmstv',0.0,0)! (maxybof-miny)/2,0 )
!call pgmtxt('B',2.7,0.5,0.5, 'Iteration')


write(bofmin_text_a,'(I8.3)') int(miny)


bofmin_text = 'BOF - min ('//trim(adjustl(bofmin_text_a))//')'
!call pgmtxt('R',7.0,0.5,0.5,trim(adjustl(bofmin_text)))


write(tit,'(I10.5)') iteration
tit=adjustl(tit)


!call pgmtxt('RV',-2.0,0.75,1.0,'Iteration = '//trim(tit)) !help
!!if (pubplot) call pgsch(1.0)


!!!!! plot the burnin information


if (iburnin(1) .gt. 0) then
xw0(1:2) = 1.*iburnin(1)
yw0(1) = 0.0
yw0(2) = maxybof+3*extra
!call pgsci(2)
!call pgline(2,xw0,yw0) !! plot the median
!call pgsci(1)
endif


do ib = 2,Nburnin

if (iburnin(1) == 0) then
xw0(1) = 0.0
xw0(2) = real(iteration)
yw0(1:2) = bofme(1) - miny
!call pgsci(2)
!call pgline(2,xw0,yw0) !! plot the median
!call pgsci(1)
else if (iburnin(ib-1) .eq. 0) then
exit

else if (iburnin(ib-1) .gt. 0) then

if (ib .lt. Nburnin) then
if (iburnin(ib) .gt. 0) cycle
endif
!read(*,*) 'before sdev ave',ib,iburnin(ib-1:ib+1),iteration - iburnin(ib-1)

if (iteration-iburnin(ib-1)+1 .le. 1) cycle
call sdevave(bofsaveplot(iburnin(ib-1):iteration),iteration-iburnin(ib-1)+1,yavg,ysdev)
xw0(1)=real(iburnin(ib-1))
xw0(2)=real(iteration)
yw0(1:2)=bofme(ib-1) - miny
!call pgsci(4)
!call pgslw(2)
yw0(1:2) = yavg
!call pgline(2,xw0,yw0) !! current burnin
yw0(1:2) = yavg+ysdev
!call pgline(2,xw0,yw0)
yw0(1:2) = yavg-ysdev
!call pgline(2,xw0,yw0)
!call pgslw(1)
!call pgsfs(4)
!call pgrect(xw0(1),xw0(2),yavg-ysdev,yavg+ysdev)


!call pgslw(2)
yw0(1)=0.0
yw0(2)=maxybof+3*extra
xw0(1:2)=iburnin(ib-1)
!call pgline(2,xw0,yw0)
!call pgslw(1)
!call pgsci(1)

endif

enddo



!call pgsci(4)
!!call pgslw(2)
!call pgsls(2)
yw0(1)=miny
yw0(2)=maxybof+3*extra

!!!!! plot the burnin info
!do ib = 1,Nburnin
!
!iburninib = iburnin(ib)
!xw0(1:2)=1.*iburninib
!!!if (pubplot .eqv. .false.) call pgline(2,xw0,yw0)
!!call pgsci(3)
!
!!if (iburninib .gt. 0) call pgline(2,xw0,yw0)
!enddo
!call pgsci(1)
!!call pgslw(1)
!call pgsls(1)


!call pgline(iteration,itarray(1:iteration),bofsaveplot(1:iteration))

!!!!!



!!! END THE pubplot if statement!

!endif

!! close pgplot if the last plot
!call pgend

enddo !end the ilc_plot for too many plots on 1 page 17th july 15

firstrun = .false.

!write(*,*) tgrid(1),tgrid(Ntgrid),minval(t),maxval(t)
!read(*,*)










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Self containined bit of code to calculate errors and mean time delays 9/11/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (iburnin2 .gt. 0) then
do i = 1,NLC
top = 0.d0
bot = 0.d0
topb = 0.d0
topt = 0.d0
botb = 0.d0
bott = 0.d0



do itau=1,Ntaugrid
tfplotcent(itau) = senv1tf(itau,i)/itsubburn
tfsd=sqrt(abs(senv2tf(itau,i) - senv1tf(itau,i)*senv1tf(itau,i)/itsubburn)/(1.*itsubburn-1.))
tfenvlo(itau)=tfplotcent(itau) - tfsd
tfenvhi(itau)=tfplotcent(itau) + tfsd
top=top+tfplotcent(itau)*taugrid(itau)
bot=bot+tfplotcent(itau)
topb=topb+tfenvlo(itau)*taugrid(itau)
botb=botb+tfenvlo(itau)
topt=topt+tfenvhi(itau)*taugrid(itau)
bott=bott+tfenvhi(itau)
enddo
taucentlo=topb/botb
taucenthi=topt/bott
taucent=top/bot !! record the position of the centroid of the time delay
respcent(i) = taucent
respsig(i) = (taucenthi-taucentlo)/2


!write(*,*) wavem(i), respcent(i), respsig(i)
!if (i .eq. NLC) then
!stop
!endif
enddo


else
do i = 1,NLC
do itau=1,Ntaugrid
top=top+y(itau,i)*taugrid(itau)
bot=bot+y(itau,i)
enddo
taucent=top/bot !! record the position of the centroid of the time delay
respcent(i) = taucent
respsig(i) = 0.1 !! fake out the time delay and errorbars for fvg plot if not converged
enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





end subroutine

















!!!! subroutine to evaluate the mean accretion disk spectrum at a single wavelength



! subroutine to outptut the spectrum of a steady irradiated disk
!! INPUTS
!hubble const (kms-1 Mpc-1) h0in,
! efficiency 0to1 effin
! bh mass embhin (m0)
! accretion rate emdotin(m0/yr)
!zhin lampost height
!alphain
!konsint h=konstin*r**alpha
!number of radii nrin
! disk albedo 0to1 (albin)
! lower adn upper wavelengths in angstroms wloin whiin number nwin

!OUTPUT wout(nwin),fout(nwin) the wavelength and fnu arrays






























!!! THIS DOES NOT WORK OPTIMISE LATER IF TIME BUT FOR NOW,
!!! SET P0 as the average of the 1st 3 frequency amplitudes

!! subroutine to estimate starting Fourier coefficients

!!! inputs T(Ntb), x(NTB), er(NTB), p(NP), p0, w0
!! array of x,y error values
!! output fourier coefficient estimates p(NP), p0, w0 first guesses for break and
!! p0 of power spectrum


!!! subroutine 2nd aril
!! subroutine to use iterated optimal scaling to estimat ethe optimum values of the starting
! fourier terms


! inputs / outputs p(1:NPF) fourier parameters


subroutine bestguess(NTb,t,x,er,w,p,p0,w0)

real t(NTb),x(NTb),pout(NPF),xtemp(NTb),er(NTb),p(NPF),perr(NPF),w(NW),pouterr(NPF)
real med
integer,parameter:: num=1000

cospar=0.1
sinpar=0.1

xtemp(1:NTb)=x(1:NTb)  !! duplicate the data
xmed=med(xtemp,NTb)

do it=1,NTb
xtemp(it)=xtemp(it)-xmed  !!! subtract the median from the xtemp array
enddo
pout(1:NPF)=p(1:NPF)


ip=1
do iw=1,NW

do it=1,num
!write(*,*) 'bestguess iteration:',it,'sin and cos',pout(ip:ip+1),perr(ip:ip+1)

pri=w(iw)*w(iw)/(w0*w0*p0*p0)

!!! sin parameter
top=0.d0
bot=0.d0
boter=0.d0

do ita=1,NTb
sw=sin(w(iw)*t(ita))
cw=cos(w(iw)*t(ita))
ersq=er(ita)*er(ita)
top=top+(x(ita)-cospar*cw)*sw/ersq
boter=boter+sw*sw/ersq
bot=boter+sw*sw*pri
!write(*,*) 'boter',ita, boter,er(ita),ersq,sw,sw*sw/ersq
enddo

sinpar=top/bot
ersin=boter
pout(ip)=sinpar
perr(ip)=ersin
!write(*,*) it,sinpar,pout(ip),perr(ip),'dgdfgddfg'
read(*,*)

!!! cosine parameter
top=0.d0
bot=0.d0
boter=0.d0

do ita=1,NTb
sw=sin(w(iw)*t(ita))
cw=cos(w(iw)*t(ita))
ersq=er(ita)*er(ita)
top=top+(x(ita)-sinpar*sw)*cw/ersq
boter=boter+cw*cw/ersq
bot=cw*cw/ersq+cw*cw*pri
enddo

cospar=top/bot
ercos=boter
pout(ip+1)=cospar
perr(ip+1)=ercos


!!! Pspec P0 parameter
top=0.d0
bot=0.d0
boter=0.d0

ip2=1.0
do iw2=1,iw
ersinsq=pouterr(ip2)*pouterr(ip2)*2
ercossq=pouterr(ip2+1)*pouterr(ip2+1)*2
z=w0/w(iw2)
top=top+(pout(ip2)*pout(ip2)/ersinsq + pout(ip2+1)*pout(ip2+1)/ercossq)*z*z
bot=bot+(1/ersinsq + 1/ercossq)*z*z*z*z
ip2=ip2+2
enddo

p0=top/bot


enddo !! end iteration loop

do it1=1,NTb  !! update the x grid
xtemp(it1)=xtemp(it1)-pout(ip)*sin(w(iw)*t(it1)) - pout(ip+1)*cos(w(iw)*t(it1))
enddo

ip=ip+2
enddo !! end frequency loop

ip=1
do iw=1,NW
write(*,*) iw,w(iw), p(ip:ip+1),'  new  ',pout(ip:ip+1)
ip=ip+2
enddo
write(*,*) p0,w0
read(*,*)

p(1:NPF)=pout(1:NPF)  !!! update the inputted fourier terms



return
end subroutine



















































!! new version. Subtracts old parameter and adds new one when calculating fourier terms
!! 10th MAY!!! Changed to make mcmcmulticont work. Now, firstcall = .false. has been moved to after first BOF check. This means the first move in parameter space is accepted should hopefully not be a problem as after running for a long time the lowest frequency step size should be small. So any move wont bring the fit too far away from the minimum and the code hsould continue to work down.
subroutine mcmcmulti_iteration(bofold,iteration,nits,x,t,er,taugrid,psigrid,&
xgrid,tgrid,Npoints,p,psafe,pscale,pversion,pAT,pRT)

real p(NP),pscale(NP),xgrid(NTgrid),tgrid(NTgrid),t(NT),x(NT),taugrid(Ntaugrid),psigrid(Ntaugrid,NLC)
real xinterp(NT),er(NT),xxrayinterp(NXray),med,dl(NLC), pold_save(NP), pnew_save(NP),&
bofnew_save(NP), bofold_save(NP),bofsave_new(5,NP), bofsave_old(5,NP)
real,save:: cosincscale,trviscscale,umbhscale
integer pAT(NP),pRT(NP),Npoints(NLC+1), pgopen
logical pversion(NP),PSAFE(NP),firstrun,exist, altertf , CREAM_th_ex, autorej
double precision xopsum,sumtau,xgridsum
real rang,cisq,dlflatm,rms,taumean(NLC),dustAGN,dustMW
integer,save::numreject=0
integer,allocatable,save:: ip_order(:),itaumax(:), itaumin(:)
character(10):: ctit
character(1000),save:: testdir
real,allocatable:: ttemp2(:,:), xtemp2(:,:), sigtemp2(:,:), wav2(:), respcent2(:),&
respsig2(:)
! run through the loop
logical,save:: firstcall =.true.,firstpsi=.true.,firstip=.true., timeon=.false.,&
rand_order = .False., step_visc_irad_together = .false., skipfvg=.True., rejth=.false.,&
rejvar = .false., rejsig=.false., rejaffine =.false., rejtrvisc=.false.,&
rejtrirad =.false., rejcosinc =.false., rejmdot = .false.
real,allocatable,save:: fracalong(:),psinorm(:),psinormold(:)
real,save:: xraymean
character(1000):: filepath_pp,t_ff
character(10):: t0,t1,t2,t3,t4
real, allocatable:: dattemp(:,:),sigtemp(:,:), t_dattemp(:,:),taumeantemp(:),tautempsig(:)
integer, allocatable:: idxlc(:), ilclo2(:), idxlc2(:)


!quick_conv_itp = .false.
iaffine_count = 1
itimemain = 0
itimebof = 0
NWres=1000
plotsave = .false.
plotscreen=.false.



if (iplotsave .le. 0) then
plotscreen = .true.
else if (mod(iteration,iplotsave) .eq. 0) then
plotsave = .True.
end if

!do iteration=itstart,nits !! nits is set arbitrarily high, should break out after convergence
account=0   !! record the acceptance fraction
call system_clock(istarttime)
iseed=istarttime
if (iteration .eq. 1) then
firstrun = .True.
else
firstrun = .False.
endif


!! Time saving initiative 29th April 2014
!!!!!! store an array indicating the fraction between two points on the finely spaced time grid are each of the data points
!if (iteration .gt. 50) pscale(1:NPF) = 1.e-8
if (firstcall) then
varerln = 0.0
varerlnxray = 0.0
varerlnold = 0.0
varerlnxrayold = 0.0
varerlnold_affine = 0.0
varerlnxrayold_affine = 0.0


!initialize pRT_save(NP) array to 0
pRT_save(1:NP) = 0

!check for samescale file here
inquire(file='../samescale.txt',exist=CREAM_th_ex)
if (CREAM_th_ex) samescale = .true.
write(*,*) 'If there are multiiple light curves at the same wavelength cream will assume'
write(*,*) 'separate stretch, offset and sigexpand factors. If you make a file that says,'
write(*,*) '"samescale.txt" in light curve directory, cream will assign the same stretch'
write(*,*) 'offset parameter for each light curve'
write(*,*) 'samescale =',samescale
write(*,*) ''

!! do not step mmdot, tridx or inc if using only top hat tf
ilccount = 0
do ilc =1,NLC
if (wavobs(ilc) == -1.0) ilccount = ilccount + 1
enddo
if (ilccount == NLC) then
pversion(NPcosincidx) = .false.
pversion(NPumbhidx) = .false.
pversion(NPtridx) = .false.
pversion(NPtridx+1) = .false.
endif
!!

!!incorporate this for -ve delays can have v tau, need to change convolution slightly
idxtaulo = floor(taugridlo/dtgrid)
interpidxmax = Ntgrid + idxtaulo

do itau = 1,Ntaugrid
if (taugrid(itau) .ge. 0) then
idxtau0 = itau
exit
endif
enddo


!!!!!! combine different light curves but allow different error bar expansion factors
!pscale(NPscaleidx:NPscaleidx+NLC-1) = 0
!pscale(1:NPF) = 1.e-8
allocate(echosave(NTgrid,NLC),sigechosave(NTgrid,NLC))

allocate(isharelc(NLC),psimeansave(NLC),bg_now(Ntgrid))
isharelc(1:NLC) = 0
do ilc = 2,NLC
do ilc2 = 1,ilc-1
if ((wavobs(ilc) .eq. wavobs(ilc2)) .and. (wavobs(ilc) .ne. -1.0)) then
isharelc(ilc) = 1
endif
enddo
enddo

do ilc = 1,NLC

if ((isharelc(ilc) == 1) .and. (samescale)) then
pscale(NPscaleidx + ilc - 1) = 0
pversion(NPscaleidx+ilc-1) = .false.
pscale(NPoffsetidx + ilc - 1) = 0
pversion(NPoffsetidx+ilc-1) = .false.

idxlo_n = NPpolyidx + (ilc-1)*NPpoly
idxhi_n = idxlo_n + NPpoly
pscale(idxlo_n:idxhi_n) = 0
pversion(idxlo_N:idxhi_n) = .false.

endif
enddo
!!!!!!


if (iteration .eq. 1) it_BIG_atstart = 1


f0mean=w0mean/2/pi
call getcwd(dirworking)
dirworking=adjustl(dirworking)
allocate(fracalong(NT),cisqnew(NLC),psinormold(NLC),psinorm(NLC),itaumax(NLC),&
fursum(Ntgrid),fracalongxray(NXray),fres(NWres),psres(NWres),wavem(NLC),&
sdparm(NP), avparm(NP),respsig(NLC), respcent(NLC),xgridop(NTgrid),itaumin(NLC),&
psinormold_affine(NLC))

xgridop(1:Ntgrid) = 0



xgridplot(1:interpidxmin-1,1:nlc)=0.0
redshiftadd1=1.+redshift
do ilc=1,NLC  !! store the redshift wavelength
wavem(ilc)=wavobs(ilc)!/(redshiftadd1)
enddo

do iw=1,NWres
if (iteration .eq. 1 .or. firstcall) then
fres(iw)=((1.*iw-1.)/NWres*(w(NW)-w(1))+w(1))*0.5/pi
endif
enddo


ddeginc = 1.0
if (nofur == 1) ddeginc = 40.0


Nstore  = ceiling(90.0/ddeginc + 1.)

allocate(parmean(NP), psistore(Ntaugrid,Nstore))

!!!! optional apply tf scalings rather than recaluate (ONLY FOR FLAT DISK WITH T(R)^3/4)
if ((pscale(NPtridx) .ne. 0) .or. (pscale(NPtridx+1) .ne. 0) &
.or. (i_tfx == 1)) lstore = .false.
if (lstore) then
!embhref = 1.e8


umdotref = 1.0
wavref  = 4000.0

!if (nofur == 1) then
do ist = 1,Nstore
degstore = (ist - 1)*ddeginc
if (degstore .eq. 90.0) degstore = 89.0
call tfbx(taugrid(idxtau0:Ntaugrid),Ntaugrid-idxtau0+1,wavref,&
-1*umdotref,-1*umdotref,p(NPtridx),p(NPtridx+1),&
embhref,degstore,urin,psistore(idxtau0:Ntaugrid,ist),redshift)

!call tfb(taugrid(idxtau0:Ntaugrid),psistore(idxtau0:Ntaugrid,ist),Ntaugrid-idxtau0+1,&
!wavref,embhref,umdotref,uhx,eff*(1.-alb),&
!urin,degstore,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,redshift,1)
psistore(1:idxtau0-1,ist) = 0

!write(*,*) wavref, embhref, umdotref
!do it = 1,Ntaugrid
!write(*,*) taugrid(it), psistore(it,ist), degstore
!enddo
!read(*,*)

enddo
!endif



end if



!!!!



if (iteration .eq. 1) then

Nburnin = 50
allocate(fursumold(NTgrid),senv1drive(Ntgrid),senv2drive(Ntgrid)) !! fursumold must be passed into the backup file and allocated there cannot therefore allocate it twice
allocate(fursumold_affine(NTgrid),ervarold_affine(NT),erxrayvarold_affine(Nxray))
allocate(senv1offset(NLC),senv2offset(NLC),senv1fur(NWres),senv2fur(NWres),cisqplot(NLC))
allocate(sumpar(NP),sum2par(NP), covsum(NP,NP),iburnin(Nburnin),bofme(Nburnin),&
ervar(NT), erxrayvar(Nxray),ervarold(NT), erxrayvarold(Nxray))


iburnin(:) = 0
bofme(:) = 0
senv1offset=0.d0
senv2offset=0.d0
senv1furconst=0.d0
senv2furconst=0.d0
senv1(:,:)=0.d0
senv2(:,:)=0.d0
senv1_bg(:,:)=0.d0
senv2_bg(:,:)=0.d0
senv1tf(:,:)=0.d0
senv2tf(:,:)=0.d0
senv1drive(:)=0.d0
senv2drive(:)=0.d0
senv1fur(:)=0.d0
senv2fur(:)=0.d0
sumpar(:) = 0.d0
sum2par(:) = 0.d0
covsum(:,:) =0.0
sumcosinc=0.d0
summmdot=0.d0


!!assign background
bg_save(:,:) = 0
bg_now(:) = 0




do ip = NPpolyidx, NPpolyidx + NPpolytot -1,1
ipoly = mod(ip-NPpolyidx,NPpoly)
power = 1.+ipoly
poly_new = p(ip)
ilc_poly_now = floor(real(ip - NPpolyidx)/NPpoly) + 1
do it = 1,Ntgrid
bg_save(it,ilc_poly_now) = poly_new*(tgrid(it)-tref_poly(ilc_poly_now))**power
enddo
enddo


endif !endif iteration = 1






allocate(greytr(6))  !! these must be allocated each firstcall
greytr(:) = 0
greytr(2) = 1
greytr(6) = 1

parmean(:) = 0



nbackuptimesiplotsave=abs(nbackup*iplotsave)

if (yesxray) then
xraymean=avg(xxray,Nxray)
do it=1,Nxray
fracalongxray(it)=(tgrid(xrayinterpidx(it))-txray(it))/dtgrid !29/05/14
enddo
endif

do it=1,NT
fracalong(it)=(tgrid(interpidx(it))-t(it))/dtgrid
enddo


!!! include this to allow the parameters ti be stepped in random order each iteration.
allocate(ip_order(NP))
do ip = 1,NP
ip_order(ip) = ip
enddo
!!






!!! force tr visc and irad to step together new jun 13 2015 !!
if ((p(NPtridx) == p(NPtridx+1)) .and. (pscale(NPtridx) == pscale(NPtridx+1))) then
step_visc_irad_together = .true.
pscale(NPtridx+1) = pscale(NPtridx)
pversion(NPtridx+1) = .false.
endif
!!!


! parameters to tell code when to alter error bars
if (yesxray) then
ipflo = NPsigexpandidx + 1
ipfhi = ipflo + NLC - 1
ipvlo = NPvarexpandidx + 1
ipvhi = ipvlo + NLC - 1
else
ipflo = NPsigexpandidx
ipfhi = ipflo + NLC - 1
ipvlo = NPvarexpandidx
ipvhi = ipvlo + NLC - 1
endif


stepchange_affine = 1.0

endif  !! end if firstcall


!!! shuffle the order in which the parms are chosen
if (rand_order .eqv. .true.) then
call shuffle(ip_order,NP)
endif
!!!


!set on until first time we alter response function (used for simulated annealing parameter)
iat_ann = 1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ip_idx=1,NP  !!! iterate through all parameters

!write(*,*) iteration,'parameter stepping',ip,NPumbhidx, NPcosincidx, NPtridx






!!! add the random parameter cycling feature as described in http://www.bayesian-inference.com/mcmcmwg (random scan feature)
if (rand_order .eqv. .true.) then
ip = ip_order(ip_idx)
else
ip = ip_idx
endif
!!!



if (ip .eq. NPcosincidx .or. ip .eq. NPumbhidx .or. ip .eq. NPdlidx .or.&
(ip .ge. NPtridx .and. ip .lt. NPtridx+NPtr) .or. firstpsi .or. &
((ip .ge. NPthcentidx) .and. (ip .lt. NPthfwhmidx + NLC))) then
altertf = .true.
else
altertf = .false.
endif






if (timeon .and. (ip .eq. 1)) call system_clock(itimelomain) ! start of main time check

!call system_clock(icountstartfourier)



!! only step the inclination after stage 1 burnin has been reached jan 2015
if ((iburnin(1) .eq. 0) .and. (pscale(NPcosincidx) .lt. 0)) then
pversion(Npcosincidx) = .false.
else
if (pscale(NPcosincidx) .lt. 0) then
pscale(NPcosincidx) = abs(pscale(NPcosincidx))
pversion(NPcosincidx) = .True.
endif
endif
!! Turn option on and off by setting the step size in the par file to be negative





!write(*,*) P(NPscaleidx),NPscaleidx,ip
if (pversion(ip) .eqv. .False.) cycle  !!! different versions use different parameters.
!! skip these if they do not vary in this version new 28th march




!!what to do if affinenow
if (iaffinenow == 1) then
if (iaffine_count == 1) then
xgridplotsave_affine(1:Ntgrid,1:NLC) = xgridplotsave(1:Ntgrid,1:NLC)
psisave_affine(1:Ntaugrid,1:NLC) = psisave(1:Ntaugrid,1:NLC)
fdisksave_affine(1:NLC) = fdisksave(1:NLC)
p(1:NP) = parsave(1:NP,iteration-1)
ervarold_affine(1:NT) = ervarold(1:NT)
varerlnold_affine = varerlnold
if (yesxray) then
erxrayvarold_affine(1:Nxray) = erxrayvarold(1:Nxray)
varerlnxrayold_affine = varerlnxrayold
endif
psinormold_affine(:) = psinormold(:)
xgridold_affine(:)   = xgridold(:)
fursumold_affine(:)  = fursumold(:)
!set all things here that are reset after a rejection. Then accept all affine moves until the last one
! if last one rejected need to reset everything to the values here
call affine_step(N_affine,iteration-1,parsave_affine_T(1:iteration-1,1:N_affine),&
covmat_affine,p_affine,1,stepchange_affine,0)
endif


!control magnitude of stepsize in affine space
write(*,*) 'iaffine now on',AT,RT, pAT_affine, pRT_affine, iteration
if (pAT_affine .eq. AT) then
stepchange_affine = frac_stepchange_affine*stepchange_affine!
pAT_affine = 0
pRT_affine = 0
else if (pRT_affine .eq. RT) then
stepchange_affine = stepchange_affine/frac_stepchange_affine
pAT_affine = 0
pRT_affine = 0
endif


p(ip) = p_affine(iaffine_count)
!write(*,*) p_affine(1:n_affine),'affine stepping', iteration,ip,p(ip), parsave(ip,iteration-1)
!read(*,*)








!if (ip .eq. NPtridx+2 .or. ip .eq. NPtridx+3) then
!!write(*,*) 'Tv1,Ti1',p(NPtridx+2),p(NPtridx+3),pscale(NPtridx+2),pscale(NPtridx+3)
!endif
!read(*,*)





else


!decide if disc parameters are stuck.
!If so then return to a previous position in parameter space and auto accept
if ((ip .eq. NPumbhidx) .or. (ip .eq. NPcosincidx) .or. (ip .eq. Nptridx) &
.or. (ip .eq. NPtridx + 1)) then
if (pRT_save(ip) .gt. 50) then
idplo = max(1,iteration - 200)
idphi = max(iteration - 50,1)
idpast = 1 !reutrn it to the first accepted parameter!iranu(real(idplo),real(idphi),iseed)
autoacc = .true.
p(ip) = parsave(ip,idpast)
pscale(ip) = 0.1
endif
endif



!!!!! deal with the memory arrays if they equal RT then halve the step size
!!!!! if they equal AT double it
if (pAT(ip) .eq. AT) then
pscale(ip) =2*pscale(ip)
pAT(ip)=0
pRT(ip)=0
elseif (pRT(ip) .eq. RT) then
pscale(ip) =0.5*pscale(ip)
pAT(ip)=0
pRT(ip)=0
endif





if ((ip .le. NPF) .and. (pscale(ip) .eq. 0)) then !! crash if step goes to 0
write(*,*) 'Parameter', ip,' has gone to 0'
write(*,*) 'inc ip', NPcosincidx
write(*,*) 'mmdot ip', NPumbhidx
write(*,*) 'Fourier terms up to ...',NPF
stop
endif

pold=p(ip) 								 !! save old value
pold_save(ip) = pold
if (psafe(ip)) p(ip)=alog10(p(ip)) !! convert to logs if stepping in logs (if looks different because psafe is a logical array)

! step parms except inclination if iteration = 1
if (((iteration .ne. 1) .OR. (ip .ne. NPcosincidx))&
.and. (firstcall .eqv. .false.)) then
p(ip)=p(ip)+rang(0.,pscale(ip),iseed)
endif

if (psafe(ip)) p(ip)=10**p(ip)      !! convert back to de-logged



endif !end affine bit






!!now add a bit to make certain top hat parameters step together
if (ip .ge. NPthcentidx .and. ip .lt. NPthcentidx +2*NLC) then
ilc = 1
do while (ilc .lt. NLC)
thtempcentscale_now = thtempcent_scale(ilc)
do ilc2 = ilc+1,NLC
if (thtempcent_scale(ilc2) == thtempcentscale_now) then
pversion(NPthcentidx+ilc2-1) = .false.
pversion(NPthcentidx+ilc2-1+NLC) = .false.
p(NPthcentidx+ilc2-1) = p(NPthcentidx+ilc-1)
p(NPthcentidx+ilc2-1+NLC) = p(NPthcentidx+ilc-1+NLC)
else
exit
endif
enddo
ilc = ilc2
enddo
!if (ip .ge. NPthcentidx .and. ip .lt. NPthcentidx +2*NLC) then
!write(*,*) 'stepping tophats together check...'
!write(*,*) p(NPthcentidx:NPthcentidx+NLC-1)
!write(*,*) p(NPthcentidx+NLC:NPthcentidx+NLC+NLC-1)
!write(*,*) 'done check'
!endif
endif
!!


!write(*,*) pversion(NPoffsetidx:NPoffsetidx+NLC-1)
!write(*,*) pscale(NPoffsetidx:NPoffsetidx+NLC-1)
!stop


!!new 20/12/2016
!! put a ceiling on the maximum step size for the extra variance parameters. These can vary wildly and not alter the likelihood function  much
!! also put a lower limit on the minimum size of the extra variance parameters
autorej = .false.
plonow = plo(ip)
phinow = phi(ip)
if (plonow .ne. -1) then
if (p(ip) .lt. plonow) autorej = .true.
endif
if (phinow .ne. -1) then
if (p(ip) .gt. phinow) autorej = .true.
endif



if ((ip .ge. Npvarexpandidx) .and. (ip .lt. NPvarexpandidx + NPvarexpand)) then
if (pscale(ip) .gt. 1) pscale(ip) = 0.5* pscale(ip)
if (p(ip) .lt. 1.e-40) rejvar =.true.!p(ip) = 10*p(ip)
endif

if ((ip .ge. NPsigexpandidx) .and. (ip .lt. NPsigexpandidx + NPsigexpand)) then
if (p(ip) .lt. 1.e-5) rejsig =.true.
endif


if (ip .eq. NPtridx) then
if (p(ip) .lt. 0.2 .or. p(ip) .gt. 10.0) rejtrvisc = .true.
endif

if (ip .eq. NPtridx+1) then
if (p(ip) .lt. 0.2 .or. p(ip) .gt. 10.0) rejtrirad = .true.
endif

if (ip .eq. NPcosincidx) then
if ((p(NPcosincidx) .lt. 0.034) .or. (p(NPcosincidx) .gt. 1.0)) rejcosinc =.true.
endif

if (ip .eq. NPumbhidx) then
if ((p(NPumbhidx) .lt. 1.e-9) .or. (p(NPumbhidx) .gt. 1.e4) ) rejmdot =.true.
endif


if (iaffinenow == 1) then
!if ((p(NPcosincidx) .lt. 0.034) .or. (p(NPcosincidx) .gt. 1.0)) rejaffine =.true.
!if (p(NPumbhidx) .lt. 1.e-6) rejaffine =.true.
if (rejtrirad .or. rejtrvisc .or. rejvar .or. rejsig .or. rejcosinc .or. rejmdot) then
rejaffine = .true.
endif
!write(*,*) 'rejecting affine',rejaffine, p(NPcosincidx)
endif
!endif


!! check if disc parameter has gotten stuck



!if (iaffinenow == 1) then
!write(*,*)'affine parameter reset check'
!do ip2 = 1,N_affine
!write(*,*) p(ip_affine(ip2))
!enddo
!endif

!!!!!! recalculate the xgrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
furconst=p(NPfurconstidx)

if (ip .le. NPF) then

psinorcos=mod(ip,2)  !! record whether the current parameter in the fourier sum is sin or cos

if (firstcall .and. iteration .eq. 1) then

do it=1,Ntgrid
xgridsum=0.0
iw=1
do ip2=1,NPF,2
xgridsum=xgridsum+p(ip2)*swgrid(iw,it)+p(ip2+1)*cwgrid(iw,it)
iw=iw+1
enddo
fursumold(it)=xgridsum
fursum(it)=xgridsum
enddo

else

dp = p(ip) - pold
if (psinorcos .eq. 1) then !! sin
iwnow=ip/2+1
do it=1,Ntgrid
fursum(it)=fursumold(it) + dp * swgrid(iwnow,it) !! take away old and update fourier sum
enddo

else !! cosine
iwnow=ip/2
do it=1,Ntgrid
fursum(it)=fursumold(it) + dp * cwgrid(iwnow,it)
enddo

endif !! end if psinorcos .eq. 1
endif !! end if iteration .eq. 1

else if (ip .gt. NPF) then
fursum(:)=fursumold(:)


endif



if (positivity) then
xgrid(:)=furconst*exp(fursum(:)) !! take the exponential of the parameters
else
xgrid(:)=furconst+fursum(:)



!write(*,*) 'not positive',furconst,ip,NP
endif


!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!   End recalculating xgrid


!!!!!! recalculate the delay function psi ONLY if we are altering the delay parameters



!!!!!!!!!!!!! psi grid
call system_clock(istarttf)

if (ip .eq. 1) then
do ilc = 1,NLC
do itau=1,Ntaugrid
psigrid(itau,ilc)=psisave(itau,ilc) !! store the updated psigrid
enddo
enddo
endif


!! Conditions requiring alteration of transfer function(psi)
if ((altertf) .and. (rejaffine .eqv. .false.) .and. &
(rejtrirad .eqv. .false.) .and. &
(rejtrvisc .eqv. .false.) .and. &
(rejcosinc .eqv. .false.) .and. &
(rejmdot .eqv. .false.) .and. &
((iaffine_count .eq. N_affine) .or. (iaffinenow == 0))) then



uinc=acos(p(NPcosincidx))*180/pi



!!input floor on response function, mdot and trslope parameter
!if (ip .eq. NPumbhidx) then
! if (pscale(ip) .lt. 1.e-5) pscale(ip) = 10*pscale(ip)
!endif
!
!if (ip .eq. NPtridx) then
! if (pscale(ip) .lt. 1.e-5) then
! pscale(ip) = 10*pscale(ip)
! if (step_visc_irad_together) pscale(NPtridx+1) = pscale(NPtridx)
! endif
!endif
!
!if (ip .eq. NPtridx+1) then
! if (pscale(ip) .lt. 1.e-5) pscale(ip) = 10*pscale(ip)
!endif




do ilc=1,NLC

if (isharelc(ilc) == 0) then

a = p(NPcosincidx)
if ( (a .gt. 1) .or. ( a .lt. 0.034) ) exit
if (lstore) then

if (ip .eq. NPcosincidx) then
if (pscale(ip) .lt. 0.5e-3) pscale(ip) = 10*pscale(ip) ! put floor on cosinc scaling
endif


!!! condition for tophat transfer function
if (wavobs(ilc) .eq. -1.0) then!condition for top hat transfer function
taucentnow = p(NPthcentidx+ilc-1)
taufwhmnow = p(NPthcentidx+NLC+ilc-1)
if ((taucentnow - taufwhmnow .lt. taugrid(1)) .or. &
(taucentnow + taufwhmnow .gt. taugrid(Ntaugrid)) .or.&
(taufwhmnow .le. 2.5*dtgrid)) then
rejth = .true.
write(*,*) 'rejecting th parm...'
write(*,*) 'ilc now',ilc
write(*,*) 'fwhmnow',taufwhmnow
write(*,*) 'taucentnow', taufwhmnow
write(*,*) 'ip, iteration and ipindenx of 1st cent parm',ip,iteration,NPthcentidx
write(*,*)
endif
call tftophat(Ntaugrid, taugrid, psigrid(1:Ntaugrid,ilc), taucentnow, taufwhmnow)
!do itau = 1,Ntaugrid
!write(*,*) taugrid(itau), psigrid(itau,ilc)
!enddo
!write(*,*) 'making new th function',ilc,taucentnow, taufwhmnow
!read(*,*)


else
ilo = floor(uinc/ddeginc) + 1
ihi = ilo + 1
alo = (ilo-1)*ddeginc
ahi = alo + ddeginc
yscale = (uinc - alo)/ddeginc
psigrid(1:Ntaugrid,ilc) = psistore(1:Ntaugrid,ilo)&
+ yscale*(psistore(1:Ntaugrid,ihi) - psistore(1:Ntaugrid,ilo))



call tfupdate(Ntaugrid-idxtau0+1, taugrid(idxtau0:Ntaugrid),&
psigrid(idxtau0:Ntaugrid,ilc), psigrid(idxtau0:Ntaugrid,ilc), umdotref, p(NPumbhidx), 1)



call tfupdate(Ntaugrid-idxtau0+1, taugrid(idxtau0:Ntaugrid),&
psigrid(idxtau0:Ntaugrid,ilc), psigrid(idxtau0:Ntaugrid,ilc), wavref, wavem(ilc), 2)



endif
!!! endif top hat function condition


!write(*,*) ip, iteration
!read(*,*)



else


if (step_visc_irad_together) p(NPtridx+1) = p(NPtridx)

! Top hat transfer function
if (wavobs(ilc) .eq. -1.0) then!condition for top hat transfer function
taucentnow = p(NPthcentidx+ilc-1)
taufwhmnow = p(NPthcentidx+NLC+ilc-1)

!dont let the step size become less than the time lag resolution
stepfwhm   = pscale(NPthcentidx+NLC+ilc-1)
fwhmnow    = p(NPthcentidx+NLC+ilc-1)
if (stepfwhm .lt. dtgrid) then
stepfwhm = dtgrid
pscale(NPthcentidx+NLC+ilc-1) = stepfwhm
endif

if ((taucentnow - taufwhmnow .lt. taugrid(1)) .or. &
(taucentnow + taufwhmnow .gt. taugrid(Ntaugrid))) then
rejth = .true.

!come back here
endif

call tftophat(Ntaugrid, taugrid, psigrid(1:Ntaugrid,ilc), taucentnow, taufwhmnow)
else
!write(*,*) 'calling tfbx',p(NPumbhidx), embhref, p(npcosincidx), p(NPtridx)

if (i_tfx .ne. 1) then
if (wavobs(ilc) .gt. 0 .and. wavobs(ilc) .lt. 100.0) then
psigrid(idxtau0:Ntaugrid,ilc) = 0
else
call tfbx(taugrid(idxtau0:Ntaugrid),Ntaugrid-idxtau0+1,wavem(ilc),&
-1*p(NPumbhidx),-1*p(NPumbhidx),p(NPtridx),p(NPtridx+1),&
embhref,uinc,urin,psigrid(idxtau0:Ntaugrid,ilc),redshift)
endif
!     call tfb(taugrid(idxtau0:Ntaugrid),psigrid(idxtau0:Ntaugrid,ilc),&
!     Ntaugrid-idxtau0+1,wavem(ilc),embhref,p(NPumbhidx),uhx,eff*(1.-alb),&
!     urin,uinc,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,redshift,1)
else
T1v = p(NPtridx+2)
T1i = p(NPtridx+3)
slope_v = p(NPtridx)
slope_i = p(NPtridx + 1)

!psigrid(1:Ntaugrid,ilc) = 0.0
if (wavobs(ilc) .gt. 0 .and. wavobs(ilc) .lt. 100.0) then
psigrid(idxtau0:Ntaugrid,ilc) = 0
else
call tfbx(taugrid(idxtau0:Ntaugrid),Ntaugrid-idxtau0+1,wavem(ilc),&
T1v,T1i,slope_v,slope_i,embhref,uinc,urin,&
psigrid(idxtau0:Ntaugrid,ilc),redshift)
endif


!write(*,*) 'finished tfbx',p(NPumbhidx), embhref, p(npcosincidx), p(NPtridx),&
!T1v,T1i,slope_v,slope_i
!write(*,*) p(NPtridx:NPtridx+Nptr)
!write(*,*) pscale(NPtridx:NPtridx+Nptr)
!read(*,*)

!do itau = idxtau0,Ntaugrid
!write(*,*) taugrid(itau),psigrid(itau,ilc)
!enddo
!read(*,*)

endif

psigrid(1:idxtau0-1,ilc) = 0
!write(*,*) 'finished tfb',altertf,rejaffine,ip,iteration
endif
! endif tophat transfer function

endif
sum=0.d0     !! new 30th april time saving initiative
top = 0.d0
do itau=1,Ntaugrid
psigridnow = psigrid(itau,ilc)
taunow = taugrid(itau)
sum=sum+psigridnow

top = top+psigridnow*taunow
enddo

!if (wavobs(ilc) == -1.0) then
!psinorm(ilc) = 1.0
!else
psinorm(ilc)=sum
!endif

!!write(*,*) ilc, wavem(ilc),psinorm(ilc), psigrid(1:Ntau,1)
!read(*,*)
psimeansave(ilc) = top/sum !use this to evaluate whether to skip or not
!write(*,*) ilc, top, sum, taugrid(1:10),'asdadafa',psigrid(1:10,ilc)




if (firstpsi) psinormold(ilc)=psinorm(ilc)
if (version .eq. 1) then
fdisksave(1:NLC)=fdisk(1:NLC)
call disksim(p(NPdlidx),0.0,0.0,0.0,embhref,p(NPumbhidx),urin,urout,nr,uinc,diskk,&
diskka,eff,alb,0.0,uhx,wavem(ilc),fdisk(ilc),1)
endif

itaumax(ilc) = Ntaugrid
do itau = 2,Ntaugrid !! new condition 3/04/2015
!find max non 0 index of tf for omp parallel convolving later
!(omp does not like exit statements within its loops so we deal with that here)
if ((psigrid(itau,ilc) .le. 0) .and. (psigrid(itau-1,ilc) .gt.0)) then
itaumax(ilc) = itau
exit
endif
enddo

itaumin(ilc) = 1
do itau = 1, itaumax(ilc) - 2
if (psigrid(itau,ilc) .gt. 0) then
exit
endif
enddo
itaumin(ilc) = itau


else !if isharelc=1 if sharing wavelength with previous light curve
itaumin(ilc) = itaumin(ilc - 1)
itaumax(ilc) = itaumax(ilc - 1)
psigrid(1:Ntaugrid,ilc) = psigrid(1:Ntaugrid,ilc-1)
psinorm(ilc) = psinorm(ilc-1)
endif


enddo !! end ilc



!change simulated annealing parameter
if (iat_ann == 1) then
T_ann = min(1.,T_ann*frac_ann)
iat_ann = 0
endif


endif !! end if altertf ...
!!



!write(*,*) ip,iteration,'end of alter tf itaumin',itaumin
!write(*,*) ip,iteration,'end of alter tf itaumax',itaumax
firstpsi=.false.
call system_clock(iendtf)


if ((noconvolve .eqv. .false.) .and. &
((iaffine_count .eq. N_affine) .or. (iaffinenow .eq. 0))) then


do ilc=1,NLC


!write(*,*) iterastion,ip,'pspec ',rms(xgrid(1:Ntgrid),Ntgrid),&
!sqrt(p(NPpspecidx)/w(1))*p(NPpspecidx+1)
!we force the stretch factors to the light curve rms if the first iteration help
!need the area of the driving light curve (echo light curve already normalised by response function area)
!sum2 = 0.d0
!do itau = 1,Ntau
!sum2 = sum2 + psi(
!enddo
!p(NPscaleidx+ilc-1)=rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)/rms1 * stretch1




!

!! nlc_shareos if > 0 then the difference between some offset parmeters
! must remain fixed. Apply this here 05/01/2018
!write(*,*) nlc_shareos, ip, NPoffsetidx, nlc, firstcall
!read(*,*)
if (((nlc_shareos > 0) .and. (ip .ge. NPoffsetidx) .and. (ip .lt. NPoffsetidx + nlc)) &
.or. firstcall ) then
do ilcs = 1,nlc_shareos
ilc_shareos_b_now = ilc_shareos_b(ilcs)
ilc_shareos_a_now = ilc_shareos_a(ilcs)
if ((ilc_shareos_b_now == ilc) .or. (ilc_shareos_a_now == ilc)) then
p_os_a = p(NPoffsetidx + ilc_shareos_a_now - 1)
p_os_b = p(NPoffsetidx + ilc_shareos_b_now - 1)
p_os_diff = shareos_diff(ilcs)
p(NPoffsetidx + ilc_shareos_b_now - 1) = p_os_a + p_os_diff
!write(*,*) ilc, ilc_shareos_a_now, ilc_shareos_b_now,&
!p_os_a,p_os_b,p_os_diff
!read(*,*)
exit
endif
enddo
endif








if (ip .ge. NPscaleidx .and. ip .le. NPscaleidx+NLC-1) then  !! do not allow -ve stretch factors
if (p(ip) .lt. 0) p(ip)=p(ip)*(-1)
endif

idxshare = 1
do ilc2 = ilc,1, -1
if ((isharelc(ilc2) == 0) .or. (samescale .eqv. .false.)) exit
idxshare = idxshare + 1
enddo

if (version .eq. 1) then
stretch=P(NPscaleidx+(ilc-idxshare))
offset=P(NPgalidx+ilc-idxshare)+fdisk(ilc)    !! the galaxy and disk contribution are degenerate
else if (version .eq. 3) then  			!! test for fake data by setting fgal to zero
stretch=P(NPscaleidx+(ilc-idxshare))
offset=p(NPoffsetidx+ilc-idxshare)
endif



!!!!!! manual convolution
!xgridop(1:Ntgrid)=0
Ntaumax = itaumax(ilc)
Ntaumin = itaumin(ilc)

!how many to skip??? dtau part 18dec 2015
psimeannow = psimeansave(ilc)

! if (dtauskip) then
! if (psimeannow .le. 2.0) then
!  idxskip = 2 !do all lookbacks
!  else if (psimeannow .le. 3.0) then
!   idxskip = 4
!  else if (psimeannow .le. 4.0) then
!   idxskip = 6
!  else
!   idxskip = 8
! endif
!else
! idxskip = 1
!
! endif


!
! write(*,*) Ntaumin,Ntaumax,ip,ilc,itaumax
Ntau_hr = min(Ntaumax,Ntaumin + Ntau_hr_in)
if (((isharelc(ilc) == 0) .or. (quick_conv_itp)) .and. &
( (wavobs(ilc) .gt. 100.0) .or. (wavobs(ilc) .lt. 0.0) )) then       !!!!!!!! Only recalculate convolution if we are dealing with a different wavelength

if (yesxray) then
!write(*,*) interpidxmin,Ntaugrid,Ntgrid,'fsdfasdfas'


if (idxskip .gt. 0) then
!write(*,*) 'where is the fault 1'
do it=interpidxmin,interpidxmax,idxskip!Ntgrid
!new as 15/6/017 re-introduce skipping of high res response function after 10 grid spacings
!start skipping points after index Ntau_hr in response function to save time (should solve resolution issue at low delays without increasing computation time)
xopsum=0.0
do itau=Ntaumin,Ntau_hr
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*(xgrid(idx) - furconst)
enddo
if (Ntau_hr .lt. Ntaumax) then
do itau=Ntau_hr+1,Ntaumax,idxskip
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*(xgrid(idx) - furconst)*idxskip
enddo
endif
xgridop(it) = xopsum/psinorm(ilc)! * stretch + offset

enddo

else

do it=interpidxmin,interpidxmax
!new as 15/6/017 re-introduce skipping of high res response function after 10 grid spacings
!start skipping points after index Ntau_hr in response function to save time (should solve resolution issue at low delays without increasing computation time)
xopsum=0.0
do itau=Ntaumin,Ntaumax
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*(xgrid(idx) - furconst)
enddo
xgridop(it) = xopsum/psinorm(ilc)! * stretch + offset
enddo
endif



!if (ilc ==4) then

!av_gr=avg(xgrid-furconst,Ntgrid)
!av_op=avg(xgridop,Ntgrid)
!write(*,*) stretch,offset,ip,xraymean,furconst,av_gr,av_op,'ararara',xgrid(40:42),'fdfdf',xgridop(40:42)
!endif



else


if (idxskip .gt. 0) then
!write(*,*) 'where is the fault 2'
!what to do if we are skipping points
do it = interpidxmin,interpidxmax,idxskip
xopsum=0.0
do itau=Ntaumin,Ntau_hr
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*xgrid(idx)
enddo
if (Ntau_hr .lt. Ntaumax) then
do itau=Ntau_hr + 1, Ntaumax, idxskip
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*xgrid(idx)*idxskip
enddo
endif
xgridop(it) = xopsum/psinorm(ilc)! * stretch + offset
enddo

else
!else if no skipping
!write(*,*) 'where is the fault 3',interpidxmin,interpidxmax,idxskip


!do it=lo(ilc),hi(ilc)
!itpidxnow = interpidx(it)
!if (ip .eq. NPcosincidx .and. ilc .eq. 1) write(*,*) xinterp(it),it,ilc !help
!xinterp(it)=xgridop(itpidxnow-1)

!quick_conv_itp: convolve grid only on points adjacent to data 02/04/2018
if (quick_conv_itp) then

!write(*,*) interpidxmax-interpidxmin,(hi(ilc) - lo(ilc))*2
!do it = interpidxmin,interpidxmax
! xgridop(it) = offset
!enddo

do it = lo(ilc),hi(ilc)
ithi = interpidx(it)
xopsum = 0.0
do itau=Ntaumin,Ntaumax
idx=ithi-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*xgrid(idx)
enddo
xgridop(ithi) = xopsum/psinorm(ilc)



itlo = ithi - 1
xopsum = 0.0
do itau=Ntaumin,Ntaumax
idx=itlo-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*xgrid(idx)
!write(*,*) idx,itau,xopsum,xgrid(idx),psigrid(itau,ilc),'crap'
enddo
xgridop(itlo) = xopsum/psinorm(ilc)

!write(*,*) tgrid(itlo),t(it),xgridop(itlo),xgridop(ithi),'quick_conv',&
!ilc,Ntaumin,Ntaumax,psinorm(ilc)
enddo
!read(*,*)

else

do it = interpidxmin,interpidxmax
xopsum=0.0
do itau=Ntaumin,Ntaumax
idx=it-(itau-1) - idxtaulo
xopsum=xopsum+psigrid(itau,ilc)*xgrid(idx)
enddo
xgridop(it) = xopsum/psinorm(ilc)
enddo


endif
!if (ip.lt.4 .and. ilc == 1)write(*,*) 'xgridop check',xgridop(Ntgrid/2:ntgrid/2+10)

endif



endif							    !!!!!!!!


!write(*,*) 'performing interpolation between the points'
!stop
! perform interpolation on all points between

!if (idxskip .gt. 0) then !filling in the gaps if skipping points
! !write(*,*) 'where is the fault 4'
! ioplo = interpidxmin
! iophi = interpidxmin + idxskip
! xoplo = xgridop(ioplo)
! xophi = xgridop(iophi)
! nextra = mod(interpidxmax - interpidxmin + 1,idxskip)
! itpmax_temp = interpidxmax - nextra
!
! do it = interpidxmin,itpmax_temp!interpidxmax
! imod_now = mod((it - idxskip),idxskip)
! if (imod_now .ne. 0) then
!  f_now = imod_now/idxskip
!  xgridop(it) = xoplo + f_now*(xophi - xoplo)
!  else
!  xoplo = xgridop(it)
!  xophi = xgridop(min(Ntgrid,it + idxskip))
! endif
! enddo
!
! if (nextra .gt. 1) then
! ifrac = 1
! xophi = xgridop(interpidxmax)
! xoplo = xgridop(itpmax_temp)
! dxop  = xophi - xoplo
! do it = itpmax_temp + 1, itpmax_temp + nextra
! xgridop(it) = xoplo + dxop*ifrac/nextra
! ifrac = ifrac + 1
! enddo
! endif
!endif

!dont convolve if a driver
else if (((wavobs(ilc) .lt. 100.0) .and. (wavobs(ilc) .ge. 0.0))) then
xgridop(1:Ntgrid) = xgrid(1:Ntgrid)

endif !end if isharelc we do not recalculate convolution if we are just dealing separate telescopes at the same wavelength



!mod october 30th 2017
!if (firstcall .and. iteration .eq. 1) then
!do ilc2 = 1,NLC
!call avgrms(interpidxmax-interpidxmin+1,xgridop(iterpidxmin:interpidxmax),&
!xres_bar, xres_std)
!call avgrms(hi(ilc2)-lo(ilc2)+1,x(lo(ilc2):hi(ilc2)),&
!xdat_bar, xdat_std)
!stretch = xdat_std/xres_std
!p(Npscaleidx+ilc2-1) = stretch
!enddo
!stretch = p(NPscaleidx+ilc-1)
!endif


if (quick_conv_itp) then
do it = lo(ilc),hi(ilc)
ithi = interpidx(it)
!write(*,*) 'testing itp before',xgrid(ithi),xgridop(ithi),stretch,offset,x(it),ilc
xgridop(ithi) = xgridop(ithi)* stretch + offset
xgridop(ithi-1) = xgridop(ithi-1)* stretch + offset
!write(*,*) 'testing itp',xgrid(ithi),xgridop(ithi),stretch,offset,x(it),ilc,psinorm(ilc)
!write(*,*) ip,iteration
!write(*,*)
enddo
else
do it = interpidxmin,interpidxmax
xgridop(it) = xgridop(it)* stretch + offset
enddo
endif

!write(*,*) ilc, stretch, offset
!write(*,*) ip,iteration,ilc, xgridop(100:110)
!read(*,*)
if (bgvary) then
write(*,*) ip,'chchch',NPpolyidx,Nppolyidx+NPpolytot,NP
if ((ip .ge. NPpolyidx) .and. (ip .lt. NPpolyidx+Nppolytot)) then


ipoly = mod(ip-NPpolyidx,NPpoly)
power = 1.+ipoly
poly_old = pold
poly_new = p(ip)
ilc_poly_now = floor(real(ip - NPpolyidx)/NPpoly) + 1

if (ilc .eq. ilc_poly_now) then
bg_now(1:Ntgrid) = bg_save(1:Ntgrid,ilc_poly_now)
do itg = 1,Ntgrid
fchange = (poly_new - poly_old)*(tgrid(itg) - tref_poly(ilc_poly_now))**power
x_bg = bg_now(itg)  + fchange
bg_save(itg,ilc_poly_now) =  x_bg
write(*,*) x_bg, fchange, poly_new,poly_old,'cehcing shit'
read(*,*)
!if (ip .ge. NPpolyidx) write(*,*) tgrid(itg),bg_now(itg),x_bg,'polycalc',poly_old,poly_new,&
!ilc_poly_now,ip,NPpolyidx,NPpoly
enddo
!read(*,*)
endif
endif
do it = interpidxmin,interpidxmax
xgridop(it) = xgridop(it) + bg_save(it,ilc)
enddo
endif
!do it = interpidxmin,interpidxmax
!write(*,*) bg_save(it,ilc),'checking bg'
!enddo

!!!!!! End of Manual Convolution




!! the next steps are if dust is involved
!!! dust corrections will affect the determination of the luminosity distance
!! this is turned off in version 3.


if (version .eq. 1) then
do it = 1,Ntgrid
xgridop(it)=dustAGN(xgrid(it),wavem(ilc),p(NPebmvidx)) !! AGN dust (Gaskel et al 2004)
xgridop(it)=xgridop(it)/(redshiftadd1*redshiftadd1) !! redshift to observed frame
xgridop(it)=dustMW(xgridop(it),wavobs(ilc),ebmvmw) !! apply milky way dust
enddo
endif







!!!!!!!
!do it = 1,Ntgrid
! xgridplot(it,ilc)=xgridop(it)
! if ( (xgridop(it) .ne. xgridop(it)) .or. (1./xgridop(it) .eq. 0) ) then
!  write(*,*) 'problem with model',it,tgrid(it),xgridop(it),ilc
!  stop
! endif
!enddo


!!! only go from +/- 4 psisig to save computation time
!if (ip.eq. 1) write(*,*)'fractest'
do it=lo(ilc),hi(ilc)
itpidxnow = interpidx(it)
!if (ip .eq. NPcosincidx .and. ilc .eq. 1) write(*,*) xinterp(it),it,ilc !help
xinterp(it)=xgridop(itpidxnow-1) + fracalong(it)*(xgridop(itpidxnow)-xgridop(itpidxnow-1))
xinterp(it)=xinterp(it)

!if (xinterp(it) .ne. xinterp(it)) then
!write(*,*) xinterp(it),it,t(it),ilc,iteration,ip,'oh no nan problem again',&
!xgridop(itpidxnow-1),xgridop(itpidxnow),fracalong(it)
!stop
!
!endif
!if (ip .eq. NPcosincidx .and. ilc .eq. 1) write(*,*) xinterp(it),it,ilc !help
enddo


!write(*,*) ilc, xinterp(lo(ilc):hi(ilc)), xgridop(interpidx(1:10)-1)
!read(*,*) !tomtom
enddo !!! end ilc loop



end if !! end if we care about echo light curves (noconvolve) (i.e if noconvolve = false, do above)





!do itau = 1,Ntaugrid
!write(*,*) taugrid(itau),ilc, 'glargeffdsfsdfd', psigrid(itau,ilc)
!enddo
!write(*,*) ip,iteration
!read(*,*)




!! another loop for xray data !!
if (yesxray) then
do it=1,Nxray
!fracalong=(tgrid(xrayinterpidx(it))-t(it))/dtgrid  !! potential PROBLEM here when including XRAY data 29th APRIL
xxrayinterp(it) = xgrid(xrayinterpidx(it)-1) + fracalongxray(it)*(xgrid(xrayinterpidx(it))-xgrid(xrayinterpidx(it)-1))
enddo
endif
!if (ip.eq. 1) read(*,*)









if (timeon) then
call system_clock(itimehimain)
itimemain = itimemain + itimehimain - itimelomain
endif




if (timeon) call system_clock(itimelobof) !! test the time for all the components of the code comeback


!!!! calculate the varerln and varlnxray if we are on one of the error bar changing parameters
if (varexpand .or. sigexpand ) then
!write(*,*) 'ipflo',ipflo,ipfhi,ipiteration
if ((yesxray) .and. ((ip .eq. NPsigexpandidx) .or. (ip .eq. NPvarexpandidx) .or. &
(firstcall .eqv. .true.))) then !adjust the X-ray varerlnxray
sum = 0.d0
fxnow = p(NPsigexpandidx)
varexnow = p(NPvarexpandidx)
do it = 1,Nxray
aaa = erxray(it)*fxnow
erxnew2 = aaa*aaa + varexnow
erxnew = sqrt(erxnew2)
erxrayvar(it) = erxnew
sum  = sum + alog(erxnew2)
enddo
varerlnxray = sum
endif


!if (yesxray) then
! ipflo = NPsigexpandidx + 1
! ipfhi = ipflo + NLC - 1
! ipvlo = NPvarexpandidx + 1
! ipvhi = ipvlo + NLC - 1
!else
! ipflo = NPsigexpandidx
! ipfhi = ipflo + NLC - 1
! ipvlo = NPvarexpandidx
! ipvhi = ipvlo + NLC - 1
!endif
!write(*,*) 'ararararara var expand'
!p(NPvarexpandidx:NPvarexpandidx+NLC-1) = 0
if (((ip .ge. ipflo) .and. (ip .le. ipfhi)) .or. &
((ip .ge. ipvlo) .and. (ip .le. ipvhi)) .or. &
(firstcall)) then


sum = 0.d0
do ilc = 1,NLC !and for other wavelengths
if (yesxray) then
fnow   = p(NPsigexpandidx+ilc)
varnow = p(NPvarexpandidx+ilc)
else
fnow   = p(NPsigexpandidx+ilc-1)
varnow = p(NPvarexpandidx+ilc-1)
endif
ilo = lo(ilc)
ihi = hi(ilc)
do it = ilo,ihi
aaa = er(it)*fnow
ernew2 = aaa*aaa + varnow
ernew = sqrt(ernew2)
ervar(it) = ernew
sum = sum + alog(ernew2)
enddo
!write(*,*) 'varerln info',ilc,alog(ernew2),sum,fnow,varnow
enddo
varerln = sum

endif

else if (firstcall .or. sigrej) then
ervar(1:Npoints(NLC+1)) = er(1:Npoints(NLC+1))
endif
!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! calculate the badness of fit
p0temp    = p(NPpspecidx)							             !! store the pspec defaults
w0temp    = p(NPpspecidx+1)
alphatemp = p(NPpspecidx+2)
betatemp  = p(NPpspecidx+3)
!!!!!!

!do it=1,NT

!if ((ip .ge. NPscaleidx) .and. (ip .lt. NPscaleidx + NLC)) then
! write(*,*) stretch, P(NPscaleidx+(ilc-1)), ip, iteration, 'MMMMLLAAAAARGGG!'
!read(*,*)
!endif

sum=0.d0

if (noconvolve .eqv. .false.) then !! do if we are interested in echo light curves

do ilc=1,NLC

Ndat=hi(ilc)-lo(ilc)+1




if (sigrej .and. iteration .gt. itmin_sigrej) then
call cisqmod(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),ervar(lo(ilc):hi(ilc)),&
Ndat,sigrej_mode(ilc),&
sigrej_k(ilc),Npass_sigrej,sigrej_mark(lo(ilc):hi(ilc)),csq)

else
csq = cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),ervar(lo(ilc):hi(ilc)),Ndat) !! ci squared
endif
!else
!csq = cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),er(lo(ilc):hi(ilc)),Ndat)
!endif

!if (ip .ge. NPsigexpandidx .and. ip .lt. NPsigexpandidx + NLC) then
!do idx1 = lo(ilc), hi(ilc)
!write(*,*) er(idx1),ervar(idx1), p(ip)
!enddo
!write(*,*) 'testing sigexpand problem', ilc, sigrej, itmin_sigrej
!write(*,*) cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),ervar(lo(ilc):hi(ilc)),Ndat)
!write(*,*) csq
!read(*,*)
!endif


cisqnew(ilc)= csq/Ndat
!write(*,*) 'sigrejtest...',Npass_sigrej,Ndat,sigrej_k(ilc),sigrej_mode(ilc)


!if ((sigexpand .eqv. .true.) .and. (varexpand .eqv. .false.)) then  !!! modify the chi square for each light curve if we expand the error bars
!if (sigexpand) then
! if (yesxray) then
!  f=p(NPsigexpandidx+ilc)
! else
!  f=p(NPsigexpandidx+ilc-1)
! endif
! !csq=csq + alog(f) * 2*Ndat !csq/(f*f) + alog(f) * 2*Ndat
! !write(*,*) 'check this here I think you have penalised for f twice. In the above line'
! !write(*,*) 'and with varerln term'
! !read(*,*)
!endif

sum=sum+ csq

enddo


!if (ip .lt. 3 .or. ip .ge. NPthcentidx) then
!write(*,*) 'prechecking-',sum, isharelc(1:NLC)
!endif
if (varexpand .or. sigexpand) then !! add on the sum ln(var + f^2sig^^2) term to likelihood
sum = sum + varerln
endif

!if (ip .lt. 3 .or. ip .ge. NPthcentidx) then
!write(*,*) 'postvarerprechecking-',sum
!endif

end if !end the loop

bof1=sum

!write(*,*) 'bof1 after expansion',bof1

!!!!! deal with x ray data
if (yesxray) then

if (varexpand .eqv. .false.) then
cisqxray=cisq(xxray,xxrayinterp,erxray,Nxray)
else
cisqxray = cisq(xxray,xxrayinterp,erxrayvar,Nxray) + varerlnxray
endif

if (sigexpand .and. varexpand .eqv. .false.) then
fx=p(NPsigexpandidx)
cisqxray=cisqxray/(fx*fx) + alog(fx) * 2*Nxray
endif





bof1=bof1+cisqxray



endif

!!!!!!!!!!!!!!!!!!!!!!!!! END OF BOF 1



!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!          !! power spec prior (bof2) and bof(bof3)
bof2=0.0
bof2temp=0.0
bof3=0.0
iw2=1



if (break) then !! if we are including a break in the power spectrum
do ip2=1,NPF,2
a = w(iw2)/w0temp
bof2temp=p0temp*dw * (a**(-1*alphatemp))  / (1. +   a**(betatemp-alphatemp)) * 0.5
bof2=bof2+(p(ip2)*p(ip2)+p(ip2+1)*p(ip2+1)) / bof2temp
bof3=bof3 + alog(bof2temp)
iw2=iw2+1
enddo
else			!! else if no break
do ip2=1,NPF,2
bof2temp=p0temp*dw/((w(iw2)/w0temp)**alphatemp) * 0.5
bof2=bof2+(p(ip2)*p(ip2)+p(ip2+1)*p(ip2+1)) / bof2temp
bof3=bof3 + alog(bof2temp)
iw2=iw2+1
enddo
endif
bof3=bof3*2

!if (bof2 .ne. bof2) then !! error check
! write(*,*) 'its all bollocks! (bof2 ne bof2)'
! stop
!endif



!if (w0temp .eq. 0) then
!write(*,*) 'bof check bof2, and 3,' , bof2,bof3,p(NPpspecidx:NPpspecidx+2)
!read(*,*)
!endif
!!!!!!                   !! alpha prior





bofnew=bof1+bof2+bof3



! Add priors on power spectrum parameters if present
!!!!!!                  !! P0
if (sigp0square .gt. 0.0) then
bof4temp=(p0temp-p0mean)
bof4=bof4temp*bof4temp/sigp0square
bofnew=bofnew+bof4
endif

!!!!!!                  !! w0
if (sigw0square .gt. 0.0) then
bof5temp=(w0temp-w0mean)
bof5=bof5temp*bof5temp/sigw0square
bofnew=bofnew+bof5
endif


!!!!!!					!! alpha
if (sigalphasquare .gt. 0) then
bof6temp=(alphatemp-alphamean)
bof6=bof6temp*bof6temp/sigalphasquare
bofnew = bofnew + bof6
endif

!!!!!!					!! beta
if (sigbetasquare .gt. 0) then
bof7temp=(betatemp-betamean)
bof7=bof7temp*bof7temp/sigbetasquare
bofnew = bofnew + bof7
endif


!!!!!!					!! cos inclination prior
if (cosinc_cut .gt. 0) then
cosinc_now = p(NPcosincidx)
bof8temp = (cosinc_now - cosinc_cut)!-2*alog( 1. / ( 1. + (cosinc_cut/cosinc_now)**alpha_cosinc ) )
bof8     = bof8temp*bof8temp / sigcosinc2
bofnew = bofnew + bof8
endif




!!!!!!                  !! mmdot prior
bof9 = 0
if (sigmmdotsquare .gt. 0) then
tempmmdot = p(NPumbhidx)
bof9temp = (tempmmdot - pri_mmdotmean)
bof9     = bof9temp*bof9temp / sigmmdotsquare
bofnew   = bofnew + bof9
endif
!write(*,*) sigmmdotsquare, tempmmdot, pri_mmdotmean,bof9,iteration, ip




!!!!                    !! top hat BOF prior terms
bof10 = 0
bof11 = 0
do ilc = 1,NLC
if (wavobs(ilc) .ne. -1.0) then
cycle
else

!top hat cent first
pricent_width = thtempcent_pri_width(ilc)
if (pricent_width .lt. 0) then
cycle
else
tempcent      = p(NPthcentidx+ilc-1)
pricent_cent  = thtempcent_pri(ilc)
pricent_width = thtempcent_pri_width(ilc)
aaa = (tempcent - pricent_cent)
bof10temp = aaa*aaa/(pricent_width*pricent_width)
bof10 = bof10 + bof10temp
endif

!top hat fwhm next
prifwhm_width = thtempfwhm_pri_width(ilc)
if (prifwhm_width .lt. 0) then
cycle
else
tempfwhm      = p(NPthfwhmidx+ilc-1)
prifwhm_cent  = thtempfwhm_pri(ilc)
prifwhm_width = thtempfwhm_pri_width(ilc)
aaa = (tempfwhm - prifwhm_cent)
bof11temp = aaa*aaa/(prifwhm_width*prifwhm_width)
bof11 = bof11 + bof11temp
endif
endif
enddo

bofnew = bofnew + bof10 + bof11


if (sigexpand) then
bof12 = 0.d0
if (yesxray) then
pfwhmnow = sigexpandpri_fwhm(1)
psenow = p(NPsigexpandidx)

if (pfwhmnow .lt. 0) then
bof12 = 0
else
aaa    = alog10(psenow/sigexpandpri_cent(1))
bof12  = aaa*aaa/(pfwhmnow*pfwhmnow)
endif

do ilc = 1,NLC
pfwhmnow = sigexpandpri_fwhm(1+ilc)
if (pfwhmnow .lt. 0) cycle
psenow = p(NPsigexpandidx+ilc)
aaa = alog10(psenow/sigexpandpri_cent(1+ilc))
bof12 = bof12 + aaa*aaa/(pfwhmnow*pfwhmnow)
enddo
else
do ilc = 1,NLC
pfwhmnow = sigexpandpri_fwhm(ilc)
if (pfwhmnow .lt. 0) cycle
psenow = p(NPsigexpandidx+ilc-1)
aaa = alog10(psenow/sigexpandpri_cent(ilc))
bof12 = bof12 + aaa*aaa/(pfwhmnow*pfwhmnow)
enddo
endif
bofnew = bofnew + bof12
endif





!!!!!!				 !! prior on the background parameters
bof13 = 0
if (bgvary) then
iplo = NPpolyidx
iphi = iplo      + NPpolytot - 1
do ip_pol = iplo,iphi
ipoly = mod(ip_pol-NPpolyidx,NPpoly) + 1
ipn = floor(real(ip_pol - NPpolyidx)/NPpoly) + 1

ppsd = polyin_pri(ipoly,ipn)
if (ppsd .gt. 0) then
bof13temp = ((p(ip_pol) - polyin_pri_mean(ipoly,ipn))/ppsd )**2
bof13 = bof13 + bof13temp
endif
enddo

bofnew = bofnew + bof13
endif



!new 18/12/17 correct bug where priors counted multiple times in pricream
!new 18/12/17 if pricream_sd(i) -ve then do not include prior
!! all other priors
if (pricream) then
pricream_bofsumnew = 0.d0



do ipc = 1,npricream
pricream_idxnow = pricream_idx(ipc)
valcream_now = p(pricream_idxnow)

pricream_sdnow  = pricream_sd(ipc)
if (pricream_sdnow .lt. 0) then
pricream_bofnow = 0

else if (firstcall .or. pricream_idxnow == ip) then
pricream_meannow = pricream_mean(ipc)
top = pricream_meannow - valcream_now
pricream_bofnow  = top*top/pricream_sdnow/pricream_sdnow
else
pricream_bofnow = pricream_bofold(ipc)
endif


pricream_bofsumnew = pricream_bofsumnew + pricream_bofnow
pricream_bofnew(ipc) = pricream_bofnow

!if ((ip .ge. NPsigexpandidx) .and. (ip .lt. NPsigexpandidx+NPsigexpand)) then
!write(*,*) 'testing pricream',npricream,ipc

!write(*,*) ip,p(ip)
!write(*,*) pricream_idx(ipc)
!write(*,*) pricream_mean(ipc)
!write(*,*) pricream_sd(ipc)
!write(*,*) pricream_bofnow,pricream_bofsumnew,'pricream_bofsum   original bof',bofnew

!endif


enddo
bofnew = bofnew + pricream_bofsumnew
!write(*,*) pricream_bofnow,pricream_bofsumnew,'pricream_bofsum   original bof',bofnew
!if ((ip .ge. NPsigexpandidx) .and. (ip .lt. NPsigexpandidx+NPsigexpand)) then
!write(*,*) pricream_idx
!write(*,*) pricream_mean
!write(*,*) pricream_sd
!write(*,*) pricream_bofnew
!write(*,*) pricream_bofsumnew,'pricream_bofsum   original bof',bofnew - pricream_bofsumnew
!!read(*,*)
!endif


endif




bofreject(IP)=bofnew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ((bofnew .ne. bofnew) .and. (rejvar .eqv. .false.) .and. (rejth .eqv. .false.)&
!.and. (rejsig .eqv. .false.)) then
!open(unit=12,file='bofnan_info.dat',access='append')
!write(12,*) 'parameter indices: mmdot,cos,pspec,sigexpand,varexpand,thcent',rejvar,rejth,rejsig
!write(12,*) NPumbhidx,NPcosincidx,NPpspecidx,NPsigexpandidx,Npvarexpandidx,NPthcentidx
!write(12,*) 'parameter',ip,p(ip),pold, pscale(ip)
!write(12,*) 'bof nan', bof1,bof1sigexpand,bof2,bof3,bof4, cisqxray
!write(12,*) 'bof2 stuff', p0temp,dw,w0mean,alphatemp
!write(12,*) 'fancy expansion stuff',bof1/(fx*fx),alog(f) * 2*(Nxray+NT),fx
!write(12,*) 'priorbofs..p0,w0,alpha,beta,cosinc,mmdot',bof4,bof5,bof6,bof7,bof8,bof9
!write(12,*) ''
!close(12)
!endif



pnew_save(ip) = p(ip)
bofnew_save(ip) = bofnew
bofold_save(ip) = bofold




if ((timeon)) then
call system_clock(itimehibof)
itimebof = itimebof + itimehibof - itimelobof
endif




!! save
!if (ip .eq. NPcosincidx) then
!open(unit=1,file='embhbofs_pre.dat',access='append')
!write(1,*) iteration, bofnew, bofold, p(ip), pold, pscale(ip),bof1, bof2, bof3, bof4, bof5, &
!bof6, bof7, bo8temp, bof9, p(NPscaleidx), psinorm(1)
!close(1)
!open(unit=1,file='embhlcs_pre.dat',access='append')
!write(1,*) t(lo(1):hi(1))
!write(1,*) x(lo(1):hi(1))
!write(1,*) xinterp(lo(1):hi(1))
!write(1,*) t(lo(2):hi(2))
!write(1,*) x(lo(2):hi(2))
!write(1,*) xinterp(lo(2):hi(2))
!write(1,*) t(lo(3):hi(3))
!write(1,*) x(lo(3):hi(3))
!write(1,*) xinterp(lo(3):hi(3))
!close(1)
!endif
!if (ip .lt. 3 .or. ip .ge. NPthcentidx) then
!
! sum = 0.d0
! sum1 = 0.d0
! do ilc = 1,NLC
! !write(*,*) cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),ervar(lo(ilc):hi(ilc)),Ndat),&
! !cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),er(lo(ilc):hi(ilc)),Ndat)
! Ndat=hi(ilc)-lo(ilc)+1
! sum = sum + cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),ervar(lo(ilc):hi(ilc)),Ndat)
! sum1 = sum1 + cisq(x(lo(ilc):hi(ilc)),xinterp(lo(ilc):hi(ilc)),er(lo(ilc):hi(ilc)),Ndat)
! enddo
! sum = sum
! sum1 = sum1
! write(*,*) ip, sum,bof1,bofnew,'checking',bofold!yesxray,varexpand,sigexpand
! !iteration,ip,NPvarexpandidx,bof1,bofold
!
!! ier = pgopen('/Xserve')
!! call pgsvp(0.3,0.9,0.55,0.9)
!! call pgswin(t(1),t(NT),0.5,3.5)
!! call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!! call pgerrb(6,hi(1)-lo(1)+1,t(lo(1):hi(1)),x(lo(1):hi(1)),er(lo(1):hi(1)),1.0)
!! call pgline(hi(1)-lo(1)+1,t(lo(1):hi(1)),xinterp(lo(1):hi(1)))
!! call pgline(NTgrid,tgrid,xgridplot(1:Ntgrid,1))
!
!! call pgsvp(0.05,0.25,0.55,0.9)
!! call pgswin(-2.,2.,0.0,1.0)
!! call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!! call pgline(Ntaugrid,taugrid,psigrid(1:Ntaugrid,1))
!
!! call pgsvp(0.3,0.9,0.1,0.5)
!! call pgswin(t(1),t(NT),24.0,28.0)
!! call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!! call pgerrb(6,hi(2)-lo(2)+1,t(lo(2):hi(2)),x(lo(2):hi(2)),er(lo(2):hi(2)),1.0)
!! call pgline(hi(2)-lo(2)+1,t(lo(2):hi(2)),xinterp(lo(2):hi(2)))
!! call pgline(NTgrid,tgrid,xgridplot(1:Ntgrid,2))
!
!! call pgsvp(0.05,0.25,0.1,0.5)
!! call pgswin(-2.,2.,0.0,1.0)
!! call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!! call pgline(Ntaugrid,taugrid,psigrid(1:Ntaugrid,2))
!
! write(*,*)
! read(*,*)
!! call pgend
!endif

!if (ip .eq. NPthcentidx+4) then
!open(unit = 1, file='thbofs.dat',access='append')
!write(1,*) iteration, bofnew, bofold, p(ip), pold, pscale(ip)
!close(1)
!endif



!if (ip .eq. NPumbhidx) then
!open(unit=1,file='embhbofs.dat',access='append')
!write(1,*) iteration, bofnew, bofold, p(ip), pold, pscale(ip),bof1, bof2, bof3, bof4, bof5, &
!bof6, bof7, bo8temp, bof9, p(NPscaleidx), psinorm(1)
!close(1)
!open(unit=1,file='embhlcs.dat',access='append')
!write(1,*) t(lo(1):hi(1))
!write(1,*) x(lo(1):hi(1))
!write(1,*) xinterp(lo(1):hi(1))
!write(1,*) t(lo(2):hi(2))
!write(1,*) x(lo(2):hi(2))
!write(1,*) xinterp(lo(2):hi(2))
!write(1,*) t(lo(3):hi(3))
!write(1,*) x(lo(3):hi(3))
!write(1,*) xinterp(lo(3):hi(3))
!
!!write(1,*) er(lo(1):hi(1))
!close(1)
!open(unit=1,file='embhlcs_mod.dat',access='append')
!write(1,*) tgrid(1:Ntgrid)
!write(1,*) xgridplot(1:Ntgrid,1)
!write(1,*) xgridplotsave(1:Ntgrid,1)
!close(1)
!open(unit=1,file='embhtf.dat',access='append')
!write(1,*) psigrid(1:Ntaugrid,1)
!close(1)
!open(unit=1,file='ihatescience.dat',access='append')
!if (yesxray) then
!write(1,*) stretch, offset, 1, idxshare, xgrid(1:Ntgrid)
!else
!write(1,*) stretch, offset, 0, idxshare, xgrid(1:Ntgrid)
!endif
!
!close(1)
!endif

!!!!!!!!!!!!!!! Compare the Badness of fits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ACCEPT!



!if (ip .ge. NPvarexpandidx) then
!write(*,*)
!write(*,*) 'test new var parms'
!write(*,*) ip,bof1, bofnew,bofold,p(ip), pold
!write(*,*) varerlnxray, varerlnxrayold
!write(*,*) varerln, varerlnold
!write(*,*) er(1:10)
!write(*,*) ervar(1:10)
!write(*,*) p(NPsigexpandidx),p(NPsigexpandidx+1)
!!read(*,*)
!endif
!
!if ((ip .eq. NPthcentidx) .or. (ip .eq. NPthfwhmidx)) then
!write(*,*)
!write(*,*) 'test new th parms'
!write(*,*) ip,bof1,bofnew,bofold,p(ip), pold
!!read(*,*)
!endif
!
!if (ip .eq. NPumbhidx) then
!write(*,*)
!write(*,*) 'test emmdot parameter'
!write(*,*) ip,bof1, bofnew, bofold, p(ip), pold
!!read(*,*)
!endif

!do ilc = 1,NLC
!if ((ip .eq. p(NPthcentidx+NLC+ilc-1))) then
!write(*,*) ilc,ip,NPthcentidx,iteration
!write(*,*) "FWHM", p(NPthcentidx+NLC+ilc-1),'old=',pold
!write(*,*) "Cent", p(NPthcentidx+ilc-1)
!write(*,*) 'bofnew, bofold', bofnew, bofold
!
!
!endif
!enddo

!if (iteration .le. 110) then!if (iaffinenow .eq. 0 .and. iteration .le. 104) then
!write(*,*) ip, bofnew, bofold
!endif


!if ((ip .eq. NPthcentidx) .and. (ip .lt. NPthcentidx + NLC)) then
!write(*,*) p(ip), pold, 'th parameter',ilc
!read(*,*)
!endif



!if (ip .ge. NPumbhidx - 10 .and. ip .le. NPcosincidx) then
!write(*,*) xgridop(Ntgrid/2:Ntgrid/2+10)
!write(*,*)xinterp(1:10)
!write(*,*) x(1:10)
!write(*,*) ervar(1:10)
!write(*,*) ip,p(ip), pold,'mbh and cos parms', bofnew,bofold,&
!NPumbhidx, NPcosincidx
!write(*,*) 'cisq',csq
!write(*,*)
!read(*,*)
!endif
if (iaffine_count .eq. N_affine .and. iaffinenow .eq. 1) then

write(*,*) 'stepping affine parameters.... checking--> stepchange_affine',stepchange_affine,&
bofnew1,bof1,'...',bofnew2,bof2,'...',bofnew3,bof3,'...'

do id =1,n_affine
write(*,*) ip_affine(id),p(ip_affine(id)),parsave(ip_affine(id),iteration-1),bofnew,bofold
enddo
endif


!if (iaffinenow .eq. 1 .and. iaffine_count .gt. 1) then
!write(*,*) bofnew,bofold,ip,'do we accept?...',NPvarexpandidx,iteration, bof1,csq,varerln,&
! bof2, bof3
!read(*,*)
!endif

!if ((ip .ge. NPthcentidx) .and. (ip .lt. NPthcentidx + NLC)) then
!write(*,*) iteration,ip,'TOP HAT INF',p(ip),pold,bofnew,bofold
!endif


!if ((ip .ge. NPsigexpandidx) .and. (ip .le. NPsigexpandidx+NPsigexpand)) then
!write(*,*) 'new sig expand test',p(ip),pold,bofnew,bofold
!endif

!do ip2 = 1,NP
!if ((pscale(ip2) .gt. 0) .or. (pversion(ip2) .eqv. .true.)) iplast = ip2
!enddo
!!write(*,*) ip,iplast
!if (ip .eq. iplast) then
!do ilc = 1,NLC
!write(ctit,'(I10.2)') ilc!! bug fix 16/09/2014 now removes old backups properly
!open(unit = 1,file = 'xinterptest_'//trim(adjustl(ctit))//'.dat')
!do it = lo(ilc),hi(ilc)
!write(1,*) t(it),x(it),er(it),xinterp(it)
!enddo
!close(1)
!enddo
!write(*,*) 'outputing interpolated model for testing to xinterptest_',iteration,ip
!!read(*,*)
!endif
!


if ( (( (( (bofnew < bofold) .or. (exp(-0.5*(bofnew-bofold)/T_ann) .gt. ran3(iseed))) &
.AND. (rejcosinc .eqv. .false.) .AND. (rejmdot .eqv. .false.) &
.AND. (rejth .eqv. .false.) .AND. (rejvar .eqv. .false.) .AND. &
(rejsig .eqv. .false.) .AND. (rejtrvisc .eqv. .false.) .AND. (rejtrirad .eqv. .false.)&
.AND. (rejaffine .eqv. .false.))&
.or. &
!exceptions for affine mode
((iaffinenow == 1) .and. (iaffine_count .lt. n_affine)   ) &
.or. &
!exceptions for chaos experiment
((chaos) .and. (altertf) .and. (iteration .gt. 100) .and. &
(rejcosinc .eqv. .false.) .AND. (rejmdot .eqv. .false.) &
.AND. (rejtrvisc .eqv. .false.) .AND. (rejtrirad .eqv. .false.) ) .or. &
(autoacc))&
.AND. (bofnew .eq. bofnew) .AND. (autorej .eqv. .false.) ) &
.or. (firstcall) ) then  !! the cosinc gt0 and lt 1 is imposed rather than reshaking the die until it is see update at top (17/7/14)






if (iaffinenow .ge. 1) then
write(*,*) ip,'affine parameter accepted...',iaffinenow,iaffine_count,n_affine,&
bofnew,bofold,p(ip),parsave(ip,iteration-1)
!write(*,*) p(NPumbhidx), p(NPcosincidx),'just checking move'
!else if (ip == NPumbhidx) then
!write(*,*) 'not affine mode mbh parameter',p(ip),pold,bofnew,bofold
endif
bofnew1=bof1
bofnew2=bof2
bofnew3=bof3   !! write down the components of the badness of fit
bofnew4=bof4
bofnew5=bof5
bofnew8=bof8

if (altertf) then  !! adjust the saved tf info
psinormold(:)=psinorm(:)
psisave(:,:) = psigrid(:,:)
endif


pRT_save(ip) = 0
autoacc = .false.
!!only update these things if not doing affine stepping or if the affine step is accepted
if (iaffinenow == 0 .or. iaffine_count .eq. n_affine) then
if (iaffine_count .eq. n_affine) pAT_affine = pAT_affine + 1
bofold=bofnew
pAT(ip)=pAT(ip)+1 ! record the move as accepted
account=account+1
!!update the rejected points
if (sigrej) then
do ilc = 1,NLC
sigrej_mark_save(lo(ilc):hi(ilc)) = sigrej_mark(lo(ilc):hi(ilc))
enddo
endif
cisqplot(1:NLC)=cisqnew(1:NLC)

!! all other priors
if (pricream) then
do ipc = 1,npricream
pricream_idxnow = pricream_idx(ipc)
if (pricream_idxnow == ip .or. firstcall) then ! .or.firstcall
pricream_bofold(ipc) = pricream_bofnew(ipc)
endif
enddo
endif

endif





do ilc=1,NLC
do it=1,Ntgrid
xgridplotsave(it,ilc)=xgridplot(it,ilc)
enddo
enddo
if (ip .eq. NPcosincidx) uinc=acos(p(ip))/pi*180
fursumold(:)=fursum(:)
xgridold(:)=xgrid(:)







if (varexpand .or. sigexpand) then !!parameters that need to be updated depending on whether step successful
!write(*,*) 'update ervar if needed',ip,NPsigexpand,NPsigexpandidx,NLC
if (((ip .ge. NPvarexpandidx) .and. (ip .lt. NPvarexpandidx + NPvarexpand)) .or.&
((ip .ge. NPsigexpandidx) .and. (ip .lt. NPsigexpandidx + npsigexpand))) then


varerlnold = varerln
if (yesxray) varerlnxrayold = varerlnxray

do ilc = 1,NLC
ilo = lo(ilc)
ihi = hi(ilc)
do it = ilo,ihi
ervarold(it) = ervar(it)
enddo
enddo

if (yesxray) then
do it = 1,NXray
erxrayvarold(it) = erxrayvar(it)
enddo
endif
endif
endif









!! update slowly varying background polynomial if accepted
if ((ip .ge. NPpolyidx) .and. (ip .lt. NPpolyidx + NPpolytot)) then
ilc_poly_now = floor(real(ip - NPpolyidx)/NPpoly) + 1
do it = 1,Ntgrid
bg_save(it,ilc_poly_now) = bg_now(it)
enddo
endif



!write(*,*) 'accept', ip,bofold,iteration
!!!REJECT!
!! reset the distribution if move rejected and ip is on the delay parms (no point otherwise)
!write(*,*) 'reject', ip,bofold,iteration
else



!if (iaffinenow .ge. 1) then
!write(*,*) 'affine parameter rejected...',iaffinenow,iaffine_count,n_affine,&
!bofnew,bofold,p(ip),parsave(ip,iteration-1)
!read(*,*)
!else if (ip == NPumbhidx) then
!write(*,*) 'not affine mode mbh parameter',p(ip),pold,bofnew,bofold
!endif


if (iaffinenow == 1) then
psinorm(:) = psinormold_affine(:)
xgrid(:)   = xgridold_affine(:)
fursum(:) = fursumold_affine(:)
psinormold(:) = psinormold_affine(:)
xgridold(:)   = xgridold_affine(:)
fursumold(:) = fursumold_affine(:)

else
psinorm(:)=psinormold(:)
xgrid(:)=xgridold(:)
fursum(:)=fursumold(:)
endif

do ilc = 2,NLC
if ((isharelc(ilc) == 1) .and. (samescale)) then
p(NPscaleidx+ilc-1) = p(NPscaleidx+ilc-2)
p(NPoffsetidx+ilc-1) = p(NPoffsetidx+ilc-2)
endif
enddo



!!! reset the error bars if they were altered
if (varexpand .or. sigexpand) then !!parameters that need to be updated depending on whether step successful
if (((ip .ge. NPvarexpandidx) .and. (ip .lt. NPvarexpandidx + NPvarexpand)).or. &
((ip .ge. NPsigexpandidx) .and. (ip .lt. NPsigexpandidx + npsigexpand))) then

if (iaffinenow == 1) then
varerln = varerlnold_affine
varerlnold = varerlnold_affine
if (yesxray) then
varerlnxray = varerlnxrayold_affine
varerlnxrayold = varerlnxrayold_affine
endif
else
varerln = varerlnold
if (yesxray) varerlnxray = varerlnxrayold
endif



do ilc = 1,NLC
ilo = lo(ilc)
ihi = hi(ilc)
do it = ilo,ihi
if (iaffinenow == 1)then
ervar(it) = ervarold_affine(it)
ervarold(it) = ervarold_affine(it)
else
ervar(it) = ervarold(it)
endif
enddo
enddo
if (yesxray) then
do it = 1,NXray
if (iaffinenow==1) then
erxrayvar(it) = erxrayvarold_affine(it)
else
erxrayvar(it) = erxrayvarold(it)
endif
enddo
endif
endif
endif








!!!!!!!!!!!!! re set the disk flux
if (iaffinenow == 1) then
do ilc=1,NLC
do it=1,Ntgrid
xgridplot(it,ilc)=xgridplotsave_affine(it,ilc)
xgridplotsave(it,ilc)=xgridplotsave_affine(it,ilc)
enddo
enddo
else
do ilc=1,NLC
do it=1,Ntgrid
xgridplot(it,ilc)=xgridplotsave(it,ilc)
enddo
enddo
endif

!!


if (iaffinenow == 1) then
do ip2 = 1,N_affine
ipnow = ip_affine(ip2)
p(ipnow) = parsave(ipnow,iteration-1)
enddo
if (iaffine_count .eq. N_affine) pRT_affine = pRT_affine + 1
else
p(ip)=pold
pRT(ip)=pRT(ip)+1 !! record the move as rejected
pRT_save(ip) = pRT_save(ip) + 1
endif


!!!! Only execute this block if the parameter affects the response function
if (altertf) then

!! reset the psi grid to the previous iteration and the lower and upper limits of the convolution
if (ip .eq. NPcosincidx) uinc=acos(p(NPcosincidx))*180/pi


do ilc=1,NLC
if (iaffinenow == 1) then
do itau=1,Ntaugrid
psigrid(itau,ilc)=psisave_affine(itau,ilc)
psisave(itau,ilc) = psisave_affine(itau,ilc)
enddo
else
do itau=1,Ntaugrid
psigrid(itau,ilc)=psisave(itau,ilc)
enddo
endif

if (version .eq. 1) then
if (iaffinenow == 1) then
fdisk(1:NLC)=fdisksave_affine(1:NLC)
fdisksave(1:NLC) = fdisksave_affine(1:NLC)
else
fdisk(1:NLC)=fdisksave(1:NLC)
endif
endif

!!update 13/12/2016
if (isharelc(ilc) .eq. 1) then
itaumax(ilc) = itaumax(ilc - 1)
else
itaumax(ilc) = Ntaugrid
do itau = 2,Ntaugrid !! new condition 3/04/2015
!find max non 0 index of tf for omp parallel convolving later
!(omp does not like exit statements within its loops so we deal with that here)
if (psigrid(itau,ilc) .le. 0 .and. psigrid(itau-1,ilc) .gt.0) then
itaumax(ilc) = itau
exit
endif
enddo
endif

if (isharelc(ilc) .eq. 1) then
itaumin(ilc) = itaumin(ilc-1)
else
itaumin(ilc) = 1
do itau = 1, itaumax(ilc) - 2
if (psigrid(itau,ilc) .gt. 0) then
exit
endif
enddo
itaumin(ilc) = itau
endif
!!update 13/12/2016



enddo

!!

!!! Only execute this block if the irad and visc  tr power law indices laws are stepped together
if ((step_visc_irad_together) .and. (ip .eq. NPtridx)) then
pscale(NPtridx+1) = pscale(NPtridx)
p(NPtridx+1) = p(NPtridx)
endif
!!!

endif
!!!!

if (rejth) rejth =.false.
if (rejvar) rejvar =.false.
if (rejsig) rejsig =.false.
if (rejtrvisc) rejtrvisc =.false.
if (rejtrirad) rejtrirad =.false.
if (rejcosinc) rejcosinc =.false.
if (rejmdot)   rejmdot =.false.
if ((rejaffine) .AND. (iaffine_count .ge. n_affine) .AND. (iaffinenow == 1)) rejaffine =.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif


bofrejectnew(ip)=bofold

firstrun=.false.
firstip=.false.

bofsave_old(1,ip) = bofnew1
bofsave_old(2,ip) = bofnew2
bofsave_old(3,ip) = bofnew3
bofsave_old(4,ip) = bofnew4
bofsave_old(5,ip) = bofnew8

bofsave_new(1,ip) = bof1
bofsave_new(2,ip) = bof2
bofsave_new(3,ip) = bof3
bofsave_new(4,ip) = bof4
bofsave_new(5,ip) = bof8







if (iaffinenow .ge. 1) iaffine_count = iaffine_count + 1

!if (iaffinenow .ge. 1) then
!write(*,*) iaffinenow,iaffine_count,ip,'affine stepping',ip_affine(1:n_affine)
!read(*,*)
!endif


enddo  !! end the ip loop



!write(*,*) pversion(NPumbhidx),pversion(NPcosincidx),'version check'
!read(*,*)

if (iaffinenow ==1) then

do ip2 = 1,N_affine
write(*,*) 'affine summary'
write(*,*) parsave(ip_affine(ip2),iteration-1),p(ip_affine(ip2)),rejaffine
write(*,*) rejtrirad,rejtrvisc,rejvar,rejsig,rejcosinc,rejmdot
enddo
!if (account .gt. 0) read(*,*)
endif


!do it = 1,Ntaugrid
!write(*,*) taugrid(it), psigrid(it,1:NLC)
!enddo
!write(*,*) 'cents', p(NPthcentidx:NPthcentidx+NLC-1)
!write(*,*) 'just before plot'
!read(*,*)










if (timeon) call system_clock(itimelo)

!!! END OF PART 2  (Now outside the iteration)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


bofsave(iteration)=bofold



!!!!  write the badness of fit components to a text file for inspection

open(unit=1,file='testbofs.dat',access='append')

write(1,*) bofnew1,bofnew2,bofnew3,bofnew4,bofnew5,bofnew8, bof9
!write(*,*) bofnew1,bofnew2,bofnew3,bofnew4,bofnew5
call getcwd(cwd)
!write(*,*) trim(cwd),': - cwd'
!read(*,*)
close(1)
!!!!
if (timeon) then
call system_clock(itimehi)
itimesavebof = itimehi-itimelo
endif



if (timeon) call system_clock(itimelo)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Burnin checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! test for burn in end after 'iburnincheck' iterations
iburnincheck=200
iburningap = 20!max(1,iburnin(1)) !leave a gap equal to size of 1st burnin before checking again
if (iteration .eq. 1) ibnow = 1
if (iteration .gt. iburnincheck ) then
!!!!

do ib = 1,Nburnin					!! identify which stage burnin we are ib at iteration iburnin with previous burn in
!  	write(*,*) ib,ibnow,'hh'							 ! at iburninprev
iburninnow = iburnin(ib)
if (ib .gt. 1) then					!! identify the previous burnin iteration
iburninprev = iburnin(ib-1)
else
iburninprev = 1
endif

if (iburnin(ib) .eq. 0) then
ibnow = ib
!write(*,*) ibnow
exit
endif
enddo
!!!!
if (ib .eq. Nburnin+1) ib = Nburnin ! make it fine for the last point


!!!! calculate the median (first burnin) or mean (subsequent burnins) (for speed)

Nme = iteration-iburninprev+1
if (ib .eq. 1) then
bofme(ib) = med(bofsave(1:iteration),Nme)
else
bofme(ib) = avg(bofsave(iburninprev:iteration),Nme)
endif
!!!!


if (iburnin(1) .gt. 0) then
iburninnow = iteration
else
iburninprev = 0
endif

if ((ib .lt. Nburnin) .and. (mod(iteration,iburningap) .eq. 0) ) then
if (ib .eq. 1) then
if (bofsave(iteration) .gt. bofme(ib)) iburnin(ib)= iteration
else
if (bofme(ib-1) .gt. bofme(ib)) iburnin(ib) = iteration
endif

!write(*,*) 'IBURN IN ARGHHHHH', ib


endif

endif



!if (iteration .gt. 4) iburnin2 = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! end of burnin checks






parsave(1:NP,iteration) = p(1:NP)

!if (iteration .gt. 4) iburnin2 = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! end of burnin checks


!! update error envelope parms
if ((ib .lt. Nburnin) .and. (iburnin(ibnow) .gt. 0)) then !reset the envelopes
senv1drive(:) = 0.d0
senv2drive(:) = 0.d0
senv1tf(:,:)  = 0.d0
senv2tf(:,:)  = 0.d0
senv1(:,:)    = 0.d0
senv2(:,:)    = 0.d0
senv1p(:)     = 0.d0
senv2p(:)     = 0.d0
senv1_bg(:,:)    = 0.d0
senv2_bg(:,:)    = 0.d0

else


do it=1,Ntgrid   !!!! driving lightcurve
a=xgrid(it)
senv1drive(it)=senv1drive(it)+a
senv2drive(it)=senv2drive(it)+a*a
enddo



itsubburn=iteration-iburninprev


!write(*,*) 'quick conv',quick_conv_itp,interpidxmin,interpidxmax
! if we are using the tf quick mode then have to recalculate model for plotting and saving
if (quick_conv_itp) then
call cream_modop_ps(parsave(1:NP,iteration:iteration),wavem,embhref,redshift,&
Nlc,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,NPscaleidx,NPthcentidx,NPtridx,&
1,Ntgrid,Ntaugrid,w,taugrid,tgrid,xgridplot,sigechosave)
!do it = interpidxmin,interpidxmax
!write(*,*) it, tgrid(it), xgridplot(it,1:NLC),'aarra'
!enddo
!
!read(*,*)
endif


do ilc=1,NLC	 !!!! Echo Lightcurve offset
do itau=1,Ntaugrid  !!! update the error envelop[es for the transfer functions
a=psigrid(itau,ilc)
senv1tf(itau,ilc)=senv1tf(itau,ilc)+a
senv2tf(itau,ilc)=senv2tf(itau,ilc)+a*a
enddo


if (bgvary) then
do it = 1,Ntgrid
b = bg_save(it,ilc)
senv1_bg(it,ilc)=senv1_bg(it,ilc)+b
senv2_bg(it,ilc)=senv2_bg(it,ilc)+b*b
enddo
endif





!do it=1,Ntgrid   !!!! Echo Lightcurves
!if (quick_conv_itp) then
! a = echosave(it,ilc)
!else
! a = xgridplot(it,ilc)
!endif

! if we have a bad iteration then just set the disk light curve to the existing mean
do it = 1,Ntgrid
a = xgridplot(it,ilc)
if (a .ne. a) then
a = senv1(it,ilc)/itsubburn
endif

senv1(it,ilc)=senv1(it,ilc)+a
senv2(it,ilc)=senv2(it,ilc)+a*a
enddo


enddo


do ip = 1,NP !! all other parameters


a = p(ip)

b1 = senv1p(ip) + a
b2 = senv2p(ip) + a*a

senv1p(ip) = b1
senv2p(ip) = b2
sdparm(ip) = sqrt( abs((b2 - b1*b1/itsubburn)/(itsubburn) ) )
avparm(ip) = b1/itsubburn
enddo




!a = p(NPoffsetidx)
!
!b1 = senv1p(NPoffsetidx) + a
!b2 = senv2p(NPoffsetidx) + a*a
!write(*,*) iteration, iburninprev,b1,b2, senv1p(NPoffsetidx)
!stop
!write(*,*) ibnow
!write(*,*) iburnin(1:ibnow)
!write(*,*) b1/itsubburn,sqrt( abs((b2 - b1*b1/itsubburn)/(itsubburn) ) )
!
!nl = 200000!numline('outputpars.dat') -1
!open(unit = 1,file='outputpars.dat')
!allocate(osnew(nl))
!
!b1new = 0.d0
!b2new = 0.d0
!idx = 0
!itestlo = 1800
!do i = 1,nl
!read(1,*) crap,crap,crap,crap,crap,crap,crap,crap,crap,osnew(i)
!
!if (i .ge. itestlo ) then
!
!idx = idx + 1
!b1new = b1new + osnew(i)
!b2new = b2new + osnew(i)**2
!if (idx .lt. 1000) write(*,*) idx, b1new, b2new
!
!endif
!
!enddo
!close(1)
!sdnew = sdev(real(osnew(iburninprev+1:iteration)),itsubburn)
!avnew = avg(real(osnew(iburninprev+1:iteration)),itsubburn)
!write(*,*) avnew,sdnew,b1new/idx,(abs(b2new/idx - (b1new/idx)**2))**0.5!sqrt( abs((b2new - b1new*b1new/itsubburn)/(itsubburn) ) )
!write(*,*) itsubburn, iteration, iburninprev, nl,idx,b2new,b1new
!
!deallocate(osnew)
!stop


!write(*,*) sqrt( abs((b2 - b1*b1/itsubburn)/(itsubburn) ) )
!write(*,*) sqrt( abs((b2 - b1*b1/(itsubburn+1))/((itsubburn+1)) ) )
!write(*,*) sdparm(Npoffsetidx), avparm(NPoffseetidx)
!stop


!special case for inclination
rln10 = alog(10.0)
radincplot=avparm(NPcosincidx)

avcosinc = avparm(NPcosincidx)
sdcosinc=sdparm(NPcosincidx)
sd_inc=abs( 1/sin(radincplot) ) * sdcosinc *180/pi
av_inc=acos(avcosinc)*180/pi!acos( avcosinc / (1.+0.5*sd_inc*sd_inc) )*180/pi
!avrad_inc = acos(sdparm(NPcosincidx))
!av_inc  = avrad_inc*180/pi
!sd_inc  = abs( 1/sin(avrad_inc) ) * sdparm(NPcosincidx) *180/pi



endif








































!!!!!!
if (timeon) then
call system_clock(itimehi)
itimesaveenv = itimehi-itimelo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!! calculate mean delay

do ilc=1,NLC
bot=0.d0
top=0.d0
if ((wavobs(ilc) .ge. 0) .and. (wavobs(ilc) .le. 100)) then
taumean(ilc) = 0
else
do itau=1,Ntaugrid
del=taugrid(itau)
tf=psigrid(itau,ilc)
top=top+del*tf
bot=bot+tf
enddo
taumean(ilc)=top/bot
endif
enddo






!!!!!!!!!!!!!!!!!!!!!!! Plot model power spectrum
!!

alphatemp = p(NPpspecidx+2)
betatemp  = p(NPpspecidx+3)
p0temp    = p(NPpspecidx)
df        = dw/twopi


if (break) then
do it=1,NWres
a = twopi*fres(it)/p(NPpspecidx+1)
psres(it)=p0temp* dw * a**(-1*alphatemp) /(1.0 +  a**(betatemp-alphatemp))  ! new change august 8 14 df = dw*0.5/pi
enddo
else
do it=1,NWres
a = twopi*fres(it)/p(NPpspecidx+1)
psres(it)=p0temp* dw /(a**P(NPpspecidx+2))  !! new change august 8 14 df = dw*0.5/pi
!write(*,*) p0temp,dw,p(NPpspecidx+1),p(NPpspecidx),p0temp* dw,(a**P(NPpspecidx+2)),'help'
enddo
endif




if (iburnin2 .gt. 0) then
do it = 1,NWres
a=psres(it)
senv1fur(it)=senv1fur(it)+a
senv2fur(it)=senv2fur(it)+a*a
enddo
endif



if (iteration .eq. 1 .or. firstcall) then
do iw=1,NWres
write(*,*) fres(iw),psres(iw)
enddo
endif




!do it=1,Ntgrid !help
!write(*,*) tgrid(it),xgrid(it)
!enddo
!read(*,*)

if (timeon) call system_clock(itimelo)


call chdir('./plots') ! change to plot directory

if (plotsave) then

!!new 10jun2016 save correlation plot and correlation matrix for offset, stretch and error bar expansion parameters
if (pythonplot) then
if (sigexpand) then
call system('cp ../../../cream_python/pycov_f_stretch_offset.py ./')
call system('python pycov_f_stretch_offset.py')
call system('rm pycov_f_stretch_offset.py')
endif

if (varexpand) then
call system('cp ../../../cream_python/pycov_f_var.py ./')
call system('python pycov_f_var.py')
call system('rm pycov_f_var.py')
endif
endif
!


!


call superfitplot(tgrid,interpidxmin,interpidxmax,xgrid,&
taugrid,psigrid,Npoints,t,x,er,p,pscale,psafe,w,iteration)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! extra bit for plotting posterior histograms
if (pscale(NPcosincidx) .ne. 0) then
ifloatinc_pp = 1
else
ifloatinc_pp = 0
endif

if (pscale(NPumbhidx) .ne. 0) then
ifloatmmdot_pp = 1
else
ifloatmmdot_pp = 0
endif

if (pscale(NPtridx) .ne. 0) then
ifloattrvisc_pp = 1
else
ifloattrvisc_pp = 0
endif



!!write(*,*) 'about to call postplot',iteration
filepath_pp = '../outputpars.dat'

!write(*,*) 'before post plot'

!call postplot(iteration, filepath_pp, ifloatinc_pp, ifloatmmdot_pp, ifloattrvisc_pp)




!stop



!! allocate the array as required by flux flux program
!


!!! plot the flux flux diagram hell
if (skipfvg .eqv. .false.) then

idxwav0 = NLC/2
Ndatmax = hi(NLC)
allocate(dattemp(Ndatmax,NLC),sigtemp(Ndatmax,NLC),idxlc(NLC),t_dattemp(Ndatmax,NLC))




do ilc = 1,NLC
Ntemp      = hi(ilc) - lo(ilc) + 1
idxlc(ilc) = Ntemp
fexp       = p(NPsigexpandidx+ilc - 1)
do it = 1,Ntemp
t_dattemp(it,ilc) = t(lo(ilc)+it-1)
dattemp(it,ilc)   = x(lo(ilc)+it-1)
sigtemp(it,ilc)   = fexp*er(lo(ilc)+it-1)
!write(*,*) it,ilc,fexp, NPsigexpandidx, p(NPsigexpandidx+ilc - 1),er(lo(ilc)+it-1),'before routine'
enddo
!read(*,*)
enddo


!! subroutine here to sort all the light curves into a temporary array that has number of wavelengths i.e 19 rather than number of telescopes (64) each with the correct expansion factor (basically gives correct expansion factor when splitting light curves up into multiple telesopes. Can then pass these arrays onto the flux flux code)
call wav_unique(wavem,NLC,Nwav_fvg)
allocate(ttemp2(Ndatmax,Nwav_fvg), xtemp2(Ndatmax,Nwav_fvg),&
sigtemp2(Ndatmax,Nwav_fvg), wav2(Nwav_fvg),idxlc2(Nwav_fvg),respcent2(Nwav_fvg),&
respsig2(Nwav_fvg),ilclo2(Nwav_fvg))

!write(*,*) 'hererererere fvg_prep', Nwav_fvg,NLC
!stop




if (NLC .ne. Nwav_fvg) then
call fvg_prep_multitel(NLC, Nwav_fvg, wavem, idxlc,Ndatmax, t_dattemp, dattemp,&
sigtemp, ttemp2, xtemp2, sigtemp2, wav2,idxlc2,ilclo2)

do iw = 1,Nwav_fvg
inow = ilclo2(iw)
respcent2(iw) = respcent(inow)
respsig2(iw)  = respsig(inow)
enddo



endif

!
if (iburnin2 .gt. 0) then
ciave = avparm(NPcosincidx)
cisd  = sdparm(NPcosincidx)
else
ciave = p(NPcosincidx)
cisd  = ciave/100
endif



open(unit=1,file='testdump.dat')
write(1,*) Ndatmax,min(NLC,Nwav_fvg)
write(1,*) avcosinc, sdcosinc

if (NLC .ne. Nwav_fvg) then
write(1,*) idxlc2(1:Nwav_fvg)
write(1,*) respcent2(1:Nwav_fvg)
write(1,*) respsig2(1:Nwav_fvg)
else
write(1,*) idxlc(1:NLC)
write(1,*) respcent(1:NLC)
write(1,*) respsig(1:NLC)
endif

if (NLC .ne. Nwav_fvg) then
write(1,*) wav2(1:Nwav_fvg)
else
write(1,*) wavem(1:NLC)
endif

!do i = 1,NLC
!write(*,*) wavem(i),respcent(i),respsig(i)
!enddo
!stop
if (NLC .ne. Nwav_fvg) then
do it = 1,Ndatmax
write(1,*) ttemp2(it,1:Nwav_fvg), xtemp2(it,1:Nwav_fvg), sigtemp2(it,1:Nwav_fvg)
enddo
close(1)
else
do it = 1,Ndatmax
write(1,*) t_dattemp(it,1:NLC), dattemp(it,1:NLC), sigtemp(it,1:NLC)
enddo
close(1)
endif

!!!! obtain all the information required for the difference spectrum fitting code !!!

! index of max and min
xmax_diff = maxval(xgrid)
xmin_diff = minval(xgrid)

do it = 1,Ntgrid
xgnow = xgrid(it)
if (xgnow .eq. xmax_diff) idxmax_diff = it
if (xgnow .eq. xmin_diff) idxmin_diff = it
enddo
tmin_diff = tgrid(idxmin_diff)
tmax_diff = tgrid(idxmax_diff)



if (itsubburn .gt. 1) then
open(unit=1,file='CREAM_diff_flux.dat')
write(1,*) NLC, tmin_diff, tmax_diff, xmin_diff, xmax_diff
do ilc = 1,NLC
xmin_diff = senv1(idxmin_diff,ilc) / itsubburn
sigmin_diff = sqrt( abs( (senv2(idxmin_diff,ilc) - senv1(idxmin_diff,ilc)&
*senv1(idxmin_diff,ilc)/itsubburn)/(itsubburn) ) )

xmax_diff = senv1(idxmax_diff,ilc) / itsubburn
sigmax_diff = sqrt( abs( (senv2(idxmax_diff,ilc) - senv1(idxmax_diff,ilc)&
*senv1(idxmax_diff,ilc)/itsubburn)/(itsubburn) ) )

write(1,*) wavem(ilc), xmin_diff, sigmin_diff, xmax_diff, sigmax_diff
enddo
close(1)
endif

write(*,*) 'before fvgplot'
write(*,*) NLC, NPsigexpandidx,NPsigexpandidx+NLC-1, Ndatmax
!call cream_fvgplot(NLC, Ndatmax, t_dattemp, dattemp, sigtemp, respcent, respsig, wavem,&
!p(NPsigexpandidx:NPsigexpandidx-1),idxwav0,ciave,cisd, idxlc,ebmvmw, redshift, omega_m, omega_l)
write(*,*) 'end of fvg plot'
deallocate(dattemp,sigtemp,t_dattemp,idxlc)
deallocate(ttemp2, xtemp2,sigtemp2, wav2,idxlc2,respsig2,ilclo2,respcent2)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



!write(*,*) 'end of plot save'




!call mypgline(par,nits,fname,npar,axlab,ymin_in,ymax_in)


else if (plotscreen) then

!write(*,*) 'before plot', maxval(psigrid(1:Ntaugrid, 1)), NLC
!read(*,*)


call superfitplot(tgrid,interpidxmin,interpidxmax,xgrid,&
taugrid,psigrid,Npoints,t,x,er,p,pscale,psafe,w,iteration)



endif









call chdir('./..') !! return to

if (timeon) then
call system_clock(itimehi)
itimeplot = itimehi - itimelo
endif






!!! update the luminosity distance
if (version .eq. 1) then
do ilc=1,NLC
dl(ilc)=dlflatm(p(NPumbhidx)*embhref,wavem(ilc),p(NPcosincidx),fvb(ilc),fvf(ilc))
enddo
else
dl(1:NLC) = 0
endif




!!!!!! backing up eveything
if (timeon) call system_clock(itimelo)



!!!!!!!!!! save the current parameter, stepsize, logstep, on or off
open (unit = 1, file = 'cream_resume.dat')
do iw = 1,NP
write(1,*) p(iw), pscale(iw), psafe(iw),pversion(iw)
enddo
close(1)




!!!!! if tfb_x on save the accretion disc parameters here
if (i_tfx == 1) then
open(unit=1,file='outputpars_tfxdisk.dat',access='append')
write(1,*) p(NPtridx:NPtridx+NPtr-1), p(NPcosincidx)
close(1)
else

!!!!!!!!!! save the iteration number, bof, MMdot, and inclination parameters
open(unit=1,file='outputpars.dat',access='append')
if (version .eq. 1) then
write(1,*) iteration, bofold, p(NPumbhidx),P(NPcosincidx),p(NPdlidx),dl,&
p(NPsigexpandidx:NPsigexpandidx+Npsigexpand-1)*sd_overflow(1:NLC),&
p(NPoffsetidx:NPoffsetidx+NLC-1)+ave_overflow(1:NLC),&
p(NPfurconstidx),taumean
else if (version .eq. 3) then
write(1,*) iteration, bofold, p(NPumbhidx),P(NPcosincidx),p(NPtridx),p(NPtridx+1),taumean

!dl,&
!p(NPfurconstidx),
endif
close(1)
!!!!!!!!!!
endif





!!!!!!!!!! save the scale factors
open(unit=1,file='outputpars2.dat',access='append')
write(1,*) p(NPscaleidx:NPscaleidx+NLC-1)*sd_overflow(1:NLC),&
p(NPoffsetidx:NPoffsetidx+NLC-1)+ave_overflow(1:NLC),&
p(NPsigexpandidx:NPsigexpandidx+NLC-1)
close(1)
!!!!!!!!!!


!!!!!!!!! save the variance expansion parameters
if (varexpand) then
open(unit=1,file='outputpars_varexpand.dat',access='append')
write(1,*) p(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1)*&
sd_overflow(1:NLC)**2
close(1)
endif


!!!!!!!!! save the top hat parameters if they are varied
ilccount = 0
do ilc = 1,NLC
if (pscale(NPthcentidx+ilc-1) .gt. 0 .or. pscale(NPthcentidx+NLC+ilc-1) .gt. 0) then
ilccount = ilccount + 1
else
exit
endif
enddo
!if (ilccount .gt. 0) then
open(unit=1,file='outputpars_th.dat',access='append')
write(1,*) p(NPthcentidx:NPthcentidx+2*NLC-1)
close(1)
!endif

!!!!!!!!!! save the scale factors
open(unit=1,file='cream_fluxflux.dat')

write(1,*) NPoffsetidx,NLC,NPscaleidx
il=1
wavold = wavem(1)
do ip = NPscaleidx,NPscaleidx+NLC-1
wavnew = wavem(il)
if ((wavnew .ne. wavold) .or. (il .eq. 1)) then
write(t0,'(F10.2)') wavnew
write(t1,'(F10.2)') avparm(ip+NLC)
write(t2,'(F10.2)') sdparm(ip+NLC)
write(t3,'(F10.2)') avparm(ip)
write(t4,'(F10.2)') sdparm(ip)
t_ff='$'//trim(adjustl(t0))//'$ & $'//trim(adjustl(t1))//'\pm'
t_ff=trim(adjustl(t_ff))//trim(adjustl(t2))//'$ & $'//trim(adjustl(t3))//'\pm'
t_ff=trim(adjustl(t_ff))//trim(adjustl(t4))//'$'

write(1,*) trim(adjustl(t_ff))
wavold = wavnew
endif

il = il + 1
enddo
close(1)
!!!!!!!!!!


open(unit=1,file='outputpars3.dat')
idx = 1
do ip = NPsigexpandidx,NPsigexpandidx+NLC-1
write(1,*) wavobs(idx), avparm(ip), sdparm(ip)
idx = idx + 1
enddo
close(1)





!!!!!!!!! write the filename to a file in output
if (iteration .eq. 1) then
open(unit=1,file='fileinfo.dat',status='replace')
write(1,*) filename2
close(1)
endif



!if (iteration .eq. 1) then
!call chdir('./../..')
!open(unit=1,file='plotpspecpar.par',status='replace')
!write(1,*) NLC
!write(1,*) wavem
!write(1,*) trim(dirworking)
!close(1)
!call chdir(trim(dirworking))
!endif







mi = mod(iteration, abs(iplotsave))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! evaluate the co-variance matrix and save all the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (mi .eq. 0) then
call covcalc(pscale,pversion,iteration,4)
endif
!!!!!!!!!!!!!!!!!!!!!! covsum is now iteration x the covariance matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!










!!!!!!!!!!!! write data to backup file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! write in the data to the backup file
!! call syst
!write(*,*) 'Dont stop me now!... saving parameters '
!write(*,*) bofold,umdot,uhx,eff,alb,urin,uinc,diskk,diskka,alphatemp,alphamean,sigalphasquare,dw,w0temp,p0temp,dtgrid



if (mi .eq. 0) then



!!!! new 27november save the mmdot and inclination parameters
open(unit = 1,file='cream_tfparms.dat')
write(1,*) avparm(NPumbhidx), sdparm(NPumbhidx)
write(1,*) avparm(NPcosincidx), sdparm(NPcosincidx)
close(1)

!!!! new 12jun save fourier parms
p0temp    = p(NPpspecidx)							             !! store the pspec defaults
w0temp    = p(NPpspecidx+1)
alphatemp = p(NPpspecidx+2)
betatemp  = p(NPpspecidx+3)


open(unit = 2, file = 'cream_gvalues.dat')
open(unit = 1,file='cream_furparms.dat')

write(2,*) 'angular_freq	sdprior	gsk	gck'
idw = 1
write(1,*) p0temp,dw,w0temp
!compute the g value based on the Fourier prior
do iw = 1,NW
a = w(iw)/w0temp
if (break) then
sdprior = p0temp*dw * (a**(-1*alphatemp))  / (1. +   a**(betatemp-alphatemp)) * 0.5
else
sdprior = p0temp*dw/(a**alphatemp) * 0.5
endif
gsk = sdprior**2 / (sdprior**2 + sdparm(idw)**2)
gck = sdprior**2 / (sdprior**2 + sdparm(idw+1)**2)
write(1,*) w(iw), avparm(idw),avparm(idw+1), sdparm(idw), sdparm(idw+1)
write(2,*) w(iw),sdprior, gsk, gck

idw = idw+2
enddo
close(1)
close(2)



!!!! new 22aug save the offset pamrs
open(unit = 1,file='cream_osparms.dat')
do ilc3 = 1,NLC
write(1,*) avparm(NPoffsetidx+ilc3-1), sdparm(NPoffsetidx+ilc3-1)
enddo
close(1)

!!!! new 22aug save the offset pamrs
open(unit = 1,file='cream_stretchparms.dat')
do ilc3 = 1,NLC
write(1,*) avparm(NPscaleidx+ilc3-1), sdparm(NPscaleidx+ilc3-1)
enddo
close(1)


!!!! new 29Jan save the reduced chi squared information
open(unit = 1,file = 'cream_redcisq.dat')
do ilc = 1,NLC
if (yesxray) then
f = p(NPsigexpandidx+ilc)
else
f = p(NPsigexpandidx+ilc-1)
endif
Ndat=hi(ilc)-lo(ilc)+1
write(1,*) wavem(ilc),cisqplot(ilc)/f/f , f, Ndat
enddo
close(1)
!!!!





!! remove old backups if have more than nbackup
!Only do the backup thing if backupon is false

if (backupon) then
if (iteration .gt. nbackuptimesiplotsave) then
write(ctit,'(I10.2)') iteration-(nbackuptimesiplotsave-iplotsave) !! bug fix 16/09/2014 now removes old backups properly
ctit=adjustl(ctit)
call system('rm outputbackup_'//trim(ctit)//'.bc')

write(*,*) 'outputbackup_'//trim(ctit)//'.bc'
endif

!! edit the creamcontpar.par file to the name of the most recent backup
write(ctit,'(I10.2)') iteration
ctit=adjustl(ctit)
call chdir('./../..')
open(unit=1,file='creamcontpar.par',status='replace')
write(1,*)trim(trim(dirworking))
write(1,*)trim('outputbackup_'//trim(ctit)//'.bc')
close(1)
call chdir(trim(dirworking))

!! Make the new backup


if (firstcall .eqv. .false.) then
open(unit=1,file='outputbackup_'//trim(ctit)//'.bc',status='replace')

write(1,*) filename2
write(1,*) nits
write(1,*) yesxray, noconvolve
write(1,*) sigexpand
write(1,*) version
write(1,*) iteration
write(1,*) iplotsave
write(1,*) positivity,pubplot, BIGsave
write(1,*) nbackuptimesiplotsave,nbackup
write(1,*) bofold
write(1,*) ebmvMW,redshift, omega_m, omega_l
write(1,*) NLC,NT,NTgrid,Ntaugrid,NPF,NPpspec,NPpspecidx,NPscaleidx,NPscale,NPdelay,NPumbhidx,NPcosincidx,NPdlidx,NW,NP,RT,&
AT,NPoffset,NPoffsetidx,NPtr,NPtridx,NPgal,NPgalidx,Nxray,NPfurconstidx,NPsigexpand,NPsigexpandidx,&
NPebmvidx,NPebmv,Nburnin
write(1,*) embhref,uhx,eff,alb,urin,ur0,uinc,diskk,diskka,alphatemp,alphamean,sigalphasquare,&
betatemp, betamean, sigbetasquare, dw,w0temp,p0temp,&
dtgrid,interpidxmin,interpidxmax,sigp0square,sigw0,p0mean,w0mean,break,NWres,tgridloextra,tgridhiextra,&
pri_mmdotmean,sigmmdotsquare,alpha_cosinc, cosinc_cut

write(1,*) iburnin
write(1,*) bofme
write(1,*) cisqplot(1:NLC)
write(1,*) fvb(1:NLC)
write(1,*) fvf(1:NLC)

write(1,*) !! break in file
write(1,*) !!

write(1,*) wavobs(1:NLC)
write(1,*) !! break in file
write(1,*) !!
write(1,*) w(1:NW)
write(1,*) !! break in file
write(1,*) !!
do i=1,NP !! read in parms
write(1,*) p(i),pscale(i),psafe(i),pversion(i)
enddo
write(1,*) !! break in file
write(1,*) !!
do i=1,NT !! read data
write(1,*) t(i),x(i),er(i)
enddo
write(1,*)
write(1,*)
write(1,*) bofsave(1:iteration)
write(1,*) !! break
write(1,*) tgrid(1:Ntgrid)
write(1,*) xgridold(1:Ntgrid)
write(1,*) fursumold(1:NTgrid)
write(1,*) !! break
write(1,*) taugrid(1:Ntaugrid)
write(1,*) !! break
write(1,*) taugrididx(1:Ntaugrid)
write(1,*) !!break
write(1,*) !!break
!!! write in the sin and cosine parameters
!!!
write(1,*) !break
write(1,*) !break
!!! now write in the psi grid
do i=1,NLC
write(1,*)psigrid(1:Ntaugrid,i)
enddo
!!!
write(1,*) !break
write(1,*) !break
!!! now write in the xgridplot save information
do ilc=1,NLC
write(1,*) xgridplotsave(1:Ntgrid,ilc)
enddo

write(1,*) !break
write(1,*) !break
write(1,*) interpidx(1:NT)
write(1,*) !break
write(1,*) !break
write(1,*) Npoints(1:NLC+1)
write(1,*) lo(1:NLC)
write(1,*) hi(1:NLC)

write(1,*) xraymed ! set to zero if using no xray data
if (yesxray) then
write(1,*) txray
write(1,*) xxray
write(1,*) erxray
write(1,*) xrayinterpidx
endif

!! write in the stuff for the error envelope plotting
do ip = 1,NP
write(1,*) senv1p(ip),senv2p(ip)
enddo

do ilc=1,nlc
write(1,*) senv1(1:Ntgrid,ilc)
write(1,*) senv2(1:Ntgrid,ilc)
write(1,*) senv1tf(1:Ntaugrid,ilc)
write(1,*) senv2tf(1:Ntaugrid,ilc)
enddo
write(1,*) senv1drive(1:Ntgrid)
write(1,*) senv2drive(1:Ntgrid)
write(1,*) senv1fur(1:NWres)
write(1,*) senv2fur(1:NWres)

close(1)
endif ! end if firstcall .eq. false
endif !endif backupon

! new 30/04/205 write all previous parameters to a file for use when evaluating covariances
if (BIGsave) then
open(unit=1,file='CREAM_allpars_BIG.dat',status='replace')
do it = 1,iteration
write(1,*) parsave(1:NP,it)
enddo
close(1)
endif





!!now save the expansion parameters

open(unit = 1,file ='cream_sigexpand_parms.dat')
write(1,*) 'columns: ffactor, sigffactor, mean error times f, error in mean error timesf'
do ilc = 1,NLC

ersum = 0.d0
nlcnow = hi(ilc) - lo(ilc) + 1
do it = lo(ilc),hi(ilc)
ernow = er(it)
ersum = ernow + ersum
enddo
ermean = ersum / nlcnow ! calculate average error bar
ersd = 0.d0
do it = lo(ilc),hi(ilc)
ernow = er(it)
a = (ernow - ermean)
ersd = ersd + a*a
enddo
ersd = sqrt(ersd/nlcnow)
ermean_sd = ersd/sqrt(real(nlcnow))


fnow   = avparm(NPsigexpandidx+ilc-1)
fnowsd = sdparm(NPsigexpandidx+ilc-1)


sig_times_f = ermean * fnow
sig_sig_times_f = sig_times_f * sqrt( (fnowsd/fnow)**2 + (ermean_sd/ermean)**2 )

write(1,*) fnow,fnowsd, sig_times_f, sig_sig_times_f

enddo

close(1)


endif  !!!! end backing up and deleting old backups

firstcall =.false.



if (timeon) then !! time check for saving parms
call system_clock(itimehi)
itimesavebackup = itimehi - itimelo
endif


!write(*,*) 'Finished backing up'!,med(xgridplot(1:Ntgrid,1),Ntgrid),offset,p(NPoffsetidx)

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

















!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
! use the new python plotting script. I think you still need to run the superfit plot to make the plotfiles
if (plotsave) then
call chdir('./plots')
call system('cp ../../../cream_python/pyplot_cream_easy_2.py ./')
call system('cp ../../../cream_python/pyplot_parm.py ./')
call system('cp ../../../cream_python/creamconvert.py ./')
call system('cp ../../../cream_python/creamfvg.py ./')
call system('cp ../../../cream_python/cream_pythonsubs.py ./')

call system('python creamconvert.py')
if (pythonplot) then
call system('python pyplot_cream_easy_2.py')
call system('python pyplot_parm.py')
call system('python creamfvg.py')
endif

call system('rm creamconvert.py')
call system('rm creamfvg.py')
call system('rm pyplot_cream_easy_2.py')
call system('rm pyplot_parm.py')
call system('rm cream_pythonsubs.py')
call chdir('../')
endif
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!




























call system_clock(iendtime)
secs_it = (iendtime-istarttime)*0.001



write(*,*) 'Iteration:',iteration, ' Secs/It:',secs_it,' It/day:',3600*24/secs_it
write(*,*) 'bof:', bofold

if (bgvary) then
write(*,*) 'Background polynomials at ref time:',tref_poly
write(*,*) p(NPpolyidx:NPpolyidx+NPpolytot-1)
write(*,*) pscale(NPpolyidx:NPpolyidx+NPpolytot-1)
endif

write(*,*)
do ilc = 1,NLC
if (wavobs(ilc) .eq. -1.0) then
write(*,*) 'parcent tf 1:NLC', p(NPthcentidx:NPthcentidx+NLC-1)
write(*,*) 'stepcen tf 1:NLC', pscale(NPthcentidx:NPthcentidx+NLC-1)
write(*,*) 'parfwhm tf 1:NLC', p(NPthfwhmidx:NPthfwhmidx+NLC-1)
write(*,*) 'stepfwhm tf 1:NLC', pscale(NPthfwhmidx:NPthfwhmidx+NLC-1)
write(*,*)
exit
endif
enddo

if (iburnin(1) .eq.0) then
if (i_tfx == 1) then
write(*,*) 'T1v (k):', p(NPtridx+2), 'step:',pscale(NPtridx+2)
write(*,*) 'T1i (k):', p(NPtridx+3), 'step:',pscale(NPtridx+3)
else
write(*,*) 'mdot (Mo/yr) :', p(NPumbhidx), 'step:',pscale(NPumbhidx)
endif
write(*,*) 'inlcination:', acos(p(NPcosincidx))*180/pi, '   cos(inc) scale', pscale(NPcosincidx)
write(*,*) 'T(r) slope, visc, rad:',p(NPtridx),p(NPtridx+1),'    stepsize',pscale(NPtridx),pscale(NPtridx+1)
else !print the average and sd of the parameters after convegence

if (i_tfx == 1) then
write(*,*) 'T1v (k):', avparm(NPtridx+2), '+/-:',sdparm(NPtridx+2)
write(*,*) 'T1i (k):', avparm(NPtridx+3), '+/-:',sdparm(NPtridx+3)
else
write(*,*) 'mdot (Mo/yr) :', avparm(NPumbhidx), '+/-:',sdparm(NPumbhidx)
endif
write(*,*) 'inlcination:', av_inc, '+/-:', sd_inc
write(*,*) 'T(r) slope, visc, rad:',avparm(NPtridx),'+/-',sdparm(NPtridx),'  ',avparm(NPtridx+1),'+/-',sdparm(NPtridx+1)

endif

write(*,*) 'Acceptance frac',1.*account/NP,account
write(*,*) 'Driving light curve offset:', p(NPfurconstidx),'      stepsize:',pscale(NPfurconstidx)
write(*,*) 'stretch factors', p(NPscaleidx:NPscaleidx+NLC-1),'       stepsize:',pscale(NPscaleidx:NPscaleidx+NLC-1)
write(*,*) 'power spectrum p0,w0 (rads-1),alpha', p(NPpspecidx),p(NPpspecidx+1)/2/3.1415,p(NPpspecidx+2),'      stepsize:'&
,pscale(NPpspecidx:NPpspecidx+2)


write(*,*) 'P0prior information',p(NPpspecidx),p0mean,sigp0square



if (timeon) then
write(*,*) ''
write(*,*) 'time check of the nobreak rouine (find the bottleneck)'
write(*,*) 'main bit', itimemain
write(*,*) 'calculate bof', itimebof
write(*,*) 'save bofs', itimesavebof
write(*,*) 'plotting routine', itimeplot
write(*,*) 'calc error envs', itimesaveenv
write(*,*) 'saving parameters', itimesavebackup
write(*,*) ''
endif


! if burnin not yet reached plot parameters and stepsizes

if (version .eq. 1) then
write(*,*) 'AGN E(B-V):', p(NPebmvidx),'     scaling',pscale(NPebmvidx)
write(*,*) 'luminosity distance:',p(NPdlidx),'       stepsize:',pscale(NPdlidx)
write(*,*) 'galaxy flux', p(NPgalidx:NPgalidx+NLC-1), '      sstepsize:',pscale(NPgalidx:NPgalidx+NLC-1)
write(*,*) 'fdisk',fdisk
else if (version .eq. 3) then
write(*,*) 'offsets', p(NPoffsetidx:NPoffsetidx+NLC-1),'   stepsize:',pscale(NPoffsetidx:NPoffsetidx+NLC-1)
write(*,*) 'T(r) parameters, visc, rad:',p(NPtridx),p(NPtridx+1),'    stepsize',pscale(NPtridx),pscale(NPtridx+1)
endif
if (sigexpand) then
write(*,*) 'Error bar expansion factor:',p(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1),'  stepsize:'&
, pscale(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1)
endif
if (varexpand) then
write(*,*) 'variance expansion parameter:',p(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1),'  stepsize:'&
, pscale(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1)
endif

write(*,*) ''



write(*,*) ''
write(*,*) ''
write(*,*) ''







!add diagnostic files to log the proposal, bofprop, bofold and stepsizes throughout the simulation
!only use this to diagnose problems as will use lots of disk space
if (stepdiag) then
open(unit = 1, file = 'stepdiag_step.dat', access = 'append')
write(1,*) pscale(1:NP)
close(1)

open(unit = 1, file = 'stepdiag_pprop.dat', access = 'append')
write(1,*) pnew_save(1:NP)
close(1)

open(unit = 1, file = 'stepdiag_psave.dat', access = 'append')
write(1,*) pold_save(1:NP)
close(1)

open(unit = 1, file = 'stepdiag_bofprop.dat', access = 'append')
write(1,*) bofnew_save(1:NP)
close(1)

open(unit = 1, file = 'stepdiag_bofsave.dat', access = 'append')
write(1,*) bofold_save(1:NP)
close(1)
endif

!!! Once this magic spell is through, what was one shall now be two.
!! introduce additional step to stop the code whenever the badness of fit remains unchanged for more than 10 cycles


if (account .lt. 1 .and. iaffinenow .eq. 0) then
numreject=numreject+1
open(unit=1,file='CREAM_crash_noparmsaccepted.dat',status='replace')
write(1,*) iteration, bofold, NPcosincidx, NPsigexpandidx, NPscaleidx, NPoffsetidx, NPtridx,&
Npthcentidx,NPthfwhmidx,NPumbhidx,NPvarexpandidx
do ip = 1,NP
if (pversion(ip) .eqv. .true.) then
if (ip .eq. NPumbhidx) write(1,*) 'embh!'
if (ip .eq. NPvarexpandidx) write(1,*) 'VAR!'
if (ip .eq. Npoffsetidx) write(1,*) 'offset'
if (ip .eq. Npthcentidx) write(1,*) 'thcent'
if (ip .eq. Npthfwhmidx) write(1,*) 'thcfwhm'
write(1,*) ip,pnew_save(ip), pold_save(ip), pscale(ip), bofold_save(ip), bofnew_save(ip),&
varerln, ervar(1), altertf
!write(1,*) bofsave_new(1:5,ip)
write(1,*)
endif
enddo
close(1)
write(*,*) 'BLARGHHHHHH! iaffine now?',iaffinenow
else
numreject=0
endif




if (numreject .eq. 50) then
write(*,*) 'more than 50 cycles have been rejected consecutively... info below'
write(*,*) 'ip, p(ip), pscael(IP), logstep, bofreject(ip),bofold=',bofold
do ip=1,NP
write(*,*) ip, p(ip),pscale(ip), psafe(ip),bofrejectbef(ip), bofreject(ip),bofrejectnew(ip)
enddo
write(*,*) 'more than 50 cycles have been rejected consecutively... info above'

stop
endif
!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*) 'end of sub'


end subroutine




































!! Plot the posterior probability distribution using the routine scatterparmsmulti
subroutine postplot(iteration, filepath, ifloatinc, ifloatmmdot, ifloattrvisc)


character(1000) filepath
character(100),allocatable:: axlab_pp(:)
character(100) fname_lp
integer, allocatable:: idxparm_pp(:),idx_lp(:)
real, allocatable:: bof_pp(:), waste_pp(:), p_pp(:,:), preal_pp(:), pall_pp(:,:),&
ymin_lp(:),ymax_lp(:)
integer iteration_pp, ifloatinc, ifloatmmdot, ifloattrvisc
integer, allocatable:: isteplog_pp(:)

filepath = adjustl(filepath)


iburninshow_pp = 1
iburnin2_pp    = 1

Nparm_pp    = 0
maxparms_pp = 2*NLC+5
Nalong_pp   = maxparms_pp + 2

if (ifloatinc == 1) Nparm_pp = Nparm_pp + 1
if (ifloatmmdot == 1) Nparm_pp = Nparm_pp + 1
if (ifloattrvisc == 1) Nparm_pp = Nparm_pp + 1







Nits_pp=iteration - iburnin2_pp - 2 !numline(trim(filepath))-1-iburnin2_pp
Nits_pp = min((iteration - iburnin2_pp - 2), numline(trim(adjustl(filepath)))-2-iburnin2_pp)

if (Nits_pp .le. 0) goto 999 !end the subroutine if it cant find the data
!write(*,*) Nits_pp, iburnin2_pp, numline(trim(filepath)), iteration
!stop

!write(*,*) 'before allocatations'
allocate(axlab_pp(Nparm_pp),p_pp(Nits_pp, Nparm_pp), bof_pp(Nits_pp), waste_pp(Nalong_pp),&
pall_pp(Nits_pp,maxparms_pp),idxparm_pp(Nparm_pp),preal_pp(Nparm_pp),&
isteplog_pp(Nparm_pp), idx_lp(NParm_pp),ymin_lp(Nparm_pp),ymax_lp(Nparm_pp))
!write(*,*) 'after allocatations'
preal_pp(:)=0

do i = 1,Nparm_pp
if (ifloatmmdot == 1 .and. i .le. 1) then
idxparm_pp(i)=1
axlab_pp(i)=' M\u\b\.\d  (M\d\(2281)\uyr\u-1\d)'
axlab_pp(i)=adjustl(axlab_pp(i))
isteplog_pp(i)=1
idx_lp(i) = NPumbhidx
else if (ifloatinc == 1 .and. i .le. 2) then
idxparm_pp(i)=2
axlab_pp(i)='Cos(i)'
axlab_pp(i)=adjustl(axlab_pp(i))
isteplog_pp(i)=3
idx_lp(i) = NPcosincidx
else if (ifloattrvisc == 1 .and. i .le. 3) then
idxparm_pp(i) = 3
axlab_pp(i) = 'T(r) \ga'
axlab_pp(i)=adjustl(axlab_pp(i))
isteplog_pp(i) = 0
idx_lp(i) = NPtridx

else
write(*,*) 'There is a problem with postplot routine'
stop
endif
enddo







!! new 18th may embed the line plots into this routine for brevity

do ip = 1,Nparm_pp
ymin_lp(ip) = 666

ymax_lp(ip)= 555

enddo

fname_lp = 'cream_lines.PS'

if (iteration .gt. 1) then
call mypgline(parsave(1:NP,1:iteration-1),iteration-1,fname_lp,idx_lp,nparm_pp,np,axlab_pp,ymin_lp,ymax_lp)
endif



!! read the parameters
!write(*,*) 'reading in the parms for hist plot'
!write(*,*) numline(trim(filepath)),'help me!'
open(unit=1,file=trim(filepath))
!write(*,*) 'reading in the parms for hist plot 2'
!do i=1,iburnin2_pp
!read(1,*) waste_pp(1:Nalong_pp)
!write(*,*) iburnin2_pp, i
!enddo
!write(*,*) Nits_pp, i, 'arararaghhhhhhhhhhIhatescience',trim(filepath)
!write(*,*) trim(adjustl(dirworking))
do i=1,Nits_pp
!write(*,*) Nits_pp, i

read(1,*) iteration_pp, bof_pp(i), pall_pp(i,1:maxparms_pp)
!write(*,*) i,iteration_pp, bof_pp(i), pall_pp(i,maxparms_pp), 'sdas'
enddo
close(1)
!write(*,*) 'finished reading in paarms for hist plot'





do i=1,Nparm_pp
p_pp(1:Nits_pp,i)=pall_pp(1:Nits_pp,idxparm_pp(i))




if (idxparm_pp(i) .eq. 1) then  !! if we have mmdot /10^7 to plot in sensible units
p_pp(1:Nits_pp,i)=alog10(p_pp(1:Nits_pp,i)) !plot the log10 of the mmdot as the steps are uniform in log  !/1.e7
endif

enddo
!write(*,*) Nits_pp,idxparm_pp,Nparm_pp,ifloatinc,ifloatmmdot
!write(*,*) p_pp(1:Nits_pp,1)
!write(*,*) p_pp(1:Nits_pp,2)


!stop

!write(*,*) isteplog_pp, 'postplot routine'
!write(*,*) 'before scatter parms routine',axlab_pp,Nparm_pp
!call scatterparmsmulti(Nparm_pp,Nits_pp,p_pp,preal_pp,axlab_pp,bof_pp,iburninshow_pp,isteplog_pp,3)
!scatterparmsmulti(n,Nits,p,preal,axlab,bof,iburnin,isteplog,Ngroups)

!write(*,*) 'after scatter parms routine',axlab_pp,Nparm_pp
999 end subroutine


























!! calculate and plot the correa
!! new jun 12 2015, idxmost, Nmost capture the most correlated parameters to feed into the
!! python plotting routine

subroutine covcalc(pscale,pversion,iteration,Nmost)
integer pgopen,Nmost,idxmost(Nmost,2),idx_all(NP),idcortemp(NP)
real cormat(NP,NP),cormat_temp(NP,NP), greytr(6), greyref(3,3),pscale(NP),&
cortemp(NP)
real,allocatable:: rdxmost(:),rdxmost_new(:)
logical dotheplot /.True./, pversion(NP)
character(1000) dirtemp
!write(*,*) 'starting thing'
do ia = 1,NP
idx_all(ia) = ia
enddo

!! calculate the voariance matrix for the parsave array

call covpar(NP,iteration+1-it_BIG_atstart,parsave(1:NP,it_BIG_atstart:iteration),covmat)
call corpar(NP,covmat,cormat)
cormat_temp = cormat
do ip = 1,NP
if (pscale(ip) .eq. 0) then
cormat(ip,:) = 1.0
cormat(:,ip) = 1.0
endif
enddo






!! now add feature to turn off some parameters based on values in correlation matrix
! e.g we only want to step the top 10 correlated parameters to the inclination and mdot
! (disk parameters) to save time. The others dont affect the disk parameter estiamtes
! so can be left alone. NEEDS TESTING
if (skipcor) then
pversion(1:NPF) = .true.
cortemp(1:NP) = cormat(idref_skipcor,1:NP)
call sort1d(cortemp,NP,cortemp,idcortemp)
pversion(1:NPF) = .false.
idx = 0
icount = 0
do while (idx .lt. Nkeep_skipcor)
if (icount .ge. NP) exit
idcnow = idcortemp(NP-icount)
if (idcnow .ne. idref_skipcor .and. idcnow .le. NPF) then
pversion(idcnow) = .True.
idx = idx + 1
endif
icount = icount + 1
enddo

endif





!! Calculate the Nmost correlated parameters new june 12 2015
do im = 1,Nmost,2

call minmax2d(cormat_temp,NP,NP,idxlo1,idxlo2,idxhi1,idxhi2)

idxmost(im,1)   = idxlo1
idxmost(im,2)   = idxlo2
idxmost(im+1,1) = idxhi1
idxmost(im+1,2) = idxhi2
cormat_temp(idxlo1,idxlo2) = 0
cormat_temp(idxhi1,idxhi2) = 0
enddo

Nmost2 = 2*Nmost
allocate(rdxmost(Nmost2),rdxmost_new(Nmost2))
call unique(real(idxmost(1:Nmost,1)),Nmost,rdxmost(1:Nmost),nuniquehi)
call unique(real(idxmost(1:Nmost,2)),Nmost,rdxmost(Nmost+1:Nmost2),nuniquehi)
call unique(rdxmost,Nmost2,rdxmost_new,nuniquehi)


open(unit=1,file='pythoncov.dat',status='replace')
do ip = 1,NP
write(1,*) covmat(:,ip)
enddo
close(1)

open(unit=1,file='pythoncor.dat',status='replace')
do ip = 1,NP
write(1,*) cormat(:,ip)
enddo
close(1)

open(unit=1,file='pythonPARS.dat',status='replace')
do its = it_BIG_atstart,iteration
write(1,*) parsave(1:NP,its)
enddo
close(1)


open(unit=1,file='paridx.dat',status='replace')
write(1,*) 'NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx,&
NPoffsetidx,NPtridx,NPsigexpandidx'
write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx,&
NPoffsetidx,NPtridx,NPsigexpandidx
close(1)

open(unit=1,file='furfreq.dat',status='replace')
write(1,*) w(1:NW)
close(1)

!!! make correlation plots ...
!...... do just inc, mmdot, trlaw (corplot_important.pdf)
!...... do just 1st 10 fourier terms (corplot_fur1to10.pdf)
!...... do power spectrum and importantterms (corplot_pspec_important.pdf)
!...... do most correlated parms (corplot_most_cor)
call getcwd(dirtemp)
call chdir(trim(adjustl(dirmain)))
!open(unit=1, file='pythonhist.dat',status='replace')
!write(1,*) trim(adjustl(dirtemp))
!write(1,*) NPcosincidx,NPumbhidx,NPtridx,NPtridx+1
!write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx, NPoffsetidx,NPtridx
!write(1,*) 'corplot_important.pdf'
!write(1,*) wavem(1:NLC)
!close(1)
!call system('python pythonhist.py')


!open(unit=1, file='pythonhist.dat',status='replace')
!write(1,*) trim(adjustl(dirtemp))
!write(1,*) NPcosincidx,NPumbhidx, NPpspecidx,NPpspecidx+1,NPpspecidx+2,NPpspecidx+3
!write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx, NPoffsetidx,NPtridx
!write(1,*) 'corplot_pspec_important.pdf'
!write(1,*) wavem(1:NLC)
!close(1)
!call system('python pythonhist.py')
!
!
!open(unit=1, file='pythonhist.dat',status='replace')
!write(1,*) trim(adjustl(dirtemp))
!write(1,*) NPcosincidx,NPumbhidx,int(rdxmost_new(1:nuniquehi))
!write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx, NPoffsetidx,NPtridx
!write(1,*) 'corplot_most_cor.pdf'
!write(1,*) wavem(1:NLC)
!close(1)
!call system('python pythonhist.py')
!
!
!open(unit=1, file='pythonhist.dat',status='replace')
!write(1,*) trim(adjustl(dirtemp))
!write(1,*) NPtridx
!write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx, NPoffsetidx,NPtridx
!write(1,*) 'corplot_trvisc.pdf'
!write(1,*) wavem(1:NLC)
!close(1)
!call system('python pythonhist.py')


!open(unit=1, file='pythonhist.dat',status='replace')
!write(1,*) trim(adjustl(dirtemp))
!write(1,*) idx_all(1:100)
!write(1,*) NLC, NPF, NP, NPcosincidx, NPumbhidx, NPpspecidx, NPscaleidx, NPoffsetidx,NPtridx
!write(1,*) 'corplot_all.pdf'
!write(1,*) wavem(1:NLC)
!close(1)
!call system('python pythonhist.py')

call chdir(trim(adjustl(dirtemp)))
call system('rm pythonPARS.dat')
call system('rm pythoncor.dat')



!! NOT DONE YET NOW RUN THE PYTHON SCRIPT TO PLOT THE CORRELATIONS THEN DELETE THE FOLDER ABOVE !!
!pythonPARS! CHECK OUT THIS it_BIG_atstart thing. I think its wrong in line 3667 above!!.



!
!!! plot the covariance matrix as a greyscale image
!if (greyscale) then
if (dotheplot) then

NPhi = NP
do iplot = 1,2
if (iplot == 1) then
NPlo = 1
!ier = pgopen('corgrey_full.PS/CPS')!.PS/CPS')! pgopen('/Xserve')! help
else
NPlo = NPumbhidx
!ier = pgopen('corgrey_main.PS/CPS')
endif


call minmax2d(cormat(NPlo:NPhi,NPlo:NPhi),NPhi-NPlo+1,&
NPhi-NPlo+1,idxlo1,idxlo2,idxhi1,idxhi2)


cmin = cormat(idxlo1+NPlo-1,idxlo2+Nplo-1)
cmax = cormat(idxhi1+NPlo-1,idxhi2+Nplo-1)

!
!
!call pgbgw()
!call pgsvp(0.1,0.9,0.1,0.9)
!




!call pgswin(real(NPlo)-0.5,real(NPhi)+0.5,real(NPlo)-0.5,real(NPhi)+0.5)
!call pgbox( 'bcnsti', 0.0, 0, 'bcmstvi',0.0,0)
!
!call pgmtxt('TV',1.5,real(NPumbhidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'MM\u\b\.\d')
!call pgmtxt('LV',1.5,real(NPumbhidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'MM\u\b\.\d')
!
!call pgmtxt('TV',2.5,real(NPcosincidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'Inc')
!call pgmtxt('LV',2.5,real(NPcosincidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'Inc')
!
!call pgmtxt('TV',1.5,real(NPpspecidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'Pspec')
!call pgmtxt('LV',1.5,real(NPpspecidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'Pspec')

!call pgmtxt('TV',1.5,real(NPscaleidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'stretch')
!call pgmtxt('LV',1.5,real(NPscaleidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'stretch')
!
!call pgmtxt('TV',2.5,real(NPoffsetidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'offset')
!call pgmtxt('LV',2.5,real(NPoffsetidx-NPlo+0.5)/(NPhi-NPlo+1),0.5,'offset')
!

greytr(:) = 0
greytr(2) = 1
greytr(6) = 1

!call pggray( cormat,NP,NP,&
!NPlo,NPhi,NPlo,NPhi,cmin,cmax,greytr )


!! reference grid
!call pgsvp(0.95,0.99,0.9,0.99)
!call pgswin(0.5,3.5,0.5,3.5)
greyref(1:3,1) = -1.0
greyref(1:3,2) = 0.
greyref(1:3,3) = 1.0
cmin = -1.0
cmax = 1.0
!call pggray( greyref,3,3,&
!1,3,1,3,cmin,cmax,greytr )

!call pgend

enddo

endif
end subroutine
!!!!!!!!!!!!!!!!!!!!!! covsum is now iteration x the covariance matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





























end module













!!!!!!!!!!!!! Subroutines Outside of Module!!!!!!!!!!





















!!!!! Subroutine to predict the fourier terms, step sizes, and power spectrum p0 using IOS
!!! Assuming power law power spectrum with no break

!p(f) = p0 (w0/w)^2

!!! INPUTS   (NOTE: y = xray data if present, else it should be the lowest frequency light curve)
! Integers NP,Ndat,NW
! Real p(NP),pscale(NP),x(Ndat),y(Ndat),sig(Ndat)


!!! OUTPUTS
!real p(NP), pscale(NP), w(NW) fourier terms and step sizes after IOS
!     p0 power spectrum term in above definition
!	  w0 reference (Not break) frequency
subroutine pspec_predict(w,NW,p,pscale,NP,x,y,sig,Ndat,w0,p0,osout)

integer,intent(in):: NW, NP, Ndat
real,intent(in):: w(NW), x(Ndat), y(Ndat), sig(Ndat)
real,intent(out):: p(NP), pscale(NP), osout, w0, p0
real,allocatable:: xplot(:),yplot(:), sigp(:), ptest(:)
real rms, avg, med
logical optscale  /.True./, pestcheck /.False./
Nav=3
w0= avg(w(1:Nav),Nav)
p0=rms(y,Ndat)
xhi = x(Ndat)
xlo = x(1)



ip=1
do iw=1,NW
a=w0/w(iw)
pscale(ip:ip+1)=p0 * a*a
ip=ip+2
enddo



sum =0.d0
do i=1,Ndat
sum = sum + y(i)
enddo
ymean = sum /Ndat








!!! Now estimated step sizes, either start fourier terms from zero, or try to use IOS to give reasonable starting values
NWest=floor(0.5*NW)
allocate(sigp(2*NWest))
if (optscale) then

!call itoptscale(Ndat,NWest,x,y-ymean,sig,w(1:NWest),p(1:2*NWest))
p(1:NP) = 0
os = 0.0
call itoptscale2(Ndat,NWest,x,y-ymean,sig,dxplot,w(1:Nwest),p(1:2*NWest),sigp(1:2*NWest),os,sigos)
!write(*,*) os,'os',sigos,'sigos',ymean,'ymean'

!! implement a check to see if the estimates are reasonable
allocate(ptest(2*Nwest))
iw = 1
do ip = 1,2*Nwest,2
ptest(ip) = abs(p(ip))
!write(*,*) w(iw),p(ip),p(ip+1)
iw = iw + 1
enddo
!do it = 1,Ndat
!write(*,*) x(it),y(it),sig(it),ymean,'sdada'
!enddo
!read(*,*)

pmed = med(ptest,2*NWest)
do ip = 1,2*Nwest
!write(*,*)ip, ptest(ip), pmed


if (ptest(ip) .gt. 1000*pmed) then
write(*,*) 'Pspec_predict: fourier coeffs too big!'
write(*,*) 'Chopping off large terms...'
p(ip:2*Nwest) = 0
exit
endif
enddo


else
p(1:NP) =0.0
endif

p(2*Nwest+1:NP) = 0.0



!stop

osout = os + ymean



!!!!!!!!!!!!!!!!! below is all is a check to see if the estimates are reasonable
if (pestcheck) then
Nplot = 10*Ndat
allocate(xplot(Nplot),yplot(Nplot))


dxplot  = (xhi - xlo)/(Nplot - 1)



do iplot = 1,Nplot
xplot(iplot) = xlo +(iplot-1)*dxplot
yplot(iplot) = os+ymean
ip=1
do iw = 1,NWest
yplot(iplot) = yplot(iplot) + p(ip)*sin(w(iw)*xplot(iplot)) + p(ip+1)*cos(w(iw)*xplot(iplot))
ip=ip+2
enddo

enddo

ylo = min(minval(y),minval(yplot))
yhi = max(maxval(y),maxval(yplot))


!ier = pgbeg(0,'/Xserve',1,1)
!call pgenv(xlo,xplot(Nplot),ylo,yhi,0,1)

!call pgerrb(6,Ndat,x,y,sig,1.0)
!call pgline(Nplot,xplot,yplot)

stop
!call pgend


end if

end subroutine






!! 23/04/2015
!! subroutine to take in a prameter,, standard deviation, mean value and turn it into
! a string

!Input p1:  the current parameter value
!	p_ave:  the average of the parameter values (usually from some burnin to current iteratio)
!    p_sd:  The sd of '''
!      dp:  Tells code how many dp and whether to just use p1(-ve dp) or p_ave  +/- p_sd (+ve dp) (use -0.5 for using just p1 to 0dp)

! Output outstring (character*10): the parameter info in the appropriate format specified by dp

subroutine turnp2str(p1,p_ave,p_sd,dp,outstring)

character*10 otemp, otemp2
character*100 outstring
real p1,p_ave,p_sd, dp



if ((p_sd .ne. p_sd) .or. (1/p_sd == 0)) p_sd = -666.0
if ((p_ave .ne. p_ave) .or. (1/p_ave == 0)) p_ave = -666


if (dp .eq. -2.0) then
write(outstring,'(F10.2)') p1
else if (dp .eq. -1.0) then
write(outstring,'(F10.1)') p1
else if (dp .eq. -0.5) then
write(outstring, '(I10)') int(p1)
else if (dp .eq. 0.0) then
write(otemp2,'(I10)') int(p_sd)
write(otemp,'(I10)') int(p_ave)
write(outstring, '(I10)') trim(adjustl(otemp))//' \(2233) '//trim(adjustl(otemp))

else if (dp .eq. 1.0) then
write(otemp2,'(F10.1)') p_sd
write(otemp,'(F10.1)') p_ave
outstring = trim(adjustl(otemp))//' \(2233) '//trim(adjustl(otemp2))
else if (dp .eq. 2.0) then
write(otemp2,'(F10.2)') p_sd
write(otemp,'(F10.2)') p_ave

outstring = trim(adjustl(otemp))//' \(2233) '//trim(adjustl(otemp2))
else if (dp .eq. 3.0) then
write(otemp2,'(F10.3)') p_sd
write(otemp,'(F10.3)') p_ave
outstring =  trim(adjustl(otemp))//' \(2233) '//trim(adjustl(otemp2))
else !Produce well mannered error
write(*,*) 'turnp2str: You silly shit! I only take values for dp of -2.0 to 3.0 in integer steps'
write(*,*) 'oh and -0.5 if you have something that you only want to write out the current parameter'
write(*,*) 'value to 0dp (dp = 0.0) is for writing out average parameter and sd to 0dp)'
stop
endif



end subroutine











!turns lo res transfer function into high res transfer function


subroutine cheat_tf(Ntaucheat,&
avmmdot, sdmmdot, avcosinc, sdcosinc,wavinput,&
taucheat,tf_cheatlo,tf_cheatmid,tf_cheathi,p)


use creammod

integer Ntaucheat
real taucheat(Ntaucheat),tf_cheatlo(Ntaucheat),tf_cheatmid(Ntaucheat),tf_cheathi(Ntaucheat),&
p(NP)



rad2deg = 180./3.141593


do itau = 1,Ntaucheat
tf_cheatlo(itau) = 0.0
tf_cheatmid(itau) = 0.0
tf_cheathi(itau) = 0.0
enddo


embh1 = avmmdot - sdmmdot
cosinc1 = avcosinc - sdcosinc
deginc1 = acos(cosinc1)*rad2deg


!call tfb(taucheat,tf_cheatlo,Ntaucheat,wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,redshift,1)

if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
tf_cheatlo(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,p(NPtridx),p(NPtridx+1),&
embh1,deginc1,urin,tf_cheatlo,redshift)
endif
!write(*,*) 'lo',wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka
!write(*,*) taucheat(1:10)
!write(*,*) tf_cheatlo(1:10)
!read(*,*)



embh1 = avmmdot
cosinc1 = avcosinc
deginc1 = acos(cosinc1)*rad2deg
!call tfb(taucheat,tf_cheatmid,Ntaucheat,wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,redshift,1)

if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
tf_cheatmid(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,p(NPtridx),p(NPtridx+1),&
embh1,deginc1,urin,tf_cheatmid,redshift)
endif

!call tfb(taucheat,tf_cheathi,Ntaucheat,wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,1)

!write(*,*) 'mid',wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka
!write(*,*) taucheat(1:10)
!write(*,*) tf_cheatlo(1:10)
!read(*,*)

!
!
!
embh1 = avmmdot + sdmmdot
cosinc1 = avcosinc + sdcosinc
deginc1 = acos(cosinc1)*rad2deg
!call tfb(taucheat,tf_cheathi,Ntaucheat,wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka,redshift,1)

if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
tf_cheathi(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,p(NPtridx),p(NPtridx+1),&
embh1,deginc1,urin,tf_cheathi,redshift)
endif

!write(*,*) 'hi',wavinput,embh1,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,p(NPtridx),p(NPtridx+1),ur0,diskk,diskka
!write(*,*) taucheat(1:10)
!write(*,*) tf_cheathi(1:10)
!read(*,*)

psimax = max(maxval(tf_cheathi),maxval(tf_cheatlo))


do itau = 1,Ntaucheat
tf_cheathi(itau) = tf_cheathi(itau) / psimax
tf_cheatlo(itau) = tf_cheatlo(itau) / psimax
tf_cheatmid(itau) = tf_cheatmid(itau) / psimax
enddo



end subroutine


















!turns lo res transfer function into high res transfer function


subroutine cheat_tf_alpha(Ntaucheat,&
avmmdot, sdmmdot, avcosinc, sdcosinc, tralpha_av, tralpha_sd, wavinput,&
taucheat,tf_cheatlo,tf_cheatmid,tf_cheathi,p)


use creammod

integer Ntaucheat
real taucheat(Ntaucheat),tf_cheatlo(Ntaucheat),tf_cheatmid(Ntaucheat),tf_cheathi(Ntaucheat),&
p(NP), psilo(Ntaucheat), psihi(Ntaucheat), psimid(ntaucheat)

rad2deg = 180./3.14159265
do itau = 1,Ntaucheat
tf_cheatlo(itau) = 0.0
tf_cheatmid(itau) = 0.0
tf_cheathi(itau) = 0.0
enddo

embhref = 1.e7
ipar = 1


embh1   = avmmdot - sdmmdot
cosinc1 = avcosinc - sdcosinc
deginc1 = acos(cosinc1)*rad2deg


x = tralpha_av + tralpha_sd
xres = 0.9
dxres = x - xres
y = tralpha_av - tralpha_sd
z = tralpha_av

if (x .gt. xres) then
x = xres
y = y - dxres
z = z - dxres
endif

!call tfb(taucheat,psihi,Ntaucheat,wavinput,embhref,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,y,y,ur0,diskk,diskka,redshift,1)
if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
psihi(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,y,y,&
embhref,deginc1,urin,psihi,redshift)
endif

x = min(0.9,tralpha_av+tralpha_sd)
if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
psilo(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,x,x,&
embhref,deginc1,urin,psilo,redshift)
!call tfb(taucheat,psilo,Ntaucheat,wavinput,embhref,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,x,x,ur0,diskk,diskka,redshift,1)
endif

if (wavinput .gt. 0 .and.wavinput .lt. 100.0) then
psimid(1:Ntaucheat) = 0
else
call tfbx(taucheat,Ntaucheat,wavinput,&
-1*umdot,-1*umdot,z,z,&
embhref,deginc1,urin,psimid,redshift)
!call tfb(taucheat,psimid,Ntaucheat,wavinput,embhref,umdot,uhx,eff*(1.-alb),&
!   urin,deginc1,z,z,ur0,diskk,diskka,redshift,1)
endif

call tfupdate(Ntaucheat, taucheat, psilo, tf_cheatlo, embhref, avmmdot-sdmmdot, ipar)
call tfupdate(Ntaucheat, taucheat, psihi, tf_cheathi, embhref, avmmdot+sdmmdot, ipar)
call tfupdate(Ntaucheat, taucheat, psimid,tf_cheatmid,  embhref, avmmdot, ipar)


psimax = max(maxval(tf_cheathi),maxval(tf_cheatlo))


do itau = 1,Ntaucheat
tf_cheathi(itau) = tf_cheathi(itau) / psimax
tf_cheatlo(itau) = tf_cheatlo(itau) / psimax
tf_cheatmid(itau) = tf_cheatmid(itau) / psimax
enddo

end subroutine


































!! version 3 uses my own (wrong transfer code to calculate tr profile of disk assuming relatino of form
!  T(r)^4 = T0^4 (R0/R)^[4b] +Tx^4 (R0/R)^[4a] x(t-tau)
!! updat 3 aug. In parameter file creaminpar.par, only put in sig expansion and step sizes up to light curve 10. Cream then uses averages of these to calculate sig and step size for remaining light curves to make the .par file tidy
!..... same for galflux parameter too!

!! version 1: Uses my own (wrong) transfer code to calculate the delay functions

!!! mcmcnew uses logarithmic stepping of fourier parameters

!!! new fitting code mcmc
!!! superfit2b_fake includes the feature to generate a fake lightcurve and convolve this to give an echo lightcurve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NTPSI 882 in the old lcfitmcmc4 but hundreds of thousands here!!!!!
!!!!FIXXXXFIFXIXXFIIXFXFX

!!! new to version 3 used the convolution subroutine convolve.f90 instead of manually
!! performing the convolution
program creamfit

use creammod

implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Number of somethings N..
integer NAGN,Nignore,NLCtemp, ilctemp,idxlo,idxhi, ilc_poly_now, ipoly
!! Strings
character(1000) dirpath,agnname,outputname,xfile

!! logical arrays
logical,allocatable:: psafe(:),pversion(:)
logical firstrun/.true./, fake/.false./, dw_sub_0 /.false./, CREAM_th_ex

!! integer arrays
integer,allocatable:: Npoints(:),pAT(:),pRT(:)
!! integers
integer i,iw,it,tgridminidx,tgridmaxidx,iteration,nits,ip,ilc,ip2,itau
integer taucutoffidx,itaucent,itaumin,itaumax,account,iw2,iw3
!! integer functions
integer gtn,ier,pgbeg,iseed,NPF_2,ilcnow,ilclo,nnow,ilchi

!! real arrays
real,allocatable:: tgrid(:),xgrid(:),psigrid(:,:),taugrid(:)
real,allocatable:: t(:),x(:),er(:),xinterp(:),p(:),pscale(:), delcenttemp(:),delnormtemp(:)
real,allocatable:: woptscale(:),poptscale(:),xline(:)
!! real vals
real wlo,whi,maxtau,mintau,tgridlo,tgridhi,psicent,psinorm,psisig,psiexp,rmsxray,start0
real ptemp,fracalong,dtaugrid,delsigproptocent,delsigproptocentscale,delcenttempscale,delnormtempscale
real diskscale,delsig,siglogw0,offset,scale, avg,poly_new, power, pscalenow, rmsnow,&
tmaxtemp,srk

real stretch,fourierscale,ratio,option,emdotscale,cosincscale,thitemp,tmintemp,dt_med_temp
real bofnew,bofold,bof2temp,alphatemp,p0temp,w0temp,bof4temp,p0rms
real delsigproptocenttemp,delsigproptocentsig_square,delsigproptocentbof,dt,len1,z,uinc,sigp0,sigw0
real onst,alpha,umbh,inc,cosinc,alphascale,wavelength,stretch1,rms1,dummy, wavobs_temp,&
betascale,osout
integer idx,numline,starttime,endtime,NWest,numreject,num_p0_av,idircheck,iprob, Ndattemp, Ntemp
integer pgopen,idsigexp,idvarexp,inl,n_lim,idscale,idos

!! real functions
real medspac,rang,med,cisq,ran3,inttemp,rmstemp,rms,trash,bofnew1,bofnew2,bofnew3,bofnew4,dl0
real kascale,dl0scale,galscale,dltemp,trv,trr,offsetscale,ur0scale,pfirst,logp0root,logpfurstep
real p0,w0,f,w0steplog,p0steplog,p_ave,w_ave,w_sum,trscalea,trscaleb,pfurstep,p0root,a,&
dtmin,dtmintemp,sdev, sigmmdot, pri_sigmmdot_log, pri_mmdotmean_log, deginc_cut,&
dt_new,dt_temp,atemp,po,pn,tnow,amps,ampc,w0test,wnow,amp0, dwn, crapcent, crapfwhm,&
pricream_par, pricream_step,overide_os,overide_osstep,overide_st,overide_ststep

real,allocatable,dimension(:):: galflux,echo_sk,echo_ck,echo_sigsk,echo_sigck, test_x
character(10) ctest_no,strcatwhi,strcatnlc
character(100)cdate,coutput,strcat, fname
!! double precision (normally sums)
double precision xgridsum,xgridopsum,bof1,bof2,bof3,bof4,bof5,sumtau,xopsum,sum,&
sumk, s1, s2


version=3
iseed=2131412
call getcwd(dirmain)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! PART 0: Choose between real or fake data (load a subroutine drw_spec
filename2='fakelc'
!filename2= adjustl (filename2)
write(*,*)'Abandon hope all ye who enter here...'
write(*,*) ''
write(*,*) ''

dummy=1.*iteration


!!! PART 1: PREPARATION: read data, assign initial parameter values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call getcwd(cwd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the input parameters

inquire(file='creaminpar_main.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open (unit=1,file="creaminpar_main.par")  !! make sure there is one extra line at bottom of file
else
open (unit=1,file="creaminpar.par")
endif

read(1,'(A)') dirpath !! path to directory
dirpath = adjustl(dirpath)

close(1)


!! load the data
call chdir(trim(dirpath))
NLC = numline('creamnames.dat') -1
call getcwd(cwd)
write(*,*) trim(adjustl(cwd)), NLC



open(unit = 1,file = 'creamnames.dat')
read(1,*) fname, wavobs_temp ! check if the first line contains x ray data
close(1)
if ((wavobs_temp .eq. 0.0) .or. (wavobs_temp .eq. -666.0)) then
yesxray = .true.
if (noconvolve .eqv. .false. .and. wavobs_temp .eq. -666.0) positivity = .true.
NLC = NLC - 1
NXray = numline(trim(fname)) - 1
allocate(txray(NXray),xxray(Nxray), erxray(Nxray))
open(unit = 3, file = trim(fname))
do it = 1,Nxray
read(3,*) txray(it), xxray(it), erxray(it)
enddo
close(3)
else
positivity = .false.
yesxray = .false.
endif


!always turn yes xray off. Just assign a top hat response centered around 0 to the driver
yesxray =.false.

allocate(wavobs(NLC), lo(NLC), hi(NLC), Npoints(NLC+1))
Npoints(1) = 0

open(unit = 1,file = 'creamnames.dat')
if (yesxray) read(1,*)


do ilc = 1,NLC
read(1,*) fname, wavobs_temp
Ntemp = numline(trim(adjustl(fname))) - 1
Npoints(ilc+1)=Ntemp+Npoints(ilc)
if (ilc .gt. 1) then
lo(ilc) = hi(ilc-1) + 1
hi(ilc) = lo(ilc) + Ntemp - 1
else
lo(ilc) = 1
hi(ilc) = Ntemp - 1
endif
if (Ntemp .le. 0) then
write(*,*) 'something wrong with light curve: ',trim(adjustl(fname))
write(*,*) 'No points read. Check file name and contents.'
stop
endif
if (wavobs_temp .eq. 0) wavobs_temp = -1.0
wavobs(ilc) = wavobs_temp
enddo
close(1)

idx = 1

!! a cheat to turn off positivity!1 new 8/9/2015
write(*,*) NLC
if (wavobs(1) .eq. -666.0) positivity = .true.

allocate(t(Npoints(NLC+1)), x(Npoints(NLC+1)), er(Npoints(NLC+1)),&
ave_overflow(NLC),sd_overflow(NLC))

call getcwd(cwd)
write(*,*) trim(adjustl(cwd)), NLC
open(unit = 1,file = 'creamnames.dat')
if (yesxray) read(1,*)

allocate(file_lc_in(NLC))

do ilc = 1,NLC
read(1,*) fname, wavobs_temp
file_lc_in(ilc) = trim(adjustl(fname))
Ndattemp = numline(trim(adjustl(fname))) - 1
open(unit = 2,file = trim(adjustl(fname)))
write(*,*) 'reading',Ndattemp,'values from ...',trim(adjustl(fname)),&
lo(ilc), hi(ilc), Npoints(ilc+1)
do it = 1,Ndattemp
read(2,*) t(idx), x(idx), er(idx)
idx = idx+1
enddo
close(2)
enddo
close(1)



!! call data_load


!! load background parameters. If file does not exist then turn off variable background
write(*,*) 'To turn on background polynonmial fits...'
write(*,*) 'Make file "creaminpar_bg.par" in target light curve directory'
write(*,*) 'structured like...'
write(*,*) 'Nppoly - integer of order of polynomial'
write(*,*) ''
write(*,*) 'polyin(1:NPpoly) the coeffs of the polynomial fit, start at zero so background'
write(*,*) 'is just flat line equal to mean of light curve.'
write(*,*) ''
write(*,*) 'polyin_step(1:Nppoly,ilc) - stepsizes set to -1 to let code guess'
write(*,*) 'polyin_pri(1:Nppoly,ilc) - Gaussian priors on polynomial coefficients'
write(*,*) ''
inquire(file='creaminpar_bg.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
bgvary = .true.
open(unit=1,file='creaminpar_bg.par')
read(1,*) NPpoly
allocate(polyin(Nppoly),polyin_step(Nppoly,NLC),polyin_pri(NPpoly,NLC),&
polyin_pri_mean(NPpoly,NLC))
read(1,*) polyin(1:NPpoly)
do ilc = 1,NLC
read(1,*) polyin_step(1:Nppoly,ilc)
read(1,*) polyin_pri(1:NPpoly,ilc)
enddo
close(1)
polyin_pri_mean(:,:) = 0
else
Nppoly = 1
allocate(polyin(Nppoly),polyin_step(Nppoly,NLC),polyin_pri(NPpoly,NLC),&
polyin_pri_mean(NPpoly,NLC))
do ilc = 1,NLC
polyin_step(1:Nppoly,ilc) = 0
enddo
endif

!






call chdir(trim(adjustl(dirmain)))









allocate(delcenttemp(NLC)) !! temporarily store the delay pamrs
allocate(delnormtemp(NLC))
allocate(galflux(NLC))

if (yesxray) then          !! allocate the arrays that store the expansion factor
allocate(sigexpandparm(NLC+1),sigexpandsteplog(NLC+1))
else
allocate(sigexpandparm(NLC),sigexpandsteplog(NLC))
endif
NLCtemp = min(NLC,10)




open (unit=1,file="creaminpar.par")
read(1,*)
read(1,*)
read(1,*) pubplot, BIGsave
read(1,'(A)')
read(1,*) iplotsave
read(1,*) nbackup
read(1,*)
read(1,*) wlo, dw, start0

wlo = 2*pi*wlo
dw  = 2*pi*dw
if (dw .lt. 0) dw_sub_0 = .true.
read(1,'(I6)') NW
read(1,*) whi
whi=2*pi*whi
read(1,'(I6)') nits
read(1,'(I6)') AT
read(1,'(I6)') RT
read(1,*) taugridlo, maxtau
read(1,*) dtaugridtemp !! if this is -ve, dtaugrid will be determined automatically

read(1,*)
read(1,*)
read(1,*) fourierscale
read(1,*)
read(1,*)


read(1,*) sigexpand


if (yesxray) then
read(1,*) sigexpandparm(1:Nlc+1)
read(1,*) sigexpandsteplog(1:Nlc+1)
else
read(1,*) sigexpandparm(1:Nlc)
read(1,*) sigexpandsteplog(1:Nlc)
endif

write(*,*)


read(1,*)
read(1,*)
read(1,*) eff
read(1,*) alb
read(1,*) diskk
read(1,*) diskka
read(1,*) urin
urout=10000.0 !! for disksim subroutine
nr=1000
read(1,*) uhx
read(1,*) umdot
read(1,*) embhref
read(1,*) pri_mmdotmean_log, pri_sigmmdot_log
pri_mmdotmean = 10**(pri_mmdotmean_log)
if (pri_sigmmdot_log .lt. 0) then
sigmmdot = 0
else
sigmmdot = pri_sigmmdot_log*alog(10.) *10**pri_mmdotmean_log
endif
sigmmdotsquare = sigmmdot*sigmmdot



read(1,*) emdotscale
read(1,*)
read(1,*)
read(1,*) cosinc !! this number is initially inputted in degrees for the ease of the user
cosinc=cos(cosinc*pi/180)
read(1,*) cosincscale
read(1,*) deginc_cut, alpha_cosinc
if (deginc_cut .gt. 0) then
cosinc_cut = cos(deginc_cut*pi/180)
alpha_cosinc = abs(sin(deginc_cut*pi/180) * alpha_cosinc*pi/180)
sigcosinc2   = alpha_cosinc*alpha_cosinc
else
cosinc_cut = -2.0
endif

read(1,*)
read(1,*)

read(1,*) break
read(1,*) p0mean
read(1,*) siglogp0
!sigp0square=sigp0*sigp0
read(1,*) p0steplog
read(1,*) f0mean
w0mean = twopi*f0mean

read(1,*) siglogw0
read(1,*) w0steplog
read(1,*) alphamean
read(1,*) alphasig
read(1,*) alphascale
read(1,*) betamean
read(1,*) betasig
read(1,*) betascale

read(1,*)
read(1,*)
read(1,*) stretchscale
read(1,*) galflux(1:NLCtemp)
if (NLCtemp .lt. NLC) then  !!! if you have any light curves, dont clog up parameter file just use data from first 10 and get mean
s1 = 0.d0
do ilctemp = 1,NLCtemp
s1 = s1 + galflux(ilctemp)
enddo
do ilctemp = NLCtemp,NLC
galflux(ilctemp) = s1/NLCtemp
enddo
endif


read(1,*) galscale
read(1,*)
read(1,*)
read(1,*) redshift, omega_m, omega_l
redshiftadd1 = 1.+redshift
read(1,*) dl0
read(1,*) dl0scale
read(1,*) ebmvmw !Milky way extinction parameter
read(1,*) ebmvagn
read(1,*) ebmvagnsteplog
read(1,*)
read(1,*)
read(1,*) trv
read(1,*) trr
read(1,*) trscalea
read(1,*) trscaleb
read(1,*) ur0
read(1,*) ur0scale
write(*,*) ebmvmw,ebmvagn,w0mean, cosincscale
read(1,*) offsetscale
read(1,*)
read(1,*) noconvolve
if (noconvolve) then ! dont bother to step the mmdot or cosine parameter if we are not using echos
cosinc_cut = -2.
sigmmdotsquare = -2.
emdotscale = 0.0
cosincscale = 0.0
trscalea = 0.0
trscaleb = 0.0
endif

close(1)




!new for export version introduce creaminpar_main.par as a black box input file containing only the important plots
!if havent made the file then just use the creaminpar.par file. Use inquire statements to allow backwards compatibility
inquire(file='creaminpar_main.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit = 1, file = 'creaminpar_main.par')
read(1,*)
read(1,*) umdot, emdotscale
read(1,*) embhref
read(1,*) cosinc,cosincscale
read(1,*) redshift
cosinc=cos(cosinc*pi/180)
close(1)
endif

dirpath=adjustl(dirpath)


!! negative delay, Only allow for negtive delays if dealing with BLR lightcurves not continuum light curves
!taugridlo = 0.0
!do ilc = 1,NLC
!if (wavobs(ilc) .lt. 0) taugridlo = -1.*maxtau
!enddo










!! make fake generate fake data for each wavelength and store it in the ./fake/ folder
!if (trim(dirpath) .eq. './fake/') then
!fake = .true.
!call chdir(trim(dirpath),idircheck)

!if (idircheck .ne. 0) then
!call system('mkdir '//trim(dirpath)) !! if a fake directory does not exist, make it.
!call chdir(trim(dirpath),idircheck)
!endif

!if (fake) call system('cp ../fakelcpar.par ./') !! copy the parameters of the fake data to the output folder


!do ilc =1,NLC
!call fakelc(wavobs(ilc),yesxray)
!enddo


!call chdir('./..')
!call getcwd(cwd)

!endif
!! change to the input directory (this is ./fake/ if using fake data)

call chdir(trim(dirpath),idircheck)
if (idircheck .ne. 0) then
write(*,*) 'failed to change to directory,:',trim(dirpath)
call getcwd(cwd)
write(*,*) 'cwd:',trim(cwd)
stop
endif







!!new option to allow for custom response function lower and upper limits
write(*,*) 'create file cream_customlag.par in target directory to customize &
the lower and upper lag limits'
inquire(file='cream_customlag.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit = 1, file = 'cream_customlag.par')
read(1,*) taugridlo,maxtau
close(1)
endif



!! Load the X-ray data
if (yesxray) then

rmsxray=rms(xxray,Nxray)
xraymed=med(xxray,Nxray)
xxray(:)=xxray(:)!-xraymed   !!! xxray should be the variable component of the (need to subtract the median disk flux)

else
xraymed = 0.0
endif


!call superfitdata(t,x,er,Npoints,Nignore) !! Load the delay light curves



!!rest of program

!write(*,*) NLC

21 do i=1,NLC      !!! save the start and end indices of t,x,err for each light curve
lo(i)=Npoints(i)+1
hi(i)=Npoints(i+1)
enddo

Nt=Npoints(NLC+1)              !! combined number of data points
dt=(t(NT)-t(1))/(NT-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(bofsave(nits))







!assign overflow/underflow values
do ilc = 1,NLC
ave_overflow(ilc) = avg(x(lo(ilc):hi(ilc)),hi(ilc) - lo(ilc) + 1)
sd_overflow(ilc) = rms(x(lo(ilc):hi(ilc)),hi(ilc) - lo(ilc) + 1)

!boundaries = 1.e5
!if ((ave_overflow(ilc) .le. 1.e3) .and. (sd_overflow(ilc) .le. 1.e3)) then
!ave_overflow(ilc) = 0.0
!sd_overflow(ilc) = 1.0
!endif

x(lo(ilc):hi(ilc)) = (x(lo(ilc):hi(ilc)) - ave_overflow(ilc))/sd_overflow(ilc)
er(lo(ilc):hi(ilc)) = er(lo(ilc):hi(ilc))/sd_overflow(ilc)
write(*,*) 'light curve ',ilc,'   average', ave_overflow(ilc),' rms ',sd_overflow(ilc)
write(*,*) 'new ave sd',avg(x(lo(ilc):hi(ilc)),hi(ilc) - lo(ilc) + 1),&
rms(x(lo(ilc):hi(ilc)),hi(ilc) - lo(ilc) + 1)
enddo
!read(*,*)
!write(*,*) 'added code to normalise data to mean 0 and rms 1 if either mean or rms'
!write(*,*) 'of a light curve is > 1.e5. Stops overflow problems when'
!write(*,*) 'calcualting chi squared.'
!stop
!'Added feature in cream to detect when average or rms of light curve exceed 10^5.
!Light curves then renormalised to mean 0 and rms 1 to avoid overflow problems when
!recalculating chi squared and prior terms. Need to convert back into input units for
!output files. ave and rms data saved in arrays ave_overflow(1:NLC) and sd_overflow(1:NLC)'




!!!!! lower and upper frequency determination wlo,whi
!! define lower , upper and frequency spacings (go from 0 to Nyquist unless otherwise {stated see par file})
f=0.5


!if (dw .lt. 0) dw=f*2*pi/(maxval(t)-minval(t))

thitemp = maxval(t)
tmintemp = minval(t)
if (wlo .lt. 0) then
wlo = twopi / (thitemp - tmintemp) / 2
dw = wlo
write(*,*) 'automatic low frequency selection... (days), min max time, (days)', twopi/wlo,&
tmintemp,thitemp
!read(*,*)
endif


if (whi .lt. 0) then
dt_med_temp = medspac(t(lo(1):hi(1)),hi(1)-lo(1)+1) !! time resolution set by median sample spacing
whi = twopi/dt_med_temp
write(*,*) 'automatic high frequency selection... (days)', twopi/whi
!read(*,*)
!wlo=0.0!2*pi/300
endif

if (dw .lt. 0) then
dw = wlo
write(*,*) 'automatic delta frequency selection... (days)', twopi/dw
!read(*,*)
endif

!if (dw .eq. 0) then
!write(*,*) 'cant have 0 wlo and -ve dw term in creaminpar.par'
!stop
!endif



!whi=NW*dw+wlo
!!! option to automatically set NW if negative (goes up to Nyquist frequency by default)

nofur = 0

!! check for the minimum separation between adjacent points (this sets the high frequency)
dtmin = 10000.0
do ilc = 1,NLC
call minspac(hi(ilc) - lo(ilc) + 1,t(lo(ilc):hi(ilc)),dtmintemp)

if (dtmintemp .lt. dtmin) dtmin = dtmintemp
enddo




if (whi .lt. 0 .and. NW .lt. 0) then



whi=2*pi/dtmin

if (1/whi .eq. 0) whi = 2*pi/(t(hi(1))-t(lo(1)))*(hi(1)-lo(1)+1)
NW=(whi-wlo)/dw

else if (whi .gt. 0 .and. NW .lt. 0) then

NW=(whi-wlo)/dw


else if (NW == 0) then !here we domnt want to use the fourier parameters bbut cream requires you to assign a few of them to avoid a crash. We force them all to zero later on indicated by the nofur variable
NW = 5
nofur = 1
whi=2*pi/dtmin
endif


!!!!!!!!!! grids T, TAU AND PSI GRIDS  TGRID AND TAUGRID HAVE SAME SPACING
if (dtaugridtemp .eq. 666.0) then
dtgrid=pi/whi!medspac(t(lo(1):hi(1)),hi(1)-lo(1)+1) !! time resolution set by median sample spacing
else
dtgrid = abs(dtaugridtemp)
endif

if (nofur == 1) dtgrid = (thitemp - tmintemp)/10

!if (dtaugridtemp .gt. 0) dtgrid=dtaugridtemp !! if read -ve dtaugridtemp from file set dt grid automatically


Ntaugrid=ceiling(((maxtau-taugridlo)*redshiftadd1)/dtgrid) !new 29/04/2015



if (yesxray) then
tgridlo=min(minval(t),minval(txray))-(Ntaugrid+1)*dtgrid        !! mintau should be -ve -dt just to be sure
tgridhi=max(maxval(t),maxval(txray)) + dt - taugridlo              !! +dt just to make sure
else
tgridlo=minval(t)-(Ntaugrid+1)*dtgrid        !! mintau should be -ve -dt just to be sure
tgridhi=maxval(t) + dt - taugridlo
endif



tgridhiextra=(tgridhi-tgridlo)/20
tgridloextra=(tgridhi-tgridlo)/20
tgridhi=tgridhi+tgridhiextra!(tgridhi-tgridlo)/4
tgridlo=tgridlo-tgridloextra
Ntgrid=ceiling((tgridhi-tgridlo)/dtgrid)+1

write(*,*) 'Ntaugrid:',Ntaugrid,'     Ntgrid:',Ntgrid,'     dt:',dt
write(*,*) 'tgridhi:',tgridhi,'      tgridlo:',tgridlo,'     dtgrid:',dtgrid

allocate(tgrid(Ntgrid))          !! t grid
allocate(xgrid(Ntgrid),xgridold(NTgrid),xgridold_affine(NTgrid))          !! x grid
xgrid(1:Ntgrid) = 0
xgridold(1:Ntgrid) = 0
xgridold_affine(1:Ntgrid) = 0
!allocate(xgridop(NTgrid))        !! convolved x grid
allocate(xgridplot(Ntgrid,NLC))  !! plotted x grid (seesuperfitplot subroutine)
allocate(xgridplotsave(Ntgrid,NLC),xgridplotsave_affine(NTgrid,NLC)) !! saved x grid (for rejections)
allocate(taugrid(Ntaugrid))       !! taugrid goes from mintau to max tau in dtau steps
allocate(taugrididx(Ntaugrid))
allocate(psigrid(Ntaugrid,NLC))
allocate(psisave(Ntaugrid,NLC),psisave_affine(Ntaugrid,NLC))
allocate(xinterp(Nt))
allocate(thtempcent(NLC), thtempfwhm(NLC),thtempcent_scale(NLC), thtempfwhm_scale(NLC),&
thtempfwhm_pri_width(NLC),thtempcent_pri_width(NLC),thtempfwhm_pri(NLC),&
thtempcent_pri(NLC))
allocate(senv1(Ntgrid,NLC),senv2(Ntgrid,NLC),senv1tf(Ntaugrid,NLC),senv2tf(Ntaugrid,NLC),& !! the information to save the error envelopes in the plot
senv1_bg(Ntgrid,NLC), senv2_bg(Ntgrid,NLC),bg_save(Ntgrid,NLC)) !! the information to save the error envelopes in the plot

senv1    = 0
senv2    = 0
senv1_bg = 0
senv2_bg = 0
senv1tf  = 0
senv2tf  = 0
!write(*,*) Ntgrid,Ntaugrid,wlo,whi,NT,maxtau,dtgrid
!stop
dtgrid_em = dtgrid / (1.+redshift)

do it=1,Ntgrid
tgrid(it)=tgridlo+(it-1)*dtgrid!(tgridhi-tgridlo)/Ntgrid
enddo


!shall I wait for user to read instructions on terminal
inquire(file='cream_sleepoff.par',exist=CREAM_th_ex)
if (CREAM_th_ex) sleepon = .false.











!! find shortest spacing between adjacent points dont let code run if you haven't gone high enough
dt_new = t(NT) - t(1)
if (yesxray) then
do it = 2, Nxray
dt_temp = txray(it) - txray(it-1)
if (dt_temp .lt. dt_new) dt_new = dt_temp
enddo
endif

do ilc = 1,NLC
do it = lo(ilc)+1,hi(ilc)
dt_temp = t(it)-t(it-1)
if (dt_temp .lt. dt_new) dt_new = dt_temp
enddo
enddo

write(*,*) 'shortest spacing is', dt_new
atemp = twopi/whi
write(*,*) 'lowest, highest, and freq sep are (days)...',twopi/wlo,twopi/whi,twopi/dw
if (atemp .gt. dt_new) then
write(*,*) 'Select higher upper frequency or rebin data... (or press enter to continue)'
!read(*,*)
endif








!do ilc = 1,NLC
!
!a = minval(x(lo(ilc):hi(ilc)))
!
!if (a .lt. 0) then
!write(*,*) 'Light curve', ilc, 'has -ve values... adding 2*minval for no reason in particular'
!!read(*,*)
!do it = lo(ilc), hi(ilc)
!x(it) = x(it) - 2*a
!enddo
!endif
!
!enddo

!!!!!! end of -ve value problem (if even a problem)






!! GRID AND FREQUENCY INFO DEFINED SIN AND COS(W*TGRID) HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! freqs



allocate(w(NW),swgrid(NW,Ntgrid),cwgrid(NW,Ntgrid))
do iw=1,NW
w(iw)=wlo+dw*(iw-1)
enddo



!write(*,*) wlo,w(1)
!stop
do it=1,Ntgrid
do iw=1,NW
swgrid(iw,it)=sin(w(iw)*tgrid(it))
cwgrid(iw,it)=cos(w(iw)*tgrid(it))
enddo
enddo


!! make crappy test plot
!amp0 = 1.0
!w0test = w(1)
!allocate(test_x(Ntgrid))
!test_x(1:Ntgrid) = 0
!do iw = 1,NW
!wnow = w(iw)
!!!a = rang(0.,1.0,iseed)
!amps = rang(0.,amp0*w0test/wnow,iseed)
!ampc = rang(0.,amp0*w0test/wnow,iseed)
!do it = 1,Ntgrid
!tnow = tgrid(it)
!test_x(it) = test_x(it) + amps*sin(wnow*tnow) + ampc*cos(wnow*tnow)
!enddo
!enddo
!
!!ier =pgopen('/Xserve')
!!call pgsvp(0.1,0.9,0.1,0.9)
!!call pgswin(tgrid(1), tgrid(NTgrid), minval(test_x), maxval(test_x))
!!call pgbox( 'bcmst', 0., 10.0, 'bcnst', 0., 10. )
!!call pgline(Ntgrid,tgrid,test_x)
!
!!call pgend
!deallocate(test_x)
!write(*,*) w(1), w(NW),dw, wlo, dw_sub_0
!stop




do itau=1,Ntaugrid
taugrid(itau)= taugridlo + (itau-1)*dtgrid_em
taugrididx(itau)=taugrid(itau)/dtgrid_em      !! gives the index you need to subtract from xgrid during the convolution

!write(*,*) taugrid(itau),taugrididx(itau)

enddo

write(*,*) 'tau lims', taugrid(1), taugrid(Ntaugrid)


allocate(interpidx(Nt))                   !! find the tgrid idx first higher than t(it) for the interpolation in the convolution

if (yesxray) then
allocate(xrayinterpidx(Nxray))
do it=1,Nxray
xrayinterpidx(it)=gtn(tgrid,txray(it),1,Ntgrid)
enddo
endif

!write(*,*) 'gtn'
do it=1,NT
interpidx(it)=gtn(tgrid,t(it),1,Ntgrid)        !! interpidx stores the grid positions to the right of each of the data points
!write(*,*) 'interpolation,',t(it),interpidx(it)!tgrid(interpidx(it)),tgrid(interpidx(it)-1)
enddo


interpidxmin = Ntaugrid
!interpidxmax = Ntgrid - idxtaulo


!!do it=1,Ntgrid
!if (tgrid(it) .gt. tgridlo+maxtau) then
!interpidxmin=it
!exit
!endif
!enddo

!write(*,*) tgrid(1), tgridlo, maxtau, tgrid(NTgrid), minval(t),interpidxmin,tgrid(interpidxmin)
!read(*,*)

!interpidxmin=minval(interpidx)-1 !! maximum and minimum indicees of tgrid (for plotting) (declared in module)
interpidxmax=maxval(interpidx)


!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! deal with monte carlo parameter setup. Optimise starting fourier terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NPF=2*NW                       			   !! idx of final fourier parms
NPfurconstidx=NPF+1

NPdelay=2
NPumbhidx=NPfurconstidx+1						   !! distribution parameters
NPcosincidx=NPfurconstidx+2

NPpspec = 4
NPpspecidx=NPcosincidx+1                   !! first idx of pspec parms


NPgal=NLC
NPgalidx=NPpspecidx+NPpspec                !! galaxy light contribution in flux units

NPscale = NLC                       !! number of stretch and scale parms
NPscaleidx=NPgalidx+NPscale        !! starting index of NLC*stretches then NLC*offsets

NPoffset=NLC
NPoffsetidx=NPscaleidx+NPscale

NPdlidx=NPoffsetidx+NPoffset                !! dl parameter

NPtr=4
NPtridx=NPdlidx+1

NPxray=1
NPxrayoffsetidx=NPtridx+NPtr

NPsigexpandidx=NPxrayoffsetidx+1
if (yesxray) then
NPsigexpand=NLC+1
else
NPsigexpand=NLC
endif



NPebmvidx=NPsigexpandidx+NPsigexpand+1
NPebmv=1


NP=NPFurconstidx+NPdelay+NPpspec+NPgal+NPscale+NPoffset+1+NPtr+NPxray+NPsigexpand+NPebmv      !! Total number of parms

!extraparameters for top hat transfer function
NPth = NLC
NPthcentidx = NP + 1
NPthfwhmidx = NPthcentidx + NLC
NP = NP + 2*NPth


!extra parameters for var0
if (yesxray) then
NPvarexpand = NLC + 1
else
NPvarexpand = NLC
endif

NPvarexpandidx = NP + 1
NP = NP + NPvarexpand

NPpolyidx = NP+1
NPpolytot = NPpoly*NLC
NP = NP + NPpolytot

allocate(p(NP))                !! allocate relevant parameter arrays
allocate(pscale(NP))
allocate(Psafe(NP))
allocate(bofreject(NP))
allocate(bofrejectnew(NP))
allocate(bofrejectbef(NP))
allocate(pversion(NP))
allocate(covmat(NP,NP))
allocate(senv1p(NP),senv2p(NP), parsave(NP,Nits),tref_poly(NLC),parsave_affine_T(Nits,NP))
allocate(sigexpandpri_cent(NPsigexpand), sigexpandpri_fwhm(NPsigexpand))
allocate(pRT_save(NP))

pRT_save(1:NP) = 0


!if (NLCtemp .lt. NPsigexpand) then  !!! if you have any light curves, dont clog up parameter file just use data from first 10 and get mean
! s1 = 0.d0
! s2 = 0.d0
! do ilctemp = 1,Nlctemp
! s1 = s1 + sigexpandparm(ilctemp)
! s2 = s2 + sigexpandsteplog(ilctemp)
! enddo
! do ilctemp = Nlctemp+1,NPsigexpand
!  sigexpandparm(ilctemp) = s1/NLCtemp
!  sigexpandsteplog(ilctemp) = s2/NLCtemp
! enddo
!endif




!!!! check if the user wants to tie the offset parameters of some light curves together
!! e.g that they step at the same time and are always seperated by a fixed amount
nlc_shareos = 0
inquire(file='shareos.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
write(*,*) 'file shareos.par present in target directory'
write(*,*) 'this ties together the offset parameters of specific light curves'
write(*,*) 'enter the light curve numerb (order that it appears in creamnames.dat)'
write(*,*) 'then the number of the light curve you want to tie it to'
write(*,*) 'then the fixed offset difference'
write(*,*) 'shareos_diff(nlc_shareos),ilc_shareos_a(nlc_shareos),ilc_shareos_b(nlc)'
nlc_shareos = numline('shareos.par') - 1
allocate(ilc_shareos_a(nlc_shareos),&
ilc_shareos_b(nlc_shareos),&
shareos_diff(nlc_shareos))
open(unit = 1, file = 'shareos.par')
do it = 1,nlc_shareos
read(1,*) ilc_shareos_a(it), ilc_shareos_b(it), shareos_diff(it)
enddo
close(1)
endif




!!!!!! check if the user wants to include priors on the errorbar expansion factors
sigexpandpri_fwhm(1:NPsigexpand) = -1
sigexpandpri_cent(1:NPsigexpand) =  0
if (sigexpand) then
inquire(file='cream_se_prior.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit=1,file='cream_se_prior.par')
read(1,*) crapcent
read(1,*) crapfwhm
close(1)
if (crapcent .lt. 0) then
sigexpandpri_cent(1:NPsigexpand) = abs(crapcent)
sigexpandpri_fwhm(1:NPsigexpand) = crapfwhm
else
open(unit = 1,file = 'cream_se_prior.par')
read(1,*) sigexpandpri_cent(1:Npsigexpand)
read(1,*) sigexpandpri_fwhm(1:Npsigexpand)
close(1)
endif
endif

write(*,*)
write(*,*)'SIGEXPAND is on. To include priors on the error bar expansion factors,'
write(*,*)'create a file "cream_se_prior.par" in the target directory of the form...'
write(*,*)'pri_cent(1:NLC)'
write(*,*)'pri_fwhm(1:NLC)'
write(*,*) 'or...'
write(*,*) '-pricent'
write(*,*) 'prifwhm'
write(*,*) 'to apply the same prior to all the light curves'
write(*,*)
write(*,*) 'current centroids...'
write(*,*) sigexpandpri_cent(1:NPsigexpand)
write(*,*) 'current fwhm (of log of prior)...'
write(*,*) sigexpandpri_fwhm(1:NPsigexpand)
write(*,*)
!call sleep(3)
endif
!!!!!!



psafe(1:NP)=.False.  !!! by default do not step in logs
pscale(1:NPF)=fourierscale


!!! do not vary all parameters in all versions of code. Sort this here.
pversion(1:NP) =.True.


if (sigexpand) then
psafe(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1) = .true.
pscale(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1) = sigexpandsteplog
p(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1) = sigexpandparm
else
pversion(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1)= .false. !! turn off error bar expansions if sigexpand set to false
endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! apply polynomial background constraints !!!!!!!!!!!!!!!!!!!!!
do ilc = 1,NLC
idxlo = NPpolyidx+(ilc-1)*NPpoly
idxhi = idxlo + NPpoly-1
pscale(idxlo:idxhi) = polyin_step(1:NPpoly,ilc)
p(idxlo:idxhi) = polyin(1:NPpoly)
enddo

!code guesses
do ip = NPpolyidx, NPpolyidx + NPpolytot -1,1
ipoly = mod(ip-NPpolyidx,NPpoly)
power = 1.+ipoly
poly_new = p(ip)
ilc_poly_now = floor(real(ip - NPpolyidx)/NPpoly) + 1
ilclo = lo(ilc_poly_now)
ilchi = hi(ilc_poly_now)
NNow = ilchi - ilclo + 1
if (polyin_step(ipoly+1,ilc_poly_now) .lt. 0) then
tmaxtemp = maxval(t)
tref_poly(ilc_poly_now) = avg(t(ilclo:ilchi),NNow)
rmsnow = rms(x(ilclo:ilchi),NNow)
pscalenow = rmsnow/4**real(power) / (tmaxtemp-tref_poly(ilc_poly_now))**power
pscale(ip) = pscalenow
endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! AGN extinction
psafe(NPebmvidx)=.True.
p(NPebmvidx)=ebmvAGN
pscale(NPebmvidx)=ebmvagnsteplog



if (version .eq. 1) then
pversion(NPtridx:NPtridx+NPtr-2)= .False. !! turn of parameter adjustments for the transfer function scaling if we are investigating the luminosity distance
pversion(NPoffsetidx:NPoffsetidx+NPoffset-1)=.False. !! turn off offset parms if version 1 as we use the luminosity distance to determine these.
else if (version .eq. 3) then
pversion(NPdlidx) = .False. !! don't vary the DL in versino 3 of the code (just interested in tr profile)
endif
if (yesxray .eqv. .False.) pversion(NPxrayoffsetidx) = .False.
!!!!


p(NPdlidx)=dl0 					!!luminosity distance
pscale(NPdlidx)=dl0scale
psafe(NPdlidx)=.True.

p(NPtridx)=trv   !! viscous parameter
p(NPtridx+1)=trr   !! irradiation parameter
psafe(NPtridx)=.True.
psafe(NPtridx+1)=.True.
pscale(NPtridx)=trscalea
pscale(NPtridx+1)=trscaleb



!check if the alternate disc parameterisation file is present this will overwrite
! mdot and slope parameters entered into creaminpar.par
inquire(file='cream_tfx.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit = 1, file = 'cream_tfx.par')
read(1,*) T1vnorm, T1vnorm_scale
read(1,*) T1inorm, T1inorm_scale
read(1,*) Trv, trscalea
read(1,*) Trr, trscaleb
close(1)
write(*,*) 'cream_tfx.par found using alternate accretion disc parameterisation...'

i_tfx = 1
p(NPtridx+2)=T1vnorm
p(NPtridx+3)=T1inorm
psafe(NPtridx+2)=.True.
psafe(NPtridx+3)=.True.
pscale(NPtridx+2)=T1vnorm_scale
pscale(NPtridx+3)=T1inorm_scale
p(NPtridx)=trv   !! viscous parameter
p(NPtridx+1)=trr   !! irradiation parameter
psafe(NPtridx)=.True.
psafe(NPtridx+1)=.True.
pscale(NPtridx)=trscalea
pscale(NPtridx+1)=trscaleb
pscale(NPumbhidx) = 0
else
p(NPtridx+2) = 0
p(NPtridx+3) = 0
pscale(NPtridx+2) = 0
pscale(NPtridx+3) = 0
i_tfx = 0
endif


!!! insert mmdot constraints
p(NPumbhidx)=umdot
psafe(NPumbhidx)=.True.  !! step bh mass logarithmically
p(NPcosincidx)=cosinc
pscale(NPumbhidx)=emdotscale
pscale(NPcosincidx)=cosincscale


!							             !! store the pspec defaults
!

!estimate fourier power spectrum properties
p(NPpspecidx+2) = alphamean
p(NPpspecidx+3) = betamean
w0mean = w(2)
wlo    = w(1)
P0mean     = 2*wlo/w0mean**2
p(NPpspecidx)   = p0mean
p(NPpspecidx+1) = w0mean
p(1:NPF) = 0

idx = 1
do iw = 1,Nw
pscale(idx:idx+1) = P0mean*(w(iw)/w0mean)**(-2)
idx = idx + 2
enddo

!if ((yesxray .eqv. .false.) .and. (iprob == 0) .and. (start0 == 0) ) then
! p(1:NPF) = p(1:NPF) *sqrt( p(NPpspecidx)*dw / (p(1)**2 + p(2)**2))
!endif
sigp0=siglogp0*p(NPpspecidx)
sigp0square=sigp0*sigp0

sigw0=siglogw0*p(NPpspecidx+1)
sigw0square=sigw0*sigw0

write(*,*) 'Finished estimating fourier coefficients'







!write(*,*) p(NPpspecidx:NPpspecidx+2),p0mean
!read(*,*)
!! for use in badness of fit priors
sigalphasquare=alphasig*alphasig           !! for use in the alpha prior

!!
pscale(NPpspecidx)   = p0steplog
psafe(NPpspecidx)    = .True.     !! the norm and break of the pspecs are stepped in log
pscale(NPpspecidx+1) = w0steplog
psafe(NPpspecidx+1)  = .True.
pscale(NPpspecidx+2) = alphascale
pscale(NPpspecidx+3) = betascale





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(fvb(NLC),fvf(NLC)) !! allocate the arrays to store the maximum and minimum freuqency
do i=1,NLC
fvb(i)=maxval(x(lo(i):hi(i)))
fvf(i)=minval(x(lo(i):hi(i)))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   STARTING PARAMETER ESTIMATES for Stretch, offset and p0!!!!!!!!
!!!!!!!!

!!!!!!!!!!!!!!!! offsets
ilc=1
do i=NPoffsetidx,NPoffsetidx+NPoffset-1
offset=med(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)!-p(NPscaleidx+ilc-1)*xraymed
p(i)=offset
pscale(i)=rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)
write(*,*) 'offset', offset, 'median',med(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1),&
'rms',pscale(i)
ilc=ilc+1
enddo


!!!!!!!!!!!!!!!! Mean disk fluxes
allocate(fdisk(NLC),fdisksave(NLC),fdisksave_affine(NLC))
!! set the mean disk fluxes in mjy
!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ilc=1,NLC
!! store the disk flux parameters !! fluxes in mjy
uinc=acos(p(NPcosincidx))*180/pi !! inclination in degrees
call disksim(p(NPdlidx),0.0,0.0,0.0,p(NPumbhidx),umdot,urin,urout,nr,uinc,diskk,diskka,eff,alb,0.0,uhx,&
wavobs(ilc),fdisk(ilc),1)
!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo														!!

if (version .eq. 1) then
p(NPgalidx:NPgalidx+NLC-1)      = galflux(1:NLC)							    !!
psafe(NPgalidx:NPgalidx+NLC-1)  = .True.
pscale(NPgalidx:NPgalidx+NLC-1) = galscale
else
p(NPgalidx:NPgalidx+NLC-1)      = 0.0
psafe(NPgalidx:NPgalidx+NLC-1)  = .false.
pscale(NPgalidx:NPgalidx+NLC-1) = 0.0
p(NPebmvidx)        = 0.0
pscale(NPebmvidx)   = 0.0
pversion(NPebmvidx) = .false.
endif


allocate(pAT(NP))
allocate(PRT(NP))
pAT(:)=0
pRT(:)=0





!!!!!!! Xray offset term

if (yesxray) then
p(NPfurconstidx)=med(xxray,Nxray)
pscale(NPfurconstidx)=0.1!p(NPfurconstidx)!rms(xxray,Nxray)/sqrt(real(Nxray))
!else if (positivity) then
!p(NPfurconstidx)=1.0
!pscale(NPfurconstidx)=0.00001
else
p(NPfurconstidx)=0.0
pscale(NPfurconstidx)=0.0
endif



!!!!!!!!!!!!!!!!!!!!!!!!! Fourier coefficients
!!!! optimise starting parameters using iterated optimal scaling  if we have xray data, use this to estimate the fourier terms,
!!!! else make it the shortest wavelength light curve
!
!write(*,*) 'Estimating Fourier coefficients...'
!pversion(NPpspecidx+1) = .False. !! skip the break frequency
!
!
!
!NPF_2 = NPF/2
!allocate(echo_sk(NPF_2), echo_ck(NPF_2), echo_sigsk(NPF_2), echo_sigck(NPF_2))
!
!
!if (yesxray) then
!!call iosecho(hi(1) - lo(1) + 1, txray, xxray, erxray, 1000.0, 0.0, 0.0,&
!! NW, w, w0mean, p(NPpspecidx), echo_sk, echo_ck, echo_sigsk, echo_sigck)
! call iosecho(Nxray, txray, xxray, erxray,&
!  NW, w, w0mean, p(NPpspecidx),  echo_sk, echo_ck, echo_sigsk, echo_sigck)
!
! if (noconvolve) then
!  call iosecho(Nxray, txray, xxray-p(NPfurconstidx), erxray,&
!  NW, w, w0mean, p(NPpspecidx),  echo_sk, echo_ck, echo_sigsk, echo_sigck)
! endif
!
!
!!write(*,*) positivity
!!stop
!
! idx = 1
! !do iw = 1,NW
! ! write(*,*) iw,w(iw), echo_sk(iw), echo_ck(iw), w0mean, p(NPpspecidx)
! !enddo
!
!
! !do it = 1,Nxray
! !write(*,*) txray(it), xxray(it), erxray(it)
! !enddo
!
!else
!
!
! call iosecho(hi(1) - lo(1) + 1, t(lo(1):hi(1)), x(lo(1):hi(1))-p(NPoffsetidx), er(lo(1):hi(1)),&
!  NW, w, w0mean, p(NPpspecidx),  echo_sk, echo_ck, echo_sigsk, echo_sigck)
!
!
!endif
!
!
!!call iosecho(hi(1) - lo(1) + 1, t(lo(1):hi(1)), x(lo(1):hi(1))-p(NPoffsetidx), er(lo(1):hi(1)), 1000.0, 0.0, 0.0,&
!! NW, w, w0mean, p(NPpspecidx), echo_sk, echo_ck, echo_sigsk, echo_sigck)
!!endif
!
!idx = 1
!do ip = 1,NPF,2
!
! p(ip) = echo_sk(idx)
! pscale(ip) = echo_sigsk(idx)
!
! p(ip+1) = echo_ck(idx)
! pscale(ip+1) = echo_sigck(idx)
! idx = idx+1
!enddo
!deallocate(echo_sk,echo_ck,echo_sigsk,echo_sigck)
!
!
!
!
!
!
!!! error check for nan parameter estimates
!idx = 1
!iprob = 0
!do iw = 1,NW
!if (p(idx) .ne. p(idx)) then
!write(*,*) 'Something has gone wrong estimating the fourier coefficients freq, sin, cos...', w(iw), p(idx), p(idx+1)
!iprob = 1
!endif
!idx = idx+2
!enddo
!if (yesxray) p(NPpspecidx) = sdev(xxray,nxray)**2
!
!
!idx=1
!do iw = 1,NW
!pscale(idx:idx+1) = sqrt(p(NPpspecidx) * dw * (w0mean/w(iw))**p(NPpspecidx+2))
!idx = idx+2
!enddo
!
!
!
!if ((iprob ==1) .or. (start0 == 1)) then
!w0mean = w(2)
! sumk = 0.d0
! do iw = 1,NW
!  sumk = sumk + 1./w(iw)**2
! enddo
!p(nppspecidx) = 2/dw/w0mean**2/sumk
!write(*,*) 'pspec predict'
!call pspec_predict(w,NW,p(1:NPF),pscale(1:NPF),NPF,t(lo(1):hi(1)),x(lo(1):hi(1))&
!,er(lo(1):hi(1)),hi(1)-lo(1)+1,w0mean,p(NPpspecidx),osout)
!
!p(NPoffsetidx+1 - 1) = osout
!p(NPpspecidx+1)=w0mean
!p(1:NPF) = 0
!endif

if (positivity) p(1:NPF)=0.0    !! start fourier coefficients at zero if enforcing positivity
!ip=1
!do iw=1,NW
!write(*,*) w(iw), p(ip:ip+1), pscale(ip:ip+1)
!ip=ip+2
!enddo









!!!!!!! End of fourier coefficient estimation


!write(*,*) siglogp0,p(nppspecidx)
!stop


!!!!!!!!!!!!!!!! Stretch Factor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rms1=rms(x(lo(1):hi(1)),hi(1)-lo(1)+1)



!!if (yesxray) then
!stretch1=rms1/sqrt(p(NPpspecidx)/2)
!!else
!!stretch1=1.0   !!if no xray set stretch1 =1 as this is the light curve used to perform ios on the fourier parms
!
! !!!! Normalise driving light curve to rms of 1 !!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! new 22nd June normalise driving light curve t rms = 1
! !sum = 0.d0
! !do iw = 1,NW
! ! sum = sum + 1./w(iw)**2
! !enddo
! !
! !po = p(NPpspecidx)
! !pn =  2./dw/p(NPpspecidx+1)**2 / sum
! !p(NPpspecidx) = pn
! !!stretch1 = po/pn * stretch1
! !write(*,*) stretch1,rms1, p(NPpspecidx)
! !stop
! !write(*,*) po,pn,stretch1
! !stop
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!endif
!
!p(NPscaleidx)=stretch1
!pscale(NPscaleidx)=p(NPscaleidx)
!do ilc=2,NLC
!
!p(NPscaleidx+ilc-1)=rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)/rms1 * stretch1
!
!write(*,*) 'scaling',  p(NPscaleidx+ilc-1), 'rms-',rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)
!write(*,*)
!
!enddo
!pscale(NPscaleidx:NPscaleidx+NLC-1) = 0.2*p(NPscaleidx:NPscaleidx+NLC-1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! set the stretch values for the light curve

rms1 = sqrt(p(NPpspecidx)/w(1))*p(NPpspecidx+1) !computed from theoretical power spectrum
do ilc = 1,NLC
p(NPscaleidx+ilc-1)=rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)/rms1
p(NPoffsetidx+ilc-1) = med(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)
enddo
pscale(NPscaleidx:NPscaleidx+NLC-1) = abs(0.5*p(NPscaleidx:NPscaleidx+NLC-1))
pscale(NPoffsetidx:NPoffsetidx+NLC-1) = abs(0.5*p(NPoffsetidx:NPoffsetidx+NLC-1))

!!
write(*,*) 'pscale',p(NPscaleidx:NPscaleidx+NLC-1)
write(*,*) 'pscale scale',pscale(NPscaleidx:NPscaleidx+NLC-1)
write(*,*) 'rms1',rms1
do ilc = 1,NLC
write(*,*) 'std ',rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!! Special case for only one light curve and no xray !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! if we only have one light curve and no x ray data, we turn off the scaling of p0.
!! to prevent degeneracy (26/10/2014)
if (yesxray .eqv. .false. ) then
pversion(NPpspecidx) = .false.
endif



do i=1,Ntgrid
xgrid(i)=0.0
iw=1

do ip=1,NPF,2
xgrid(i)=xgrid(i)+swgrid(iw,i)*p(ip)+cwgrid(iw,i)*p(ip+1)
iw=iw+1
enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End of parameter estimates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



do ip=1,NP
if (pscale(ip) .eq. 0) pversion(ip) = .false.
enddo

pversion(npcosincidx) = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!! PART 2  INSIDE THE ITERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!! WHERE TO STORE THE OUTPUT OF THE SUBROUTINE
!! make an output directory inside the current folder (./fake/ if using fake data)
!! there might be more than one run of the code for each set of data
! in this case we use the following directory structure
! ./fake/output_date_testno/   where testno is the number of tests started on that date

call date_and_time(cdate)
cdate=adjustl(cdate)

i=1
idircheck=0


do while (idircheck .eq. 0)

write(*,*) idircheck
write(ctest_no,"(I10.3)") i
ctest_no=adjustl(ctest_no)
coutput='./output_'//trim(cdate)//'_'//trim(ctest_no)
coutput=adjustl(coutput)





call chdir(trim(coutput),idircheck)




if (idircheck .ne. 0) then


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!    Save a brief summary of the file name and the varying parameters
!!!!!!!!!!!!!!!
!!! new 14/8/2014 make a catalogue of the output directories. i.e
! output_20140814_001  NLC,xray,mmdot,cosinc,sigexpand,whi
strcat=trim(coutput)//','
write(strcatnlc,'(I10.1)') NLC
strcatnlc=adjustl(strcatnlc)
write(strcatwhi,'(f10.1)') whi/2/pi
strcatwhi=adjustl(strcatwhi)
strcat=trim(strcat)//', NLC='//trim(strcatnlc)
if (yesxray) strcat=trim(strcat)//', yesxray '
if (pscale(NPumbhidx) .gt. 0) strcat=trim(strcat)//', mmdot '
if (pscale(NPcosincidx) .gt. 0) strcat=trim(strcat)//', cosinc '
if (sigexpand) strcat = trim(strcat)//', sigexpand '


strcat=trim(strcat)//', whi ='//trim(strcatwhi)

strcat=adjustl(strcat)



open(unit=1,file='cream_log.txt',access='append')
write(1,*) trim(strcat)
close(1)
!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call system('mkdir '//trim(coutput))
call system('cp fakelcpar_mc.par ./'//trim(coutput))
call chdir(trim(coutput),idircheck)
exit
else
call chdir('./..')
endif

i=i+1

enddo



!! check if we should turn on or off the python plots. By default these are on but if you are running very fast simulations (e.g just a constant offset), the performance bottleneck can be the time taken to produce the python plots.
inquire(file='../cream_nopythonplot.par',exist=CREAM_th_ex)
if (CREAM_th_ex) pythonplot = .false.

!! check if we should show the warning (information messages about the inquire(file='../cream_nopythonplot.par',exist=CREAM_th_ex)
if (CREAM_th_ex) pythonplot = .false.


!! copy creaminpar file to current directory so we have a record of the input simulation parameters
call system('cp ../../creaminpar.par ./')

!! copy the renormalisation consntants into the  the current directory so we have a record to
!transformback later
open(unit = 1,file = 'cream_renorm_const.par')
do ilc = 1,NLC
write(1,*) ave_overflow(ilc),sd_overflow(ilc)
enddo
close(1)

call chdir('./plots',idircheck)   !!! check for a plots directory
if (idircheck .ne. 0) then
call system('mkdir ./plots')
else
call chdir('./..')  !! if a plots directory does exist then you have just sank into it. GET out
endif

call getcwd(cwd)
write(*,*)'output dir: ', trim(cwd)!,trim('./output_'//trim(cdate)//'_'//trim(ctest_no))



!call system('mkdir '//trim('./output_'//trim(cdate)//'_'//trim(ctest_no)))
!! copy the fake input parameter file to the output directory
if (fake) call system('cp ../fakelcpar.par ./')!//trim('./output_'//trim(cdate)//'_'//trim(ctest_no)))



!write(*,*) 'EMERGENCY START', uhx,umdot,urin,ur0



!write(*,*) sum / NPF, p(1), p(NPF), maxval(p(1:NPF)), minval(p(1:NPF))
!stop


if (sigexpand .eqv. .false.) then
p(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1) = 1.0
pscale(NPsigexpandidx:NPsigexpandidx+NPsigexpand-1) = 0.0
endif

do ip = 1,NPF
write(*,*) p(ip), pscale(ip)
enddo




!apply outlier rejection dont start until model > itmin_sigrej iterations in
!to give code time to fit something sensible (stops crashes)

itmin_sigrej = 100
inquire(file='../sigrej.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
sigrej = .true.
allocate(sigrej_k(NLC),sigrej_mark(Npoints(NLC+1)),sigrej_mode(NLC),&
sigrej_mark_save(Npoints(NLC+1)))
sigrej_mark(1:Npoints(NLC+1)) = 0
sigrej_mark_save(1:Npoints(NLC+1)) = 0

open(unit = 1, file = '../sigrej.par')
do it = 1, Nlc
read(1,*) sigrej_mode(it), srk

!sigrej_k(it)

!! new 30march 2018 RM589, 767 and 789
!this is a more robust way of deciding where to put the outlier rejection threshold k
!the formula below estimates the k abs((D-M)/sig) of the largest outlier
if (srk .lt. 0) then
nnow = hi(it) - lo(it) + 1
sigrej_k(it) = 2.5 * ( alog10(1.*nnow) / alog10(100.0) )**(2./3) !!
write(*,*) 'lc',it,' auto sigrej=',sigrej_k(it)
else
sigrej_k(it) = srk
endif


enddo
close(1)

open(unit=1,file='sigrej_op.par')
do it = 1,NLC
write(1,*) sigrej_mode(it),sigrej_k(it)
enddo
close(1)

endif




!new 18/12/2017 pricream altered to allow one to set the parameters and step sizes as well as priors
!apply additional gaussian priors to cream parameters
inquire(file='../pricream.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
pricream = .true.
npricream = numline('../pricream.par') - 1
allocate(pricream_mean(npricream), pricream_idx(npricream), pricream_sd(npricream),&
pricream_bofnew(npricream), pricream_bofold(npricream))
open(unit = 1, file = '../pricream.par')

idsigexp = 0
idvarexp = 0
idos = 0
idscale = 0
do i = 1,npricream
read(1,*) pricream_idxnow, pricream_par, pricream_step, pricream_mean(i), pricream_sd(i)

!since user may not know the exact parameter number of each parameter, include some code values
!e.g if pricream_idx_now = -1 then mdot, , -2 inclination, -3 trviscslope, -4 tr iradslope


if (pricream_idxnow .eq. -1) then
pricream_idx(i) = NPumbhidx
if (pricream_par .ne. -1) p(NPumbhidx) = pricream_par
if (pricream_step .ne. -1) pscale(NPumbhidx) = pricream_step

else if (pricream_idxnow .eq. -2) then
pricream_idx(i) = NPcosincidx
if (pricream_par .ne. -1) p(NPcosincidx)       = pricream_par
if (pricream_step .ne. -1) pscale(NPcosincidx)  = pricream_step

else if (pricream_idxnow .eq. -3) then
pricream_idx(i) = NPtridx
if (pricream_par .ne. -1) p(NPtridx) = pricream_par
if (pricream_step .ne. -1) pscale(NPtridx) = pricream_step

else if (pricream_idxnow .eq. -4) then
pricream_idx(i) = NPtridx + 1
if (pricream_par .ne. -1) p(NPtridx+ 1) = pricream_par
if (pricream_step .ne. -1) pscale(NPtridx+1) = pricream_step

else if (pricream_idxnow .eq. -5) then

pricream_idx(i) = NPsigexpandidx + idsigexp
if (pricream_par .ne. -1) p(NPsigexpandidx+idsigexp) = pricream_par
if (pricream_step .ne. -1) pscale(NPsigexpandidx+idsigexp) = pricream_step
idsigexp = idsigexp + 1

else if (pricream_idxnow .eq. -6) then

pricream_idx(i) = NPvarexpandidx + idvarexp
if (pricream_par .ne. -1) p(NPvarexpandidx+idvarexp) = pricream_par
if (pricream_step .ne. -1) pscale(NPvarexpandidx+idvarexp) = pricream_step
idvarexp = idvarexp + 1

else if (pricream_idxnow .eq. -7) then

pricream_idx(i) = NPscaleidx + idscale
if (pricream_par .ne. -1) p(NPscaleidx+idscale) = pricream_par
if (pricream_step .ne. -1) pscale(NPscaleidx+idscale) = pricream_step
idscale = idscale + 1

else if (pricream_idxnow .eq. -8) then

pricream_idx(i) = NPoffsetidx + idos
if (pricream_par .ne. -1) p(NPoffsetidx+idos) = pricream_par
if (pricream_step .ne. -1) pscale(NPoffsetidx+idos) = pricream_step
idos = idos + 1

else

pricream_idx(i) = pricream_idxnow
if (pricream_par .ne. -1) p(pricream_idxnow) = pricream_par
if (pricream_step .ne. -1) pscale(pricream_idxnow) = pricream_step

endif

enddo
close(1)


endif









! only use quick mode if continuum fitting
!do ilc = 1,Nlc
!if (wavobs(ilc) .lt. 0) then
!quick_conv_itp = .false.
!endif
!enddo

!load information on tophat priors in cream_th.par
inquire(file='../cream_th.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit = 1,file = '../cream_th.par')
do ilc = 1,NLC
read(1,*) thtempcent(ilc), thtempcent_scale(ilc), thtempcent_pri(ilc),&
thtempcent_pri_width(ilc),&
thtempfwhm(ilc), thtempfwhm_scale(ilc), thtempfwhm_pri(ilc), thtempfwhm_pri_width(ilc)
enddo
close(1)
else
thtempcent(1:NLC)       = 0.1
do ilc = 1,NLC
thtempcent_scale(ilc) = 1.0 + (ilc-1)*0.1!3*dtgrid + (ilc-1)*0.2*dtgrid
enddo
thtempcent_pri(1:NLC)   = 0
thtempcent_pri_width(1:NLC)   = -1.0
thtempfwhm(1:NLC)       = 2.0!3*dtgrid
thtempfwhm_scale(1:NLC)   = 0.0
thtempfwhm_pri(1:NLC)     = 0.0
thtempfwhm_pri_width(1:NLC)   = -1.0
endif

!minimum fwhm of box car
do i =1,NLC
if (thtempfwhm(i) .lt. 3*dtgrid) thtempfwhm(i) = 3*dtgrid
enddo


!give a warning message to user
do ilc = 1,NLC
if (wavobs(ilc) == -1.0) then
write(*,*)
write(*,*) 'YOU have turned on the TOP HAT TF FUNCTION'
write(*,*) '(by setting one or more of the wavelengths to -1 in "creamnames.dat"),'
write(*,*) 'cent             1:NLC -->',thtempcent(1:NLC)
write(*,*) 'cent step        1:NLC -->',thtempcent_scale(1:NLC)
write(*,*) 'cent prior       1:NLC --> ',thtempcent_pri(1:NLC)
write(*,*) 'cent prior width 1:NLC --> (set = -1.0 to turn off)',&
thtempcent_pri_width(1:NLC)
write(*,*) 'fwhm             1:NLC -->',thtempfwhm(1:NLC)
write(*,*) 'fwhm step        1:NLC -->',thtempfwhm_scale(1:NLC)
write(*,*) 'fwhm prior       1:NLC --> ',thtempfwhm_pri(1:NLC)
write(*,*) 'fwhm prior width 1:NLC --> (set = -1.0 to turn off)',&
thtempfwhm_pri_width(1:NLC)
write(*,*)
write(*,*) 'above are the stats. You can change them by making a file called "cream_th.par"'
write(*,*) 'and putting it in the target directory along with the "creamnames.dat" file'
write(*,*) 'in the format ...'
write(*,*) 'cent, centstep, centprior (cent), centprior (width) ... and same for fwhm'
write(*,*) '1 wavelength per line. you might need an extra space at bottom of file'
!write(*,*) 'fortran is bollocks. Im doing this on Sunday and its christmas.. how sad is that!'
write(*,*)
if (sleepon) call sleep(3)
exit
endif
enddo


!!if the first wavelength is a driver 0.0 then make it a top hast with zero
!! width and centered on zero
!thtempcent(1) = 0.0
!thtempcent_scale(1) = 0.0
!thtempfwhm(1) = 0.1
!thtempfwhm_scale(1) = 0.0

!!



!!




!! turn on or off top hat transfer function parameters
ilcnow = 1
do ip = NPthcentidx, NPthcentidx + NLC - 1
!write(*,*) 'setting shit',ilcnow,wavobs(ilcnow),thtempcent(ilcnow),thtempcent_scale(ilcnow)
if (wavobs(ilcnow) .eq. -1) then !then we turn on top hat tf
p(ip) = thtempcent(ilcnow)
pscale(ip) = thtempcent_scale(ilcnow)
pversion(ip) = .true.
psafe(ip) = .false.
p(ip+NLC) = thtempfwhm(ilcnow)
pscale(ip+NLC) = thtempfwhm_scale(ilcnow)
pversion(ip+NLC) = .true.
psafe(ip+NLC) = .true.
else
p(ip) = 0.0!thtempcent(ilcnow)
pscale(ip) = 0.0!thtempcent_scale(ilcnow)
pversion(ip) = .false.
p(ip+NLC) = thtempfwhm(ilcnow)
pscale(ip+NLC) = 0.0!thtempfwhm_scale(ilcnow)
pversion(ip+NLC) = .false.
endif
ilcnow = ilcnow + 1
enddo











write(*,*) 'start offsets'
write(*,*) p(NPoffsetidx:NPoffsetidx+NLC-1)
write(*,*) pscale(NPoffsetidx:NPoffsetidx+NLC-1)
write(*,*)
write(*,*) 'start stretches'
write(*,*) p(NPscaleidx:NPscaleidx+NLC)
write(*,*) pscale(NPscaleidx:NPscaleidx+NLC)
write(*,*)
write(*,*) 'start cents th'
write(*,*) p(NPthcentidx:NPthcentidx+NLC-1)
write(*,*) pscale(NPthcentidx:NPthcentidx+NLC-1)
write(*,*) pversion(NPthcentidx:NPthcentidx+NLC-1)
write(*,*) psafe(NPthcentidx:NPthcentidx+NLC-1)
write(*,*)
write(*,*) 'start fwhm th'
write(*,*) p(NPthcentidx+NLC:NPthcentidx+2*NLC-1)
write(*,*) pscale(NPthcentidx+NLC:NPthcentidx+2*NLC-1)
write(*,*) pversion(NPthcentidx+NLC:NPthcentidx+2*NLC-1)
write(*,*) psafe(NPthcentidx+NLC:NPthcentidx+2*NLC-1)
write(*,*)



!! varexpand parameters
psafe(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) = .True.
write(*,*)
write(*,*) 'INFO FOR TURNING ON EXTRA QUADRATURE VARIANCE PARAMETER...'
write(*,*) 'put file "cream_var.par" in target directory to include variance scaling parameter'
write(*,*) 'syntax line1: varexpandparm(1:NPvarexpand)'
write(*,*) '       line2: varexpandsteplog(1:NPvarexpand) (set = -1 to just put in one number'
write(*,*) 'for each of these and apply it to all light curves or -2 to use average error'
write(*,*) 'bar size for each light curve)'
write(*,*)
inquire(file='../cream_var.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
varexpand = .true.
pversion(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) =.true.
allocate(varexpandparm(NPvarexpand),varexpandsteplog(NPvarexpand))
open(unit = 1,file = '../cream_var.par')
read(1,*) varexpandparm(1)
close(1)
if (varexpandparm(1) .eq. -1) then

if (yesxray) then
varexpandparm(1) = rms(xxray,Nxray)**2/50
varexpandsteplog(1) = 0.01
else
do ilc = 1,NLC
varexpandparm(ilc) = rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2/50
varexpandsteplog(ilc) = 0.01
enddo
endif

else
open(unit = 1,file = '../cream_var.par')
do ilc = 1,NLC
read(1,*) varexpandparm(ilc), varexpandsteplog(ilc)
if (varexpandparm(ilc) .gt. 0) then
varexpandparm(ilc) = varexpandparm(ilc) * rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2
else
varexpandparm(ilc) = abs(varexpandparm(ilc))
endif

if (varexpandsteplog(ilc) .eq. 0) varexpandsteplog(ilc) =0
enddo
close(1)
endif

! if (yesxray) then
!  varexpandparm(1) = avg(erxray,Nxray)**2 / 50
!  do ilc = 1,NLC
!  varexpandparm(ilc+1) = rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2/50
!  enddo
! else
!  do ilc = 1,NLC
!  varexpandparm(ilc) = rms(x(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2/50
!  enddo
! endif
! varexpandsteplog(1:NPvarexpand) =1.e-7
!else if (varexpandsteplog(1) .eq. -1) then !so dont have to fill file with numbers, can just put 1 in and set all light curves to that value
! varexpandparm(1:NPvarexpand) = varexpandparm(1)
! varexpandsteplog(1:NPvarexpand) = abs(varexpandsteplog(1))
!
!!if 3 then use the default starting values for variance expansion parametrers but set your own stepsizes
!else if (varexpandparm(1) .eq. -3) then !so dont have to fill file with numbers, can just put 1 in and set all light curves to that value
! if (yesxray) then
!  varexpandparm(1) = avg(erxray,Nxray)**2 / 50
!  do ilc = 1,NLC
!  varexpandparm(ilc+1) = avg(er(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2/50
!  enddo
! else
!  do ilc = 1,NLC
!  varexpandparm(ilc) = avg(er(lo(ilc):hi(ilc)),hi(ilc)-lo(ilc)+1)**2/50
!  enddo
! endif
! open(unit = 1,file = '../cream_var.par')
! read(1,*) varexpandsteplog(1)
! read(1,*) varexpandsteplog(1:NPvarexpand)
! close(1)
!
!
! !pversion(NPvarexpandidx:NPvarexpandidx) = .true.
!
!
!else
! open(unit = 1,file = '../cream_var.par')
! read(1,*) varexpandparm(1:NPvarexpand)
! read(1,*) varexpandsteplog(1:NPvarexpand)
! close(1)
!endif


idx = 1
do ilc = NPvarexpandidx, NPvarexpandidx+NPvarexpand-1
if ((varexpandsteplog(idx) == 0) .and. (varexpandparm(idx) == 0)) then
varexpandparm(idx) = 0.0
pversion(ilc) = .false.
endif
idx = idx + 1
enddo

p(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) = varexpandparm(1:NPvarexpand)
pscale(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) = varexpandsteplog(1:NPvarexpand)
write(*,*) 'variance expansion ON'
write(*,*) 'variance expansion parms...',p(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1)
write(*,*) 'variance expansion scale...',&
pscale(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1)
else
varexpand = .false.
pversion(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) =.false.
pscale(NPvarexpandidx:NPvarexpandidx+NPvarexpand-1) = 0
write(*,*) 'variance expansion OFF'


endif


!call sleep(4)




pscale(NPfurconstidx) = pscale(NPfurconstidx)/20.
pscale(1:NPF) = pscale(1:NPF)/20
if (yesxray) then
p(NPpspecidx) = p(NPpspecidx)/10
pscale(NPpspecidx) = p(NPpspecidx)/2
psafe(NPpspecidx) = .true.
pversion(NPpspecidx) = .true.
endif
!pscale(NPoffsetidx:NPoffsetidx+NLC-1) = pscale(NPoffsetidx:NPoffsetidx+NLC-1)/10.
!p(NPscaleidx:NPscaleidx+NLC-1) = p(NPscaleidx:NPscaleidx+NLC-1)/10.
!pscale(NPscaleidx:NPscaleidx+NLC-1) = pscale(NPscaleidx:NPscaleidx+NLC-1)/10.



!! turn off the fourier parameters if nofur
if (nofur == 1) then

pversion(NPumbhidx) = .false.
pscale(NPumbhidx) = 0
pversion(NPcosincidx) = .false.
pscale(NPcosincidx) = 0
do ip = 1,NPF
pversion(ip) = .false.
p(ip) = 0
pscale(ip) = 0
enddo

do ip = NPscaleidx,NPscaleidx+NPscale-1
pversion(ip) = .false.
p(ip) = 0
pscale(ip) = 0
enddo

endif


!! new feature 18 sept 2017 check if this is a continuation run or a fresh run
inquire(file='../../cream_resume.dat',exist=CREAM_th_ex)
if (CREAM_th_ex) then
write(*,*) 'This is a continuation run.. loading cream_resume.dat'
if (sleepon) call sleep(2)
open(unit = 1, file = '../../cream_resume.dat')
do iw = 1,NP
read(1,*) p(iw), pscale(iw), psafe(iw), pversion(iw)
write(*,*) iw, NP, p(iw), pscale(iw), psafe(iw), pversion(iw)
enddo
close(1)
call system('rm ../../cream_resume.dat')
endif




!new feature to force parameters to take certain values and stay there
!apply outlier rejection
inquire(file='../parfix.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
nparfix = numline('../parfix.par') - 1
open(unit = 1, file = '../parfix.par')
iosfix = 0
do it = 1, nparfix
read(1,*) ipnow, parfixnow
if (ipnow .eq. -1) then
if (parfixnow .ne. -1) then
pscale(NPoffsetidx+iosfix) = 0
p(NPoffsetidx+iosfix) = parfixnow
endif
iosfix = iosfix + 1
endif
enddo
close(1)

endif



!new feature to force parameters to take certain values and stay there
!apply outlier rejection
!skipcor,idref_skipcor,Nkeep_skipcor

inquire(file='../skipcor.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
skipcor = .true.
nparfix = numline('../skipcor.par') - 1
open(unit = 1, file = '../skipcor.par')
read(1,*) idref_skipcor, Nkeep_skipcor
close(1)

if (idref_skipcor == -1) then
idref_skipcor = NPumbhidx
else if (idref_skipcor == -1) then
idref_skipcor = NPcosincidx
endif

endif





! check for chaos file if it exists then allow all steps for transfer function parameters
inquire(file='../cream_chaos.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
chaos = .true.
endif



! check for simulated annealing parameter. If file not present then keep old acceptance probability.
! acceptance probability exp(-|dBOF|/(2*T_ann))
inquire(file='../cream_anneal.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit = 1, file = '../cream_anneal.par')
read(1,*) T_ann,frac_ann
close(1)
else
T_ann = 1.0
frac_ann = 1.0
endif





!!! creat a file to specify lower and upper limits of parameters dont use if set to -1
inquire(file='../cream_parlim.par',exist=CREAM_th_ex)
allocate(plo(NP),phi(NP))
plo(1:NP) = -1
phi(1:NP) = -1
if (CREAM_th_ex) then
n_lim = numline('../cream_parlim.par') - 1
do inl = 1,n_lim
read(1,*) idlim,plonow,phinow

if (idlim .eq. - 1) then
plo(NPumbhidx) = plonow
phi(NPumbhidx) = phinow
else if (idlim .eq. -2) then
plo(NPcosincidx) = plonow
phi(NPcosincidx) = phinow
else if (idlim .eq. -3) then
plo(NPtridx) = plonow
phi(NPtridx) = phinow
else if (idlim .eq. -4) then
plo(NPtridx + 1) = plonow
phi(NPtridx + 1) = phinow
else if (idlim .eq. -5) then
plo(NPtridx + 2) = plonow
phi(NPtridx + 2) = phinow
else if (idlim .eq. -6) then
plo(NPtridx + 3) = plonow
phi(NPtridx + 3) = phinow
else if (idlim .gt. 0) then
plo(idlim) = plonow
phi(idlim) = phinow
endif

enddo
endif





!! if stepping trvisc and irad slopes together,
!must not put both parameters in cream_affine.par, (just specify the viscous one)
!other wise the code will scream if you put a parameter matrix into jacobi that
!has two identical columns

!! allocate affine stepping array
iaffinenow = 0
inquire(file='../cream_affine.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
affine = .true.

n_affine = numline('../cream_affine.par') - 2
allocate(ip_affine(n_affine),pversion_affine(NP),p_affine(N_affine),&
covmat_affine(N_affine,N_affine))
pversion_affine(1:NP) = pversion(1:NP)
open(unit = 1, file = '../cream_affine.par')
idaff_sigexp = 0
idaff_varexp = 0
do it = 1,n_affine
read(1,*) idaffine
if (idaffine .eq. -1) then
ip_affine(it) = NPumbhidx
else if (idaffine .eq. -2) then
ip_affine(it) = NPcosincidx
else if (idaffine .eq. -3) then
ip_affine(it) = NPtridx
else if (idaffine .eq. -4) then
ip_affine(it) = NPtridx + 1
else if (idaffine .eq. -5) then
ip_affine(it) = NPsigexpandidx + idaff_sigexp
idaff_sigexp = idaff_sigexp + 1
else if (idaffine .eq. -6) then
ip_affine(it) = NPvarexpandidx + idaff_varexp
idaff_varexp = idaff_varexp + 1
else
ip_affine(it) = idaffine
endif
enddo
read(1,*) i_affinestart
close(1)

!introduce a parameter frac_stepchange_affine to control if we alter the step sizes in affine mode
if (i_affinestart .lt. 0) then
frac_stepchange_affine = 2.0
else
frac_stepchange_affine = 1.0
endif
i_affinestart = abs(i_affinestart)


endif











!copy the cream altering files to the current directory so we can resume later if needed
!in cream_resume.dat
call system('cp ../sigrej.par ./')
call system('cp ../pricream.par ./')
call system('cp ../cream_th.par ./')
call system('cp ../skipcor.par ./')
call system('cp ../cream_parlim.par ./')
call system('cp ../cream_affine.par ./')



!overide the standard step sizes for the background and offset parameters
inquire(file='../offsetstretch_fix.par',exist=CREAM_th_ex)
if (CREAM_th_ex) then
open(unit=1,file='../offsetstretch_fix.par')
do ilc=1,NLC
read(1,*) overide_os,overide_osstep,overide_st,overide_ststep
if (overide_os .ne. -1.0) then
p(NPoffsetidx+ilc-1) = overide_os
pscale(NPoffsetidx+ilc-1) = overide_osstep
end if
if (overide_st .ne. -1.0) then
p(NPscaleidx +ilc -1) = overide_st/rms1
pscale(NPscaleidx+ilc-1) = overide_st/rms1 * overide_ststep
end if
end do
end if



do iteration=1,nits

!if (iteration .gt. 100) pversion(NPF/2:NPF) = .true.

call getcwd(cwd)
cwd=adjustl(cwd)



!all steps to take if affine mode turned on step
!affine parameters every even iteration, step others every odd iteration

if (affine) then

if ((mod(iteration,2) .eq. 0) .and. (iteration .gt. i_affinestart)) then
iaffinenow = 1
pversion(1:NP) = .false.
do it = 1,n_affine
ipnow = ip_affine(it)
pversion(ipnow) = .true.
enddo
else
iaffinenow = 0
pversion(1:NP) = pversion_affine(1:NP)
do it = 1,n_affine
ipnow = ip_affine(it)
if (iteration .gt. i_affinestart) pversion(ipnow) = .false.
enddo
endif


endif




do ip = 1,N_affine
idaffine = ip_affine(ip)
do ip2 = 1,N_affine
idaffine2 = ip_affine(ip2)
covmat_affine(ip,ip2) = covmat(idaffine,idaffine2)
enddo
do it = 1,iteration
parsave_affine_T(it,ip) = parsave(idaffine,it)
enddo
enddo



call mcmcmulti_iteration(bofold,iteration,nits,x,t,er,taugrid,psigrid,&
xgrid,tgrid,Npoints,p,psafe,pscale,pversion,pAT,pRT)





enddo
if (nofur == 1) then
call system('cp ../../cream_python/python_clean.py ./')
call system('python python_clean.py')
call system('rm python_clean.py')
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program




!program convolvetest

!integer Nt
!real t(100),x(100),xop(100),psi(1000),taugrid(1000)

!Nt=100
!Npsi=1000

!!x
!x(1:Nt)=0
!x(10:11)=0.5

!!! t grid
!mint=0
!maxt=10
!dt=1.*(maxt-mint)/Nt
!do i=1,Nt
!t(i)=1.*(i-1)*dt
!x(i)=sin(t(i))
!write(*,*) t(i),dt

!enddo

!stop
!!! psi grid
!psi(1:Npsi)=0
!do i=1,Npsi
!taugrid(i)=(i-1)*dt
!enddo
!psi(10:11)=0.5


!! convolve
!do i=1,10
!call convolvelc(Nt,Npsi,t,x,taugrid,psi,xop)
!enddo
!write(*,*) 't',maxval(t)
!write(*,*) t(1:10)
!write(*,*) 'x',maxval(x)
!write(*,*) x(1:10)
!write(*,*) 'tau',maxval(taugrid),maxval(xop)
!write(*,*) taugrid(1:10)
!write(*,*) xop(1:10)




!!! plot the result

!!ier = pgbeg(0,'/XSERVE',1,3)
!!call pgsch(2.0)
!!call pgenv(t(1),t(Nt),minval(x),maxval(x),0,1)
!!call pgline(Nt,t,x)
!!call pgpt(n,t,x,1)

!!call pgenv(taugrid(1),taugrid(Npsi),minval(psi),maxval(psi))
!!call pgline(Npsi,taugrid,psi)

!!call pgenv(t(1),t(Nt),minval(xop),maxval(xop))
!!call pgline(Nt,t,xop)
!stop
!!call pgend

!endprogram






!! New version 12 th march. Fixed bug: Now the convlved xgrid op reduces to 0 at the start
!... (where the edge effects occur) This was previously returning nans

!! UPDATE 1/8/14 removed the firstcal feature which fails if subsequent calls use different length input arrays
!! subroutine to convolve an input xgrid with a kernel psigrid each with equal grid spacing
!! op xop(Ntgrid)
!SPACING OF INTPUT TIME SERIES AND THE KERNEL MUST BE EQUAL
!! INPUTS Ntgrid, Ntaugrid size of input time series and kernel respectively
!! 		  tgrid,xgrid input time series.      taugrid,psigrid input kernel

!! OUTPUTS xgridop the convolution of the input time series with kernel

subroutine convolve(Ntgrid,Ntaugrid,tgrid,xgrid,taugrid,psigrid,xgridop)

implicit none

integer Ntgrid,Ntaugrid,idx,itau,i,it
real tgrid(Ntgrid),xgrid(Ntgrid),taugrid(Ntaugrid),psigrid(ntaugrid),xgridop(Ntgrid)
logical,save:: firstcall=.true.
real,save:: deet
double precision xgridopsum,taunorm
!integer,save,allocatable:: tauidx(:)
integer tauidx(Ntaugrid)

!if (firstcall) then
!write(*,*)'yay'
!allocate(tauidx(Ntaugrid))
deet=(tgrid(Ntgrid)-tgrid(1))/(Ntgrid-1.)

do i=1,Ntaugrid
tauidx(i)=taugrid(i)/deet
enddo

!firstcall = .false.

!endif



!!!!!! perform the convolution
do it=1,ntgrid  !! go from the index + the min tau to grid idx Nt - max tau
xgridopsum=0.d0

idx=1
itau=1
!write(*,*) tauidx
!stop
taunorm=0
do itau=1,ntaugrid
!if (it.eq.1) taunorm=taunorm+psigrid(itau)!*deet !! only need to do this once (normalise convolution kernel)

idx=it-tauidx(itau)
!write(*,*) itau,ntaugrid,idx,xgridopsum,taugrid(itau),psigrid(itau)
if (idx .lt. 1 .or. idx .gt. ntgrid) cycle
taunorm=taunorm+psigrid(itau)
xgridopsum=xgridopsum+psigrid(itau)*xgrid(idx)
!itau=itau+1
!stop

enddo
xgridop(it)=xgridopsum/taunorm!deet/taunorm
!write(*,*) tgrid(it),xgrid(it),xgridop(it),taunorm,xgridopsum
if (xgridop(it) .ne. xgridop(it)) then
xgridop(it) =0.0 !! set to zero if outside range of input
!write(*,*) 'problem with convolve.f90'
!write(*,*) xgridopsum,taunorm,itau
!stop
endif

enddo

!do itau=1,Ntaugrid
!write(*,*) taugrid(itau),psigrid(itau)
!enddo

return
endsubroutine







!!!!! subroutine to a light cure and convolution kernel and convolve them (simpler faster version than convovlve)

!!! UNTESTED!!!!!!!

!Notes light curve and kernel must be equally spaced in time before running this subroutine
!If not, should interpolate beforehand

!! inputs
! t(NT),x(NT) light curve, spacing dt
! psitau(Ntau) conolution kernel (transfer function)



!! Outputs
! echo(Ntau+NT)  !! first elements will be zero until idx - Ntau > 0


subroutine convolve2(x,t,dt,psitau,NT,Ntau)

integer NT,Ntau
real x(NT),t(NT),psitau(Ntau),dt,echolc(NT+Ntau)

Necho=NT+Ntau
echolc(:)=0.0

do idx=Ntau+1,Necho
sum=0.d0


do itau=1,Ntau
ilb=idx-itau
sum=sum+psitau(itau)*x(ilb)
enddo

echolc(idx)=sum*dt

enddo



end subroutine


















!! update 4th may14 added the feature to output temperature radius information to a text file called
! 'disksimtrop.op' as columns in the form r(rs), tvisc^4, tirad^4 (set outputtr to true to generate the output file)
! AGN parameters are saved in 'disksimtropparms.dat'



!version 1: input luminosity distance return fnu
!version 2: input H0 ,omega_m,omega_lam, return fnu


! faster version of diskim. Uses unique radius positions only on the radius grid
! Note: 29th April 14, using 1e5 and 1e6 nr produces similar results differing by
! ~ 1e-4 mJy at 1000-3000 Angstroms (Can prob get away with using 1e5 radii for speed)


! not over all azimuths
!! outputs are in mJy
!Ud ouput if version 2, input if version 1 new 11/08/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine disksim(uD,h0,omega_m,omega_lam,um,umdot,urlo,urhi,nr,uinc,k,alpha,eff,alb,z,uhx,w,fnu,iversion)

real cgs,ranu,bnu,k,dl_zoo,uD ! functions
logical outputtr /.True./, fastmode/.false./

iseed=13123414
!!! sub parms
!urlo=3.0
!urhi=10000.0
!nr=1000


fastmode = .false.
!if ((eff .eq. 0) .and. (alpha .eq. 0) .and. (k .eq. 0)) fastmode = .true.




vkms = 1.e5
pc = cgs( 'PC' )
fmjy = cgs( 'MJY' )


day = 24. * 3600.
year = cgs( 'YR' )
emsun = cgs( 'MSUN' )
angstrom = cgs( 'ANGSTROM' )


! physics constants
pi = 4. * atan2( 1., 1. )
twopi=pi*2
c = cgs( 'C' )
hplanck = cgs( 'H' )
bolt = cgs( 'K' )
gnewt = cgs( 'G' )
sigma = cgs( 'SIGMA' )



! cosmo parms
hubcon = h0*1*10**3/(30.857e21) !70kms-1 Mpc-1= 70*100000ms-1 Mpc-1 =70*100000/(pc*1e6)  ! checked and correct
!D = 1.*c * 1.*z /(1.*hubcon)  !! luminosity distance in centimeters
if (iversion .eq. 1) then
D=uD*1e6*pc
elseif (iversion .eq. 2) then
D= dl_zoo( z, omega_m, omega_lam )*c/hubcon
uD = D/1.e6/pc
endif


!! units
r1 = c * day
d1 = 1.e6 * pc
h1 = vkms / d1
emsunyr = emsun / year
rs=2*gnewt*um*emsun/c/c


!write(*,*) rs/r1, um, umdot
!stop

!write(*,*) D, omega_m,omega_lam,hubcon,z,um,umdot
!write(*,*) h0,omega_m,omega_lam,um,umdot,urlo,urhi,nr,uinc,k,alpha,eff,alb,z,uhx,w
!stop
radinc=uinc*pi/180
ci=cos(radinc)
si=sin(radinc)
hx=uhx*rs
t0= alog10(3*gnewt*umdot*um) - alog10(8*pi*sigma) + alog10(emsunyr) + alog10(emsun)   !! constant T terms in log10



dlum= eff*umdot*c*c ! driving luminosity in solar mass units
t1= alog10(1-alb)-alog10(4*pi*sigma)+alog10(dlum)+alog10(emsunyr)! +alog10(hx)  !! cinstant viscous terms


!!!!! radius spacing
udr= (urhi-urlo)/nr !!spacing between adjacent radii in rs
udlnr=(alog(urhi)-alog(urlo))/(1.*nr) ! spacing in ln r
!dr=udr*rs			!! '' in cm
r0=exp(alog(urlo))*rs

add=0.d0

!! optional bit (saving parameters)
if (outputtr) then
open(unit=1,file='disksimtropparms.dat')
write(1,*) rs,uinc,um,umdot,eff
close(1)
open(unit=1,file='disksimtrop.dat')
open(unit = 2, file = 'disksimsave.dat')
endif
!! end of optional bit

do ir=1,nr   !! radius loop
uran=ranu(0.0,udlnr,iseed)
ur=exp(alog(urlo)+udlnr*(ir-1.)+uran)
r=ur*rs
dr=r-r0
twopidr = twopi*dr
r0=r
rsem=exp(alog(urlo)+udlnr*(ir-0.5))*rs !! half radius between grid positions
!da=dr/rsem
!dca=cos(da)
!na=twopi/da


!uran=ranu(0.0,udr,iseed)
r=ur*rs    !! random grid pos

!!! for calculating viscous temperature
f=1.0-sqrt(urlo/ur)
tv4ln = t0-3*alog10(r)
!write(*,*) 'tc',tv4ln,ir,ia,r/rs,t0
tv4 = f*10**tv4ln     !!! viscous temperature ^4


if (fastmode .eqv. .false.) then
dh = k*(r/rs)**alpha*rs
!
dhd = alpha*dh/r
cpsi=sin(atan(dhd)-atan((dh-hx)/r)) ! angle of illumination
x=hx - dh
rstar = sqrt(x*x + r*r)   !!! lampost to surace element length
ti4log=t1- 2*alog10(rstar)+alog10(cpsi)  !! irradiating temperature in log10
ti4=10**ti4log
temp = sqrt(sqrt(ti4 + tv4)) !! temperature at radius
else
temp = sqrt(sqrt(tv4))
endif





if (outputtr) write(1,*) r/r1,tv4,ti4, dh, cpsi, x/rstar
idx=1



plank=bnu(w,temp)   !! plank function


if (fastmode .eqv. .false.) then
dsa=rsem*(ci-dhd*si)*2*pi*dr  !! solida angle element
else
dsa = rsem*ci*twopidr
endif

if (outputtr) write(2,*) w, r/r1, Temp, ci, dsa/r1/r1, plank


add =add+dsa*plank
!!! do away with azimuth loop
!do ia=1,na   !! azimuth loop

!dsa=rsem*(ci-dhd*si)*da*dr  !! solida angle element

!plank=bnu(w,temp)   !! plank function

!add =add+dsa*plank

!enddo !! end ia

enddo !! end ir

if (outputtr) close(1)
fnu=add/D/D/fmjy !! the fnu at the desired wavelength in mjy


if (outputtr) close(2)

return





end subroutine






















































!! update 4th may14 added the feature to output temperature radius information to a text file called
! 'disksimtrop.op' as columns in the form r(rs), tvisc^4, tirad^4 (set outputtr to true to generate the output file)
! AGN parameters are saved in 'disksimtropparms.dat'



!version 1: input luminosity distance return fnu
!version 2: input H0 ,omega_m,omega_lam, return fnu


! faster version of diskim. Uses unique radius positions only on the radius grid
! Note: 29th April 14, using 1e5 and 1e6 nr produces similar results differing by
! ~ 1e-4 mJy at 1000-3000 Angstroms (Can prob get away with using 1e5 radii for speed)


! not over all azimuths
!! outputs are in mJy
!Ud ouput if version 2, input if version 1 new 11/08/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! note ranoff option turns off the random radius position This may yield less accurate disk but is good for testing mcmc diffdisk_mcmc code

subroutine disksim_dp(uD,h0,omega_m,omega_lam,um,umdot,urlo,urhi,nr,uinc,k,alpha,eff,alb,z,uhx,w,fnu,iversion)

implicit none


double precision uD,h0,omega_m,omega_lam,um,umdot,urlo,urhi,&
uinc,k,alpha,eff,alb,z,uhx,w,vkms,pc ,fmjy,day,year,emsun,angstrom,pi,twopi,c,hplanck,&
bolt,gnewt,sigma,add,hubcon,D,r1,d1,h1,emsunyr,rs,radinc,ci,si,hx,t0,dlum,t1,&
udr,udlnr,r0,uran,ur,r,dr,twopidr,rsem,f,tv4ln,tv4,dh,dhd,cpsi,x,rstar,ti4log,ti4,&
temp,plank,dsa,fnu






real dl_zoo,cgs,ranu,bnu

integer iseed,nr,ir,iversion,idx
logical outputtr /.False./, fastmode/.false./, ranoff/.true./

iseed=13123414
!!! sub parms
!urlo=3.0
!urhi=10000.0
!nr=1000


fastmode = .false.
!if ((eff .eq. 0) .and. (alpha .eq. 0) .and. (k .eq. 0)) fastmode = .true.


write(*,*) uD,h0,omega_m,omega_lam,um,umdot,urlo,urhi,nr,uinc,&
k,alpha,eff,alb,z,uhx,w,fnu,iversion

vkms = dble(1.e5)
pc = dble(cgs( 'PC' ))
fmjy = dble(cgs( 'MJY' ))
day = dble(24. * 3600.)
year = dble(cgs( 'YR' ))
emsun = dble(cgs( 'MSUN' ))
angstrom = dble(cgs( 'ANGSTROM' ))


! physics constants
pi = dble(4. * atan2( 1., 1. ))
twopi=dble(pi*2)
c = dble(cgs( 'C' ))
hplanck = dble(cgs( 'H' ))
bolt = dble(cgs( 'K' ))
gnewt = dble(cgs( 'G' ))
sigma = dble(cgs( 'SIGMA' ))



! cosmo parms
hubcon = h0*1*10**3/(30.857e21) !70kms-1 Mpc-1= 70*100000ms-1 Mpc-1 =70*100000/(pc*1e6)  ! checked and correct
!D = 1.*c * 1.*z /(1.*hubcon)  !! luminosity distance in centimeters
if (iversion .eq. 1) then
D=uD*1e6*pc
elseif (iversion .eq. 2) then
D= dl_zoo( real(z), real(omega_m), real(omega_lam) )*c/hubcon
uD = D/1.e6/pc

endif


!! units
r1 = c * day
d1 = 1.e6 * pc
h1 = vkms / d1
emsunyr = emsun / year
rs=2*gnewt*um*emsun/c/c

!write(*,*) D, omega_m,omega_lam,hubcon,z,um,umdot
!write(*,*) h0,omega_m,omega_lam,um,umdot,urlo,urhi,nr,uinc,k,alpha,eff,alb,z,uhx,w
!stop
radinc=uinc*pi/180
ci=cos(radinc)
si=sin(radinc)
hx=uhx*rs
t0= dlog10(3*gnewt*umdot*um) - dlog10(8*pi*sigma) + dlog10(emsunyr) + dlog10(emsun)   !! constant T terms in log10



dlum= eff*umdot*c*c ! driving luminosity in solar mass units
t1= dlog10(1-alb)-dlog10(4*pi*sigma)+dlog10(dlum)+dlog10(emsunyr)! +dlog10(hx)  !! cinstant viscous terms


!!!!! radius spacing
udr= (urhi-urlo)/nr !!spacing between adjacent radii in rs
udlnr=(dlog(urhi)-dlog(urlo))/(1.*nr) ! spacing in ln r
!dr=udr*rs			!! '' in cm
r0=exp(dlog(urlo))*rs

add=0.d0

!!! optional bit (saving parameters)
!if (outputtr) then
!open(unit=1,file='disksimtropparms.dat')
!write(1,*) rs,uinc,um,umdot,eff
!close(1)
!open(unit=1,file='disksimtrop.dat')
!endif
!!! end of optional bit

do ir=1,nr   !! radius loop
if (ranoff) then
uran=udlnr
else
uran =dble(ranu(0.0,real(udlnr),iseed))
endif

ur=real(exp(dlog(urlo)+udlnr*(ir-1.)+uran))
r=ur*rs
dr=r-r0
twopidr = twopi*dr
r0=r
rsem=exp(dlog(urlo)+udlnr*(ir-0.5))*rs !! half radius between grid positions
!da=dr/rsem
!dca=cos(da)
!na=twopi/da


!uran=ranu(0.0,udr,iseed)
r=ur*rs    !! random grid pos

!!! for calculating viscous temperature
f=1.0-sqrt(urlo/ur)
tv4ln = t0-3*dlog10(r)
!write(*,*) 'tc',tv4ln,ir,ia,r/rs,t0
tv4 = f*10**tv4ln     !!! viscous temperature ^4



if (fastmode .eqv. .false.) then
dh = k*(r/rs)**alpha*rs
!
dhd = alpha*dh/r
cpsi=sin(atan(dhd)-atan((dh-hx)/r)) ! angle of illumination
x=hx - dh
rstar = sqrt(x*x + r*r)   !!! lampost to surace element length
ti4log=t1- 2*dlog10(rstar)+dlog10(cpsi)  !! irradiating temperature in log10
ti4=10**ti4log
temp = sqrt(sqrt(ti4 + tv4)) !! temperature at radius
else
temp = sqrt(sqrt(tv4))
endif



!if (outputtr) write(1,*) ur,tv4,ti4, dh, cpsi, x/rstar
idx=1



plank=dble(bnu(real(w),real(temp)))   !! plank function

if (fastmode .eqv. .false.) then
dsa=rsem*(ci-dhd*si)*2*pi*dr  !! solida angle element
else
dsa = rsem*ci*2*pi*dr
endif

add =add+dsa*plank
!!! do away with azimuth loop
!do ia=1,na   !! azimuth loop

!dsa=rsem*(ci-dhd*si)*da*dr  !! solida angle element

!plank=bnu(w,temp)   !! plank function

!add =add+dsa*plank

!enddo !! end ia

enddo !! end ir

if (outputtr) close(1)
fnu=add/D/D/fmjy !! the fnu at the desired wavelength in mjy


write(*,*) fnu,'out'
return

end subroutine































































! faster version of diskim. Uses unique radius positions only on the radius grid
! not over all azimuths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! random angles version i.e rnadom sampling of radi over all azimuths and radi
!! slower version but more accurate
!! enter observed wavlength, get observed spectrum
subroutine disksimobs(h0,um,umdot,urlo,urhi,nr,uinc,k,alpha,eff,alb,z,uhx,wobs,fnu)

real cgs,ranu,bnu,k ! functions
logical fastmode

fastmode = .false.
if ((eff .eq. 0) .and. (alpha .eq. 0) .and. (k .eq. 0)) fastmode = .true.


iseed=13123414
!!! sub parms
!urlo=3.0
!urhi=10000.0
!nr=1000

!!! convert between observed and emitted wavelength
w=wobs/(1.+z)

vkms = 1.e5
pc = cgs( 'PC' )
fmjy = cgs( 'MJY' )
day = 24. * 3600.
year = cgs( 'YR' )
emsun = cgs( 'MSUN' )
angstrom = cgs( 'ANGSTROM' )


! physics constants
pi = 4. * atan2( 1., 1. )
twopi=pi*2
c = cgs( 'C' )
hplanck = cgs( 'H' )
bolt = cgs( 'K' )
gnewt = cgs( 'G' )
sigma = cgs( 'SIGMA' )

! cosmo parms
hubcon = h0*1*10**3/(30.857e21) !70kms-1 Mpc-1= 70*100000ms-1 Mpc-1 =70*100000/(pc*1e6)  ! checked and correct
D = 1.*c * 1.*z /(1.*hubcon)  !! luminosity distance in centimeters

!! units
r1 = c * day
d1 = 1.e6 * pc
h1 = vkms / d1
emsunyr = emsun / year
rs=2*gnewt*um*emsun/c/c


radinc=uinc*pi/180
ci=cos(radinc)
si=sin(radinc)
hx=uhx*rs
t0= alog10(3*gnewt*umdot*um) - alog10(8*pi*sigma) + alog10(emsunyr) + alog10(emsun)   !! constant T terms in log10



dlum= eff*umdot*c*c ! driving luminosity in solar mass units
t1= alog10(1-alb)-alog10(4*pi*sigma)+alog10(dlum)+alog10(emsunyr)  !! cinstant viscous terms


!!!!! radius spacing
udr= (urhi-urlo)/nr !!spacing between adjacent radii in rs
udlnr=(alog(urhi)-alog(urlo))/(1.*nr) ! spacing in ln r
!dr=udr*rs			!! '' in cm
r0=exp(alog(urlo))*rs

add=0.d0
do ir=1,nr   !! radius loop
uran=ranu(0.0,udlnr,iseed)
ur=exp(alog(urlo)+udlnr*(ir-1.)+uran)
r=ur*rs
dr=r-r0
r0=r
rsem=exp(alog(urlo)+udlnr*(ir-0.5))*rs !! half radius between grid positions
da=dr/rsem
dca=cos(da)
na=twopi/da


uran=ranu(0.0,udr,iseed)
r=ur*rs    !! random grid pos

!!! for calculating viscous temperature
f=1.0-sqrt(urlo/ur)
tv4ln = t0-3*alog10(r)

tv4ln = f*10**tv4ln     !!! viscous temperature ^4


if (fastmode .eqv. .false.) then
dh = k*(r/rs)**alpha*rs
!
dhd = alpha*dh/r
x=hx - dh
rstar = sqrt(x*x + r*r)   !!! lampost to surace element length
cpsi=sin(atan(dhd)-atan((dh-hx)/r)) ! angle of illumination
ti4=t1- 2*alog10(rstar)+alog10(cpsi)  !! irradiating temperature in log10
temp = sqrt(sqrt(10**ti4 + tv4ln)) !! temperature at radius
else
temp = sqrt(sqrt(tv4ln))
endif



do ia=1,na   !! azimuth loop

if (fastmode .eqv. .false.) then
dsa=rsem*(ci-dhd*si)*da*dr  !! solida angle element
else
dsa = rsem*da*dr
endif

plank=bnu(w,temp)   !! plank function

add =add+dsa*plank

enddo !! end ia
!write(*,*) plank,temp,bnu(w,temp),r/rs,w,dsa,rsem,da*dr*rsem
enddo !! end ir

fnu=add/D/D/fmjy/(1.+z) !! the fnu at the desired wavelength in mjy


return

end subroutine















!!!!!!!! Subroutine to clculate the disk spectrum for a T(r) law like...
!! no lamppost component
!! T = T0 (r/r0)^-[3/4] ( ( 1-(rin/r)^0.5) / (1-(rin/r0)^0.5) )^1/4


subroutine disksim_2(uD,ch0,omega_m,omega_lam,r0_in,t0,um,urlo_in,urhi_in,nr,uinc,ak,&
alpha,eff,alb,z,uhx,w,fnu,iversion)

double precision fnuinside

!nr = 100000

eff = 0
vkms = 1.e5
pc = cgs( 'PC' )
fmjy = cgs( 'MJY' )
day = 24. * 3600.
year = cgs( 'YR' )
emsun = cgs( 'MSUN' )
angstrom = cgs( 'ANGSTROM' )


! physics constants
pi = 4. * atan2( 1., 1. )
twopi=pi*2
c = cgs( 'C' )
hplanck = cgs( 'H' )
bolt = cgs( 'K' )
gnewt = cgs( 'G' )
sigma = cgs( 'SIGMA' )


!! units
r1 = c * day
d1 = 1.e6 * pc
h1 = vkms / d1
emsunyr = emsun / year
rs=2*gnewt*um*emsun/c/c


if (r0_in .lt. 0) then
ur0 = abs(r0_in)*r1/rs
else
ur0 = r0_in
endif

if (urlo_in .lt. 0) then
urlo = abs(urlo_in)*r1/rs
else
urlo = urlo_in
endif

if (urhi_in .lt. 0) then
urhi = abs(urhi_in)*r1/rs
else
urhi = urhi_in
endif



T02 = t0*t0
T04 = t02*t02
a = 8*pi*sigma/(3*gnewt)
b = a*T04*ur0**3/emsunyr

f = (1. - (urlo/ur0)**0.5)

umdot = b*rs/um/emsun*rs**2*f

!write(*,*) b,rs,um,emsun,rs,f



call disksim_dp(dble(uD),dble(ch0),dble(omega_m),dble(omega_lam),dble(um),&
dble(umdot),dble(urlo),dble(urhi),nr,dble(uinc),dble(ak),&
dble(alpha),dble(eff),dble(alb),dble(z),dble(uhx),dble(w),&
fnuinside,iversion)

!write(*,*) dble(uD),dble(ch0),dble(omega_m),dble(omega_lam),dble(um),&
!dble(umdot),dble(urlo),dble(urhi),nr,dble(uinc),dble(ak),&
!!dble(alpha),dble(eff),dble(alb),dble(z),dble(uhx),dble(w),&
!fnuinside,iversion,'jfdhskfhdskfhsw disksim_2'

!read(*,*)




fnu = real(fnuinside)


end subroutine





!! UPDATE: 22nd April 2014: Removed interpolating onto evenly spaced grid. The iterated optimal
! scaling now fits using just the input time series x(NT), sampled at times t(NT), with errors
! er(NT).
! The interpolation was completely unnecessary, produced errors and slowed the code.


!!!!! FUROPT optimise starting fourier terms
!INPS t(NT),x(NT),er(NT),p(2*NW),w(NW),dtgridtemp=spacing between time points (i.e t(2)-t(1) ALL POINTS should be evenly spaced in time)
!OPS p(2*NW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine itoptscale(NT,NW,t,x,er,dtgridtemp,w,p,sigp)
integer NP
integer i,ip,it,ier,pgbeg
! parameter( maxt = 1000 )
real w(NW), p(2*NW),sigp(2*NW)
real x(NT), t(NT), er(NT), sigpone(2), pone(2),med,xdup(NT)
! real xdup(maxt), xtemp(maxt)
NP=2*NW
tgridtemplo=t(1)-0.01 !! start just a little to the left of the first fdata point


do i=1,NT
xdup(i)=x(i)
enddo


ip=1

do i=1,NW    !!!! furoptscal works for 1 freq, this loop performs optscal on
!! one frequency, subtracts this and continues
!! pone are the sin and cosine parameters for 1 frequency
!!furoptscal has a few lines to deal with 0 frequency input(i.e a constant)

ip1=ip+1

call furoptscal(Nt,w(i),t,xdup,er,pone,sigpone)
p(ip)=pone(1)
sigp(ip)=sigpone(1)
p(ip1)=pone(2)
sigp(ip1)=sigpone(2)
write(*,*) w(i),pone

do it=1,Nt
wt = w(i) * t(it)
fit = p(ip) * sin(wt) + p(ip1) * cos(wt)
xdup(it)=xdup(it) - fit
enddo



!write(*,*)i,Nw,w(i),pone(1),pone(2)
ip=ip+2     !! update the indicies of the p array


enddo      ! end NW loop




return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine















!! UPDATE: call furoptscal2 which incorporates a constant offset (doesn't require you to subtract the mean first)
!! UPDATE: 22nd April 2014: Removed interpolating onto evenly spaced grid. The iterated optimal
! scaling now fits using just the input time series x(NT), sampled at times t(NT), with errors
! er(NT).
! The interpolation was completely unnecessary, produced errors and slowed the code.


!!!!! FUROPT optimise starting fourier terms
!INPS t(NT),x(NT),er(NT),p(2*NW),w(NW),os,sigos (offset initial guesses and uncertainties)
!
!OPS p(2*NW)
! os,sigos offset and uncertainty
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine itoptscale2(NT,NW,t,x,er,dtgridtemp,w,p,sigp,os,sigos)
integer,intent(in):: NT,NW
real,intent(in):: t(NT), x(NT), er(NT), dtgridtemp, w(NW)
real,intent(out)::  p(2*NW),sigp(2*NW), os, sigos
integer i,ip,it,ier,pgbeg,NP

real sigpone(3), pone(3),med,xdup(NT)

NP=2*NW
tgridtemplo=t(1)-0.01 !! start just a little to the left of the first fdata point


do i=1,NT
xdup(i)=x(i)
enddo


ip=1

do i=1,NW    !!!! furoptscal works for 1 freq, this loop performs optscal on
!! one frequency, subtracts this and continues
!! pone are the sin and cosine parameters for 1 frequency
!!furoptscal has a few lines to deal with 0 frequency input(i.e a constant)

ip1=ip+1

call furoptscal2(Nt,w(i),t,xdup,er,pone,sigpone)
p(ip)=pone(2)
sigp(ip)=sigpone(2)
p(ip1)=pone(3)
sigp(ip1)=sigpone(3)

if (i .eq. 1) then
os=pone(1)
sigos=sigpone(1)
!write(*,*) 'os info', os,sigos
!read(*,*)
endif

!write(*,*) w(i),pone

do it=1,Nt
wt = w(i) * t(it)
fit = p(ip) * sin(wt) + p(ip1) * cos(wt) +pone(1)
xdup(it)=xdup(it) - fit
enddo



!write(*,*)i,Nw,w(i),pone(1),pone(2)
ip=ip+2     !! update the indicies of the p array


enddo      ! end NW loop




return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine














!! Tested and works 21 feb

!! iterated optimal scaling code to fit only a sin and cosine parameter
!! the constant is fitted by inputting a zero frequency
! data t(N),x(N),sig(N) t,x,uncertainty
! p(1) sin and p(2) cosine amplitudes with uncertainty sigp(1) and sigp(2)
! angular frequency w

!! outputs p(1) and p(2) changed values based on the updated parameter values

!!! OUTPUTS sigp(NP), p(NP) (P(NP) is also an input of default values).

subroutine furoptscal(N,w,t,x,sig,p,sigp)

implicit none
integer N,NP,iteration,nits,ip,ip2,ix,i
real x(N),sig(N),p(2),sigp(2),w,t(N)
real sw(N),cw(N),sigsq(N),b
double precision top,bot,top1


sigp(1:2)=0.0 !! set starting errors to 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! if you assign a zero frequency,(just fitting constant to stop /0)
if (w .eq. 0) then
top=0.0
bot=0.0
do i=1,N
b=1/(sig(i)*sig(i))
top=top+x(i)*b
bot=bot+b
enddo
p(2)=top/bot
sigp(2)=1/bot
p(1)=0.0
goto 21
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!! assign the sin and cosine parms
do i=1,N
sw(i)=sin(w*t(i))
cw(i)=cos(w*t(i))
sigsq(i)=sig(i)*sig(i)
enddo


nits=100   !!! change this in a minute

do iteration=1,nits                                               !!!! start iteration loop





!!!! sin parameter p(1)
top=0.0
bot=0.0
do i=1,N
b=sw(i)/sigsq(i)
top=top+(x(i)-p(2)*cw(i))*b
bot=bot+b*sw(i)
enddo
p(1)=top/bot
sigp(1)=1/bot


!!!! cos parameter p(2)
top=0.0
bot=0.0
do i=1,N
b=cw(i)/sigsq(i)
top=top+(x(i)-p(1)*sw(i))*b
bot=bot+b*cw(i)
enddo
p(2)=top/bot
sigp(2)=1/bot




enddo !!! end iteration  loop

!! uncertainty is  the square root of the variance
sigp(1)=sqrt(sigp(1))
sigp(2)=sqrt(sigp(2))
21	   return

end subroutine
















!! New 11th June alter to include a constant offset term in the ios x(t) =x0 + Ssin(wt) +Ccos(wt)

!! Tested and works 21 feb

!! iterated optimal scaling code to fit only a sin and cosine parameter
!! the constant is fitted by inputting a zero frequency
! data t(N),x(N),sig(N) t,x,uncertainty
! p(1) sin and p(2) cosine amplitudes with uncertainty sigp(1) and sigp(2)
! angular frequency w

!! outputs p(1) and p(2) changed values based on the updated parameter values

!!! OUTPUTS sigp(NP), p(NP) (P(NP) is also an input of default values).

subroutine furoptscal2(N,w,t,x,sig,p,sigp)

implicit none
integer N,NP,iteration,nits,ip,ip2,ix,i
real x(N),sig(N),p(3),sigp(3),w,t(N)
real sw(N),cw(N),sigsq(N),b
double precision top,bot,top1


sigp(1:3)=0.0 !! set starting errors to 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! if you assign a zero frequency,(just fitting constant to stop /0)
if (w .eq. 0) then
top=0.0
bot=0.0
do i=1,N
b=1/(sig(i)*sig(i))
top=top+x(i)*b
bot=bot+b
enddo
p(2)=top/bot
sigp(2)=1/bot
p(1)=0.0
goto 21
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!! assign the sin and cosine parms
do i=1,N
sw(i)=sin(w*t(i))
cw(i)=cos(w*t(i))
sigsq(i)=sig(i)*sig(i)
enddo


nits=1000   !!! change this in a minute

do iteration=1,nits                                               !!!! start iteration loop



!!  Fit the offset
top=0.0
bot=0.0
do i=1,N
b=1/sigsq(i)
top=top+(x(i)-p(2)*sw(i)-p(3)*cw(i))*b
bot=bot+b
enddo
p(1)=top/bot
sigp(1)=1/bot
!	   write(*,*) p(1),sigp(1),top,bot
!	   read(*,*)





!!!! sin parameter p(1)
top=0.0
bot=0.0
do i=1,N
b=sw(i)/sigsq(i)
top=top+(x(i)-p(1)-p(3)*cw(i))*b
bot=bot+b*sw(i)
enddo
p(2)=top/bot
sigp(2)=1/bot


!!!! cos parameter p(2)
top=0.0
bot=0.0
do i=1,N
b=cw(i)/sigsq(i)
top=top+(x(i)-p(1)-p(2)*sw(i))*b
bot=bot+b*cw(i)
enddo
p(3)=top/bot
sigp(3)=1/bot



!	   write(*,*) iteration,p,sigp
enddo !!! end iteration  loop
!	   stop
!! uncertainty is  the square root of the variance
sigp(1)=sqrt(sigp(1))
sigp(2)=sqrt(sigp(2))
sigp(3)=sqrt(sigp(3))


21	   return

end subroutine


































!! function to compute the average of a set of dat 12/04/14

!INPUT
!n			integer
!input x(n)  real*4

!! OUTPUT
!avg          real*4

function avg(x,n)

real x(n),avg
double precision sum

sum=0.d0
do i=1,n
sum=sum+x(i)
enddo

avg=sum/n

return
end function





!! function to compute the standard deviation of a dataset 12/04/14

!INPUT
!n			integer
!input x(n)  real*4

!! OUTPUT
!avg          real*4

function sdev(x,n)

real x(n),sdev
double precision sum

sum=0.d0
do i=1,n
sum=sum+x(i)
enddo

avg=sum/n


sum=0.d0
do i=1,n
a=x(i)-avg
sum=sum+a*a
enddo

if (n.gt. 1) then
sdev=sqrt(sum/(n-1))
else
sdev=0.0
write(*,*) 'entered only one point into sdev...'
stop

endif







return
end function













!! function to compute the standard deviation and averageof a dataset 12/02/15

!INPUT
!n			integer
!input x(n)  real*4

!! OUTPUT
! avg         real*4
! sdev
subroutine sdevave(x,n,avg,sdev)

real x(n),sdev,avg
double precision sum

sum=0.d0
do i=1,n
sum=sum+x(i)
enddo

avg=sum/n


sum=0.d0
do i=1,n
a=x(i)-avg
sum=sum+a*a
enddo

if (n.gt. 1) then
sdev=sqrt(sum/(n-1))
else
sdev=0.0
write(*,*) 'entered only one point into sdev...'
stop

endif







return
end subroutine



!subroutine to compute 1 iteration of the cream model


subroutine cream_model(Ntgrid,Ntau,Nwav,NK,tau,tgrid,wav,w,sk,ck,&
T1v,T1i,slope_v,slope_x,deginc,embh,emdot,&
os,st,thcent,thwide,xmod_d,xmod,psiop,redshift,&
ur0_in,diskk_in,diskka_in,alb_in,eff_in,rinsch_in,uhx_in)

integer,intent(in):: Ntgrid,Nwav,Nk,Ntau
real,intent(in):: tgrid(Ntgrid), tau(Ntau), wav(Nwav), w(Nk),sk(Nk),ck(Nk),T1v,T1i,slope_v,&
slope_x,deginc,embh,os(Nwav),st(Nwav),thcent(Nwav),thwide(Nwav)
real,intent(out):: xmod(Ntgrid,Nwav),xmod_d(Ntgrid), psiop(Ntau,Nwav)
real,allocatable,dimension(:,:),save::sinw,cosw
logical,save:: i_1stcall = .true.
double precision:: sum, sum2
real ur0_in,diskk_in,diskka_in,rinsch_in,alb_in,eff_in,uhx_in
logical diag /.False./


if(rinsch_in .ne. -1)then
rinsch = rinsch_in
else
rinsch = 3.0
endif

if(diskka_in .ne. -1)then
diskka = diskka_in
else
diskka = 0.0
endif

if(diskk_in .ne. -1)then
diskk = diskk_in
else
diskk = 0.0
endif


if(uhx_in .ne. -1)then
uhx = uhx_in
else
uhx = 3.0
endif

if(ur0_in .ne. -1)then
ur0 = ur0_in
else
ur0 = 3.0
endif

if(alb_in .ne. -1)then
alb = alb_in
else
alb = 0.0
endif

if(eff_in .ne. -1)then
eff = eff_in
else
eff = 0.1
endif


dtgrid = (tgrid(Ntgrid) - tgrid(1))/(ntgrid - 1)

!!evaluate sin and cosine grid (no need to do this more than once per target)
if (i_1stcall) then
allocate(sinw(NK,Ntgrid),cosw(NK,Ntgrid))
do it = 1,Ntgrid
tnow = tgrid(it)
!write(*,*) it, tnow
do iw = 1,Nk
wnow = w(iw)
sinw(iw,it) = sin(wnow*tnow)
cosw(iw,it) = cos(wnow*tnow)
enddo
enddo
i_1stcall = .false.
endif


!write(*,*) 'driving light curve computation',Ntgrid,Nk,T1v,T1i,slope_v,slope_x,embh,deginc,rinsch
!compute the driving light curve
do it = 1,Ntgrid
sum = 0.d0
do iw = 1,Nk
sum = sum + sk(iw)*sinw(iw,it) + ck(iw)*cosw(iw,it)
enddo
xmod_d(it) = sum
enddo


!since we have to look back the distance of the time grid, the first Ntau elements of
!tgrid must be ignored
itlo = Ntau+1


wavprev = -0.01
do ilc = 1,NWav


wavnow = wav(ilc)
osnow  = os(ilc)
stnow  = st(ilc)
thcentnow = thcent(ilc)
!write(*,*) thcent(1:Nwav),thcentnow,'thcentpars', T1v,T1i

if ((ilc .gt. 1) .and. (wavnow .eq. wavprev) .and. (wavnow .ne. -1.0)) then
psiop(1:Ntau,ilc) = psiop(1:Ntau,ilc-1)
if (diag) write(*,*) 'no convolve (repeat wavelength)',wav(ilc),ilc,Nwav
else
if (diag) write(*,*) 'convolution computation wav',wav(ilc),ilc,Nwav,&
thcentnow,thwide(ilc)

if ((wavnow .gt. 0) .and. (wavnow .lt. 100.0)) then
psiop(1:Ntau,ilc) = 0
else if ((wav(ilc) .eq. -1)) then
call tftophat(Ntau, tau, psiop(1:Ntau,ilc), thcentnow, thwide(ilc))
else if ((T1v .eq. 0) .and. (T1i .eq. 0)) then
!call tfb(tau,psiop(1:Ntau,ilc),&
!Ntau,wavnow,embh,emdot,uhx,eff*(1.-alb),&
!rinsch,deginc,slope_v,slope_x,ur0,diskk,diskka,redshift,1)
call tfbx(tau,Ntau,wavnow,&
-1*emdot,-1*emdot,slope_v,slope_x,&
embh,deginc,rinsch,psiop(1:Ntau,ilc),redshift)
!write(*,*) 'option 2 ',wavnow,T1v,T1i,slope_v,slope_x,embh,emdot,deginc,rinsch,redshift
!do itau = 1,Ntau
!if (tau(itau) .ge. 0) then
!idzero = itau
!exit
!endif
!enddo
!write(*,*) tau(idzero-1:idzero+10)
!write(*,*) psiop(idzero-1:idzero+10,ilc)
!write(*,*) 'TFBX CANNOT DEAL WITH -VE LAG GRID IF IN DISK MODE! FIX'
!stop
else
!write(*,*) 'option 3 ',wavnow,T1v,T1i,slope_v,slope_x,embh,deginc,rinsch,redshift
call tfbx(tau,Ntau,wavnow,T1v,T1i,slope_v,slope_x,embh,deginc,rinsch,psiop(1:Ntau,ilc),&
redshift)
endif
endif





!dont convolve if have driver
if ((wavnow .gt. 0) .and. (wavnow .lt. 100.0)) then

xmod(1:Ntgrid,ilc) = xmod_d(1:Ntgrid)*stnow + osnow

else

!no need to reconvolve if dealing with light curves at same wavelength
taulo = tau(1)
tauhi = tau(ntau)
dtau = (tauhi - taulo)/(Ntau-1)

itlo = max(1,floor(tauhi/dtau))
ithi = min(Ntgrid,Ntgrid - floor(abs(taulo)/dtau))

!write(*,*) itlo,ithi,Ntgrid
!stop

if ((ilc .gt. 1) .and. (wavnow .eq. wavprev) .and. (wavnow .ne. -1)) then
osprev = os(ilc-1)
stprev = st(ilc-1)
do it = itlo,Ntgrid
xmod(it,ilc) = (xmod(it,ilc-1) - osprev)*stnow/stprev + osnow
enddo
else
sum2 = 0.d0
do it = itlo,ithi
sum = 0.d0
psiprev = 0
do itau = 1,Ntau
taunow = tau(itau)
idtau = floor((taunow)/dtau)
idnow = max(1,it-idtau)
idnow = min(idnow,Ntgrid)
psinow = psiop(itau,ilc)
if (it .eq. itlo) sum2 = sum2 + psinow
if ((psinow .eq. 0) .and. (psiprev .gt. 0) ) exit
!write(*,*) idtau,it,idnow,'dsfsd',tgrid(it),tau(max(1,idtau))
!sum = sum + xmod_d(it-idtau+1)*psinow
sum = sum + xmod_d(idnow)*psinow
psiprev = psinow
enddo
!read(*,*)
xmod(it,ilc) = sum/sum2*stnow + osnow
enddo
endif


endif !end if wav < 100 > 0

wavprev = wavnow
enddo

if (diag) then
!do itau = 1,Ntau
!write(*,*) tau(itau), psiop(itau,1:Nwav)
!enddo
write(*,*) 'T1v, T1i, wavnow,embh,emdot,uhx,eff*(1.-alb),&
rinsch,deginc,slope_v,slope_x,ur0,diskk,diskka,redshift'

write(*,*) T1v, T1i, wavnow,embh,emdot,uhx,eff*(1.-alb),&
rinsch,deginc,slope_v,slope_x,ur0,diskk,diskka,redshift

write(*,*) 'os'
write(*,*) os(1:Nwav)
write(*,*) 'st'
write(*,*) st(1:Nwav)
write(*,*) 'wav'
write(*,*) wav(1:Nwav)
write(*,*) 'thcent'
write(*,*) thcent(1:Nwav)
write(*,*) 'thwide'
write(*,*) thwide(1:Nwav)

endif


!if diag true then plto things
if (diag) then
open(unit = 1, file = 'cream_model_diag.dat')
do it = 1,Ntgrid
write(1,*) tgrid(it),xmod_d(it),xmod(it,1:Nwav)
enddo
close(1)

open(unit = 1, file = 'cream_model_psi.dat')
do it = 1,Ntau
write(1,*) tau(it),psiop(it,1:nwav)
enddo
close(1)
write(*,*) 'diagnostic file in cream_model.f90 examine and the enter to resume code'
read(*,*)
endif


end subroutine























!! interface with cream to calculate the model
!! for a bunch of iterations
subroutine cream_modop(cov,ave,wav,embh,z,Nwav,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,&
NPscaleidx,NPthcentidx,NPtridx,Nits,Ntgrid,Ntau,w,tau,tgrid,xmodop,sigxmod)!psiop,sigpsi)


integer,intent(in):: Nwav,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,NPscaleidx,NPthcentidx,&
Nits,NP,Ntau
real,intent(in):: ave(NP),cov(NP,NP),wav(Nwav),tgrid(Ntgrid),tau(Ntau),embh,&
w(NPF/2),z
real,intent(out):: xmodop(Ntgrid,Nwav),sigxmod(NTgrid,Nwav)!, psiop(Ntau,Nwav),&
!sigpsi(Ntau,Nwav)
!,xmod_d(Ntgrid),xmod(Ntgrid,Nwav)
real p(NP,NP), u_gaus_0_v(NP,NP), eval(NP), evec(NP,NP), gaus_0_v(NP),pnew(NP)
real T1v,T1x,slope_v,&
slope_x,deginc,os(Nwav),st(Nwav),thcent(Nwav),thwide(Nwav)
real sk(NPF/2),ck(NPF/2),xmod(Ntgrid,Nwav),xmod_d(Ntgrid), psi(Ntau,Nwav)
double precision:: sum2_xmod(Ntgrid,Nwav),sum_xmod(Ntgrid,Nwav)

pi = 3.141592653589793238462643383
rad2deg = 180/pi
!,&
!sum2_psi(Ntau,Nwav),sum_psi(Ntau,Nwav)
!calculate the eigen values and vectors from the covariance matrix
call jacobi(cov,np,np,eval,evec,nrot)
do idx = 1,np
evnow = eval(idx)
eval(idx)= max(0.0,evnow)
enddo


!sum_psi(1:Ntau,1:Nwav) = 0
!sum2_psi(1:Ntau,1:Nwav) = 0
sum_xmod(1:Ntgrid,1:Nwav) = 0
sum2_xmod(1:Ntgrid,1:Nwav) = 0


!write covariance matrix out to file
!open(unit = 1,file = 'filecov.dat')
!do ip = 1,NP
!write(1,*) cov(1:NP,ip)
!enddo
!close(1)

iseed = -1
do iter = 1,Nits

write(*,*) 'creammodel evaluating steps',iter
!calculate  new parameter step from averages and covariance matrix evals and evecs
do ipc = 1,NP
gaus_0_v(ipc) = rang(0.,sqrt(eval(ipc)),iseed)
enddo
do ipr = 1,NP
evt = gaus_0_v(ipr)
do ipc = 1,NP
u_gaus_0_v(ipc,ipr) = evt*evec(ipc,ipr)
enddo
enddo
do ipc = 1,NP
sum = 0.d0
do ipr = 1,NP
sum = sum + u_gaus_0_v(ipc,ipr)
enddo
pnew(ipc) = sum + ave(ipc)

if (ipc .eq. NPcosincidx) then
do while ((pnew(ipc) .lt. 0) .or. (pnew(ipc) .gt. 1))
write(*,*) 'cosinc less than 0 or gt 1 recasting',ave(ipc),pnew(ipc)
pnew(ipc) = rang(0.5,0.25,iseed)
enddo
endif

if (ipc .ge. Nptridx .and. ipc .lt. NPtridx + Nwav) write(*,*) 'trparm',ave(ipc),pnew(ipc)
enddo


umdot = pnew(NPumbhidx)
deginc = acos(pnew(NPcosincidx))*rad2deg
idxp = 1
do ip = 1,NPF,2
sk(idxp)   = pnew(ip)
ck(idxp)   = pnew(ip+1)
idxp = idxp + 1
enddo

idx = 1
do ip = NPscaleidx,NPscaleidx+Nwav-1
st(idx) = pnew(ip)
idx = idx + 1
enddo

idx = 1
do ip = NPoffsetidx,NPoffsetidx+Nwav-1
os(idx) = pnew(ip)
idx = idx + 1
enddo

idx = 1
do ip = NPthcentidx,NPthcentidx+Nwav-1
thcent(idx) = pnew(ip)
thwide(idx) = pnew(ip+Nwav)
idx = idx + 1
enddo

! slope parameters
sv  = pnew(NPtridx)
si  = pnew(NPtridx+1)
T1v = pnew(NPtridx+2)
T1i = pnew(NPtridx+3)

!write(*,*) 'calling cream model'

call cream_model(Ntgrid,Ntau,Nwav,NPF/2,tau,tgrid,wav,w,sk,ck,&
T1v,T1i,sv,si,deginc,embh,umdot,&
os,st,thcent,thwide,xmod_d,xmod,psi,z,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)

!write(*,*) 'cream model done'



do ilc = 1,Nwav
do it = 1,Ntgrid
xm = xmod(it,ilc)
sum_xmod(it,ilc) = sum_xmod(it,ilc) + xm
sum2_xmod(it,ilc) = sum2_xmod(it,ilc) + xm*xm
!if (sum_xmod(it,ilc) .ne. sum_xmod(it,ilc)) then
!write(*,*) xm,'xmprob!!!!',sv,si,T1v,T1i,z,embh,deginc
!read(*,*)
!endif
enddo
!do itau = 1,Ntau
!xm = psi(itau,ilc)
!sum_psi(itau,ilc) = sum_psi(itau,ilc) + xm
!sum2_psi(itau,ilc) = sum2_psi(itau,ilc) + xm*xm
!enddo
enddo

enddo !end iteration


do ilc = 1,Nwav
do it = 1,Ntgrid
avexmod = sum_xmod(it,ilc)/Nits

sxm = sqrt( abs(sum2_xmod(it,ilc)/Nits - (avexmod*avexmod)) )

!if (sxm .ne. sxm) then
!write(*,*) 'problem with sigxmod',tgrid(it),wav(ilc),sum2_xmod(it,ilc)/nits,avexmod,avexmod**2
!endif
sigxmod(it,ilc) = sxm
xmodop(it,ilc) = avexmod
enddo
enddo

!do ilc = 1,Nwav
!do it = 1,Ntau
!avepsi = sum_psi(it,ilc)/Nits
!sigpsi(it,ilc) = sqrt( sum2_psi(it,ilc)/Nits - avepsi*avepsi )
!psiop(it,ilc)  = avepsi
!enddo
!enddo


!too few iterations to calculate error snakes
if (Nits .le. 2) then
sigxmod(1:Ntgrid,1:Nwav) = 0.0
endif

end subroutine
























!! interface with cream to calculate the model
!! for a bunch of iterations
subroutine cream_modop_ps(parsave,wav,embh,z,Nwav,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,&
NPscaleidx,NPthcentidx,NPtridx,Nits,Ntgrid,Ntau,w,tau,tgrid,xmodop,sigxmod)!psiop,sigpsi)


integer,intent(in):: Nwav,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,NPscaleidx,NPthcentidx,&
Nits,NP,Ntau
real,intent(in):: parsave(NP,Nits),wav(Nwav),tgrid(Ntgrid),tau(Ntau),embh,&
w(NPF/2),z
real,intent(out):: xmodop(Ntgrid,Nwav),sigxmod(NTgrid,Nwav)!, psiop(Ntau,Nwav),&
!sigpsi(Ntau,Nwav)
!,xmod_d(Ntgrid),xmod(Ntgrid,Nwav)
real pnew(NP)
real T1v,T1x,slope_v,&
slope_x,deginc,os(Nwav),st(Nwav),thcent(Nwav),thwide(Nwav)
real sk(NPF/2),ck(NPF/2),xmod(Ntgrid,Nwav),xmod_d(Ntgrid), psi(Ntau,Nwav)
double precision:: sum2_xmod(Ntgrid,Nwav),sum_xmod(Ntgrid,Nwav)
logical diag /.False./

pi = 3.141592653589793238462643383
rad2deg = 180/pi
!,&
!sum2_psi(Ntau,Nwav),sum_psi(Ntau,Nwav)

!sum_psi(1:Ntau,1:Nwav) = 0
!sum2_psi(1:Ntau,1:Nwav) = 0
sum_xmod(1:Ntgrid,1:Nwav) = 0
sum2_xmod(1:Ntgrid,1:Nwav) = 0


!write covariance matrix out to file
!open(unit = 1,file = 'filecov.dat')
!do ip = 1,NP
!write(1,*) cov(1:NP,ip)
!enddo
!close(1)

iseed = -1
do iter = 1,Nits

do ip = 1,NP
pnew(ip) = parsave(ip,iter)
enddo

umdot = pnew(NPumbhidx)
deginc = acos(pnew(NPcosincidx))*rad2deg

if (diag) then
write(*,*) umdot,deginc,wav,embh,z,Nwav,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx,&
NPscaleidx,NPthcentidx,NPtridx,Nits,Ntgrid,Ntau
write(*,*) 'umdot,deginc,wav,embh,z,Nwav,NP,NPF,NPumbhidx,NPcosincidx,NPoffsetidx'
write(*,*) 'NPscaleidx,NPthcentidx,NPtridx,Nits,Ntgrid,Ntau'
!stop
endif

idxp = 1
do ip = 1,NPF,2
sk(idxp)   = pnew(ip)
ck(idxp)   = pnew(ip+1)
idxp = idxp + 1
enddo

idx = 1
do ip = NPscaleidx,NPscaleidx+Nwav-1
st(idx) = pnew(ip)
idx = idx + 1
enddo

idx = 1
do ip = NPoffsetidx,NPoffsetidx+Nwav-1
os(idx) = pnew(ip)
idx = idx + 1
enddo

idx = 1
do ip = NPthcentidx,NPthcentidx+Nwav-1
thcent(idx) = pnew(ip)
thwide(idx) = pnew(ip+Nwav)
idx = idx + 1
enddo

! slope parameters
sv  = pnew(NPtridx)
si  = pnew(NPtridx+1)
T1v = pnew(NPtridx+2)
T1i = pnew(NPtridx+3)

if (diag) write(*,*) 'calling cream model',iter,umdot,deginc,sv,si,T1v,T1i

call cream_model(Ntgrid,Ntau,Nwav,NPF/2,tau,tgrid,wav,w,sk,ck,&
T1v,T1i,sv,si,deginc,embh,umdot,&
os,st,thcent,thwide,xmod_d,xmod,psi,z,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)

if (Nits .le. 2) then
xmodop(1:Ntgrid,1:Nwav) = xmod(1:Ntgrid,1:Nwav)
exit
endif

if (diag) write(*,*) 'cream model done'

Nitsplot = Nits

do ilc = 1,Nwav
do it = 1,Ntgrid
xm = xmod(it,ilc)

if (xm .ne. xm) then
Nitsplot = max(1,Nitsplot - 1)
exit
endif

sum_xmod(it,ilc) = sum_xmod(it,ilc) + xm
sum2_xmod(it,ilc) = sum2_xmod(it,ilc) + xm*xm
!if (sum_xmod(it,ilc) .ne. sum_xmod(it,ilc)) then
!write(*,*) xm,'xmprob!!!!',sv,si,T1v,T1i,z,embh,deginc
!read(*,*)
!endif
enddo
!do itau = 1,Ntau
!xm = psi(itau,ilc)
!sum_psi(itau,ilc) = sum_psi(itau,ilc) + xm
!sum2_psi(itau,ilc) = sum2_psi(itau,ilc) + xm*xm
!enddo
enddo

enddo !end iteration

if (Nits .ge. 2) then
do ilc = 1,Nwav
do it = 1,Ntgrid
avexmod = sum_xmod(it,ilc)/Nitsplot

sxm = sqrt( abs(sum2_xmod(it,ilc)/Nitsplot - (avexmod*avexmod)) )

!if (sxm .ne. sxm) then
!write(*,*) 'problem with sigxmod',tgrid(it),wav(ilc),sum2_xmod(it,ilc)/nits,avexmod,avexmod**2
!endif
sigxmod(it,ilc) = sxm
xmodop(it,ilc) = avexmod
enddo
enddo
endif

!too few iterations to calculate error snakes
if (Nits .le. 2) then
sigxmod(1:Ntgrid,1:Nwav) = 0.0
endif
!do ilc = 1,Nwav
!do it = 1,Ntau
!avepsi = sum_psi(it,ilc)/Nits
!sigpsi(it,ilc) = sqrt( sum2_psi(it,ilc)/Nits - avepsi*avepsi )
!psiop(it,ilc)  = avepsi
!enddo
!enddo

end subroutine










!! 9.915 Turns out it was right before, putting it back to old t0v4 t0x4
!! 3/9/15 fixed bug in parameterisation of T0v4 and T0x4. now, same values of apower and bpower will give the same temperature regardless of reference radius ur0 (if in lightdays) as they should
!! change 11th may 2014, altered so that the argument of dtdl =10**argument is incremented by a float (antiunderflow) to prevent underflows (antiunderflow is around 15.0) I then divide u by 10^15
!! fix 2.3 10th May 2014 bug fix solid angle element no multiplies by cos(inc) instead of (wrongly) previously dividing.
!! fix 2.2 28th march incorporate all paramaterisations of transfer function into a single code
!! fix 2.1 27th march 2014 do not multiply integrand by dtdl0 and dbdt0 until end of angle loop to prevent underflows
!! evaluation of dTdL saved until after integrand evaluation of dBdT split into before and after the integrand (cosh terms before dBdT0 after)


!! version 1: T(r) profile with parameteriosed as T^4 = T0v^4(r0/r)^(4b) +  T0x^4(r0/r)^(4a) Does not do flaring!!

!! version 2: T9r) profile with T paramaterised as t^4 = ([3G M Mdot] / [8pisigma r^3]) + L/ [8pisig ((h-z)^2+r^2)^1.5] [(alpha -1)h +z] / [(alpha^2h^2/r^2 +1 )]
!! version 2 works with flaring of the form



!! NOTE PROBLEM if fitting T=T0*(r/r0)**b, what is dT/dL (need this to work out psi)
!... can just do what I was previously doing and put in fiull model to get dT/dL but why would I then only inlcude the simpler T0(r/r0)**b T_r pofile.
!! simpler version of the transfer code delaysist3.
!... This code obtains a transfer function based on a power law temperature radius profile
!... of this form. T(r)=T0(r/r0)^b
! T0 can be anything, but advisable to stick to...
!. T0=3/(8pi stefbolt r0^3) G M Mdot + hx Lx /(4pi (r0^2 + hx^2)^3/2) i.e flat disk with irradiation

!! subroutine to evaluate delay distributions for input grid of evenly spaced time delay grid

!!!INPUTS inclination, lampost height rin rout delta ang etc see par file
!6000.0        !! wavelength (ANGSTROMS)
!1.0e8         !! embh    (m0)
!1.0           !! emdot (m0 / yr)
!uzx            !! lampost height (schqarzchild radii)
!alpha         !! disk height radius parameter (0 = flat disk) alpha
!k             !! constant of proportionality in above h,r law konst
!rin           !! inner radius
!alb           !! disk albedo
!eff           !! disk efficiency
!inc           !! disk inclination (inc)
! apower       !! lamppost component of tr profile
! bpower       !! viscous component of tr profile
!f0            !! the covering fraction at the reference radius of r0
!f0			   !! the power law slope change of the covering fraction with radius
!redshift
!Ntau number of time delas
!tau(Ntau) input delays at which to evaluate psi
!!! OUTPUTS psi(Ntau)

subroutine tfb(tau,psiop,Ntau,lamda_prez,embh,emdot_prez,uzx,eta,urin,inc,bpower,apower,ur0,&
ak,alpha,z,version)
implicit none

integer i,nr,ia,naz,iseed,idx,itau,Ntau,lo,hi,istartclock,iendclock
real lamda_prez,lamda,emdot,embh,emdot_prez,zx,eta,rin,k,alphainc,apower,bpower
real pi,twopi,fourpi,c,bk,h,s,emsun,day,yr,angst,cc,rssun,rs,um,uw,udmdt
real urr,pi4s,uTT,uT4d,uT4x,ux,deg2rad,tf0,oldi,rold,inc1,si,ci,zp1,z
real,save,allocatable::rsave(:),azsave(:),hsave(:)
real gnewt
real az,azhalf,aznew,daz,dr,dlnr,rnew,halfdr,rhalf,rrhalf,ddhalf,dhalf
real f,T4d,dd,d,cpsi,hh,hhh,ddd,T,TT,TTT,T4x,cgs,caz,cday,cdtau,dBdT,dBdT0
real dTdL,dTdL0,integrand,rr,rrr,tdel,ran3,r
double precision sum,gsum, top, bot
real gtemp,g(Ntau),psi(Ntau),tdelsubtau,sig_g,sig_g_sq,tau(Ntau)
real uzx,urin,urold,urnew,udr,udlnr,uhalfdr,urhalf,x,dtau,azold
real g_sdnum,logT4d0,logT4x0,logt4x,logt4d, psiop(Ntau), tauz(Ntau)
real logdTdL0,logdTdL,ucdtau,dBdT0_b,dsa,htemp,alphak,y1,xstop,add,diff
real zsubh,u,inc,uld,fnow
real plank,b,T0,r0,ur0,T0x4,T0v4,logTx4,logT0v4,logT0x4
integer version,choplo,chophi
real alpha,ak,ahor,ahorsq,zxsubh,ab,logab,antiunderflow,argument
logical info /.false./

call system_clock(istartclock)
!perform tests on redshift correction to check this is correct
zp1 = 1. + z
lamda = lamda_prez/zp1
emdot = emdot_prez!*zp1**3

if (info) then
open(unit = 5,file='tfb_trprof.dat')
endif

tauz(1:Ntau) = tau(1:Ntau)/zp1

psi(:)=0.0
nr=10000
udlnr=0.01       !! logarithmic step size set 0.001 for display purposes e.g papers etc
xstop=15         !! introduce a cutoff size of hc/(lamda kT) to stop evaluating psi
g_sdnum=5.0     !!! cut off tau summation after g_sdnum gaussian widths



!if( firstcall ) then      !!! calculate the constant parameters, only in 1st ieration
! physics constants
pi = 4. * atan2( 1., 1. )
twopi = 2. * pi
fourpi = 4. * pi
c = cgs( 'C' )
bk = cgs( 'K' )
gnewt = cgs( 'G' )
h = cgs( 'H' )
s = cgs( 'SIGMA' )
! astro constants
emsun = cgs( 'MSUN' )
day = 24. * 3600.
yr = cgs( 'YR' )
angst = cgs( 'ANGSTROM' )
! units
cc = c * c
rssun =2* gnewt * emsun / cc
rs = rssun*embh
!rs = rs*embh
um = emsun
uw = angst

udmdt = emsun / yr
urr = rs * rs
pi4s = fourpi * s
uTT = sqrt( udmdt / ( pi4s * urr * rs ) )



!disk temps

ux = h * c / ( uw * bk * lamda)

dBdT0=bk/(lamda*angst)**2   !!!! Rayleigh (hot) approximation to plank derivative
dBdT0_b=(h*c/(lamda*lamda*angst*angst))/bk*(h*c/(lamda*lamda*angst*angst))

! time delay
cday=c*day
dtau=(tau(Ntau)-tau(1))/(Ntau-1)
cdtau=c*dtau*day
ucdtau=cdtau/rs
!write(*,*) ucdtau,tau(2)-tau(1),dtau,c,day,dtau*c*day,rs,dtau*c*day/rs
!stop

uld=c*day


!! convolution parameters
sig_g=dtau/2
sig_g_sq=sig_g*sig_g



inc1=inc*pi/180
si=sin(inc1)
ci=cos(inc1)
tf0 = pi * bk * c / ( s * uw * uw * rs )
oldi = -1000.

udr=0!cdtau/rs
urold=urin
zx=uzx*rs         !! convert to cgs units
rin=urin*rs     !! ''

!!!!!!

!! more version specific terms (version 1)
!! calculate the T0v and T0i (evaluated at r0) according to flat disk theory
if (version .eq. 1) then

if (ur0 .lt. 0) then
r0 = abs(ur0)*uld
else
r0=ur0*rs
endif



T0v4=3/(8*pi)*10**(alog10(gnewt*embh)+alog10(emdot)&
-3*alog10(r0)+alog10(um)+alog10(udmdt)-alog10(s))



T0x4=1/(2*fourpi)*10**(alog10(zx)+alog10(eta*c*c)&
+alog10(emdot)+alog10(udmdt)-3*alog10(r0)-alog10(s))

!T0x4=1/fourpi*10**(alog10(zx)+alog10(eta*c*c)&
!+alog10(emdot)+alog10(udmdt)-3*alog10(r0)-alog10(s))
!! dtdl0 information
dTdL0=1./(16*pi)*10**(alog10(zx)-alog10(s)+(4*apower-3.)*alog10(r0)+15.0)
logdTdL0=(alog10(zx)-alog10(s)+(4*apower-3.)*alog10(r0)-alog10(16*pi))



else if (version .eq. 2) then
logT0v4=alog10(3/(pi*8))+(alog10(gnewt*udmdt)&
+alog10(embh)+alog10(emsun*emdot)-alog10(s))


logT0x4=alog10(eta/(pi*8))+alog10(zx)+2*alog10(c)&
+alog10(udmdt)+alog10(emdot)-alog10(s)

logdTdL0=alog10(1/(pi*32*s))
endif

!write(*,*) r0/rs, T0V4, T0X4, 'tfb info'

i=1   !! radius counter
idx=1 !! overall counter (counts all points)

!! start radius grid
!open(unit = 2, file='tfb_diag2.dat')
!write(*,*) ucdtau



if (info) then
write(*,*)
write(*,*) 'tfb input parms'
write(*,*) 'M (M0/1e7)', embh/1e7
write(*,*) 'Mdot (M0/yr)', emdot
write(*,*) 'Wavelength (Ang)', lamda
write(*,*) 'inclination (deg)', inc
write(*,*) 'lampost height (rs)', uzx
write(*,*) 'eta, urin,', eta, urin
write(*,*) 'T0v,T0x: ',sqrt(sqrt(T0v4)), sqrt(sqrt(T0x4))

if (version .eq. 1) then   !! version specific information
write(*,*)  'lampost power, viscous power, ur0',apower, bpower,ur0
else if (version .eq. 2) then
write(*,*) 'flare parms ak,alpha (h(r)=ak*r^alpha): ',ak,alpha
endif
write(*,*) 				   !!
endif


do i =1,nr !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Choose log or linear steps
!!! radius grid
if (udr .lt. ucdtau) then  !! log
urnew=urold*exp(udlnr)
udr=urnew-urold        !! save the radius spacing when stepping in log
dr=udr*rs

else                       !! linear
udr=ucdtau
urnew=urold+udr
!write(*,*) i, 'liner takeover',rs
!stop
end if

uhalfdr=udr/2              !! evaluate solid angles half way between 2 grid radii
urhalf=urold+uhalfdr
rhalf=urhalf*rs

rrhalf = rhalf*rhalf
!write(*,*) 'radius info',urnew,i,ucdtau,udr

!write(*,*) i,urold,urnew,udr,ucdtau
urold=urnew !! record the new old radius


!fnow = f0*(rhalf/r0)**f0slope

!!! azimuth grid       !!!! see below
daz=udr/urhalf
naz=twopi/daz          !!!! this can all be dealt with in the first call
daz=twopi/(naz)
!!! define solida angle elemnt prior to entering az loop
!write(*,*) naz,twopi,daz,nr,'number az tfb'
!! for flat disk, solid angle term can be evaluated outside the azimuth loop
if (version .eq. 1) dsa=daz*urhalf*udr*ci  !! solid angle element






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (version .eq. 1) then
zsubh=zx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! temperature radius info
!! viscous disk

T=sqrt(sqrt( T0v4*(r0/rhalf)**(bpower*4)+T0x4*(r0/rhalf)**(apower*4) ))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate log dtdl ... dtdL=(1-alb)/(8 pi sigma d^(3/2)) * ((alpha-1)h + zx) / sqrt((dh/dr)^2+1) / 4 T^3
argument=logdTdL0-3*alog10(T) -4*apower*alog10(rhalf)
!! new 11th may: prevent underflows
antiunderflow=40.0
dTdL=10**(argument+antiunderflow) !! can cause underflows for large R0 values
else if (version .eq. 2) then
h=ak*rhalf**alpha
zsubh=zx-h
ddhalf=zsubh*zsubh+rrhalf
dhalf=sqrt(dd)
ahor=alpha*h/rhalf
ahorsq=ahor*ahor
ab=((alpha-1.)*h+zx)/(ahorsq+1.) !!! collection of terms for the flared tr profile
logab=alog10(ab)
logTx4=logab-3*alog10(dhalf)
T=sqrt(sqrt(10**(logT0v4-3*alog10(rhalf)) + 10**logTx4))
if (T.ne.t)then
write(*,*) 'problem with T tfb',t0v4,rhalf,t0v4/rhalf**3,ahorsq+1.,ahorsq+1
stop
endif
dTdL=10**(logdTdL0-3*(alog10(dhalf)+alog10(T))+logab)
!!! define solida angle elemnt prior to entering az loop
htemp=ak*(rhalf)**alpha
dsa=urhalf*sqrt((htemp*alpha/(rhalf))**2+1)*udr*daz  !! solid angle element
endif !! end the version specific loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! evaluate plank function
x=ux/T        !! record h c /(lamda K T)

if (x .lt. xstop) then
integrand = fourpi*x*x/(cosh(x)-1.)    ! the xx/(cosh(x)-1) is part of dBdT !! the integrand of the transfer function


else
exit
endif      !!! stop the loops when plank becomes negligible


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! evaluate different radii and azimuths by taking random steps
azold=0.0

if (info) then
write(5,*)urold*rs/uld,sqrt(sqrt( T0v4*(r0/rhalf)**(bpower*4))),&
sqrt(sqrt(T0x4*(r0/rhalf)**(apower*4)))
endif

do ia=1,naz

aznew=azold+daz
azold=aznew
az=aznew+ran3(iseed)*daz
azhalf=aznew+daz/2
caz=cos(az)

r=(urold+ran3(iseed)*udr)*rs  !!! choose a random element between rold and rnew

rr=r*r
dd=(zsubh*zsubh+rrhalf)
d=sqrt(dd)
tdel=( zsubh*ci - r*caz*si + d )/cday !! calculate time delay in days


if (dTdL.eq.0) then
write(*,*)'dtdl is 0',i, logdTdL0, d,dd,rr,T,-3*(alog10(d)+alog10(T))+logab,sqrt(dd*d)*((alpha-1.)*h+zx)
exit
endif


!!!!!!!!!!


!!!!!the convolution gaussian approximation

!!! convolution function

add=g_sdnum* sig_g
lo=max(1,floor((tdel - add)/dtau)) !!lower index of tau
hi=min(Ntau,ceiling((tdel + add)/dtau)) !!upper index of tau

u=integrand*dsa*dTdL
do itau=lo,hi
tdelsubtau=tdel-tau(itau)

gtemp=exp(-0.5*tdelsubtau*tdelsubtau/sig_g_sq)
psi(itau)=psi(itau)+u*gtemp
!g(itau)=gtemp
!gsum=gsum+gtemp
enddo





!*fnow!*dBdT0                            !!!!! version 1.2 27th March 2014,
!if (argument .lt. -40.0) u=u/10**(antiunderflow)





idx=idx+1





!if (tdel .le. 5) write(2,*) r/uld, az, tdel!, x, t, dTdL, dbdt


enddo !! end azimuth




if (dTdL .eq.0) then
write(*,*) dTdL ,'dtdl bad'
exit
endif
!write(*,*) Ntau,tdel,r/rs,i,dtdl,x,integrand*dsa,'indicees'
if (x .gt. xstop) exit !! exit if we have reached the cutoff x
enddo !! end radius


!!!!!!! the edges are turned into nans, just set to zero if this happens
do itau=1,Ntau
if (psi(itau) .ne. psi(itau)) psi(itau)=0
enddo
!!!!!


!!!!! PROBLEM Long delay info is all being shoved into last bin, solve this
!! for now by making it equal to the previous datapoint. Could also linearly interpolate but should not matter,
!! it is very small
psi(Ntau)=psi(Ntau-1)

!close(2)
if (info) then
call system_clock(iendclock)
write(*,*)'Run time...',iendclock-istartclock
top = 0.d0
bot = 0.d0
do itau = 1,Ntau
top = top + psi(itau)*tau(itau)
bot = bot + psi(itau)
enddo

write(*,*) r0/rs, T0V4, T0X4, 'tfb info. Taumean =',top/bot
write(*,*)
!stop
!read(*,*)
endif
!return


!write(*,*) 'done before norm', psi(3)
!! normalise so area = 1
sum = 0.d0
do itau = 1,Ntau
sum = sum + psi(itau)
enddo
psi(1:Ntau) = psi(1:Ntau)/sum

!write(*,*) 'done', psi(3)
!do itau = 1,Ntau
!write(*,*) tau(itau), psi(itau)
!enddo

!interpolate to account for redshift effect
if (z > 0.0) then
call itp1d(psiop,choplo,chophi,tau,tauz,psi,Ntau,Ntau,0)
else
psiop(1:Ntau) = psi(1:Ntau)
endif




if (info) close(5)

end subroutine










!!gflags tfb_redshifttest_prog.f90 tfb_redshifttest.f90 ../cgs.f90 ../ran_new2.f90 ../itp2.for
!
!program testtf14
!integer pgopen
!real,allocatable:: psi(:),tau(:),psiz(:),wavsave(:),taumeansave(:)
!real x2(2), y2(2)
!
!!rm393
!!0.083311345429
!!62185670.0381
!
!taulo  = 0.0
!tauhi  = 100.0
!dtau   = 0.1
!urin   = 3.0
!urout  = 1000.0
!uembh  = 1.e8!6.218e7
!uemdot = 1.0!0.083311345429
!uemdotref = 1.0
!deginc = 0.0
!uhx    = 3.0
!eta    = 0.1
!alb    = 0.0
!wav    = 4770.0
!z = 0.0!0.584504203112
!
!slopevisc = 0.75
!slopeirad = slopevisc
!!slopeirad = 0.75
!Ntau = (tauhi - taulo)/dtau + 1
!
!
!wavlo = 1000.0
!wavhi = 10000.0
!dwav = 100.0
!Nwav = ceiling((wavhi - wavlo)/dwav + 1)
!
!allocate(tau(Ntau),psi(Ntau),psiz(Ntau),wavsave(Nwav),taumeansave(Nwav))
!
!do itau = 1,Ntau
!tau(itau) = (itau-1)*dtau
!enddo
!
!
!
!
!
!do iw = 1,Nwav
!wavnow = wavlo + (iw - 1)*dwav
!wavsave(iw) = wavnow
!write(*,*) 'calculating response for wavelenth',wavnow
!call tfb(tau,psi,Ntau,wavnow,uembh,uemdot,uhx,eta,urin,deginc,slopevisc,slopevisc,urin,0.0,0.0,0.0,1)
!
!bot = 0.d0
!top = 0.d0
!do itau = 1,Ntau
!top = top + psi(itau)*tau(itau)
!bot = bot + psi(itau)
!enddo
!
!taumeansave(iw) = top/bot
!enddo
!
!open(unit = 1,file = 'tf_mean.dat')
!do iw = 1,Nwav
!write(1,*) wavsave(iw), taumeansave(iw)
!enddo
!close(1)
!
!
!
!call system_clock(istart)
!call tfb(tau,psi,Ntau,wav,uembh,uemdot,uhx,eta,urin,deginc,slopevisc,slopevisc,urin,0.0,0.0,0.0,1)
!call system_clock(iend)
!call tfb(tau,psiz,Ntau,wav,uembh,uemdotref,uhx,eta,urin,deginc,slopevisc,slopeirad,urin,0.0,0.0,z,1)
!
!
!pmax = maxval(psi)
!pzmax = maxval(psiz)
!psi(1:Ntau) = psi(1:Ntau)/ pmax
!psiz(1:Ntau) = psiz(1:Ntau)/ pzmax
!
!open(unit=1,file='tfbtest.dat')
!do itau = 1,Ntau
!write(1,*) tau(itau),psi(itau)
!enddo
!close(1)
!write(*,*) 'run time', iend - istart, pmax
!
!
!!ier =pgopen('/Xserve')
!!call pgsvp(0.1,0.9,0.1,0.9)
!!call pgswin(0.0,tau(Ntau),0.0,1.0)
!!call pgbox( 'bcmst', 0., 10.0, 'bcnst', 0., 10. )
!
!!mean lags
!botpsi = 0.d0
!toppsi = 0.d0
!botpsiz = 0.d0
!toppsiz = 0.d0
!do itau = 1,Ntau
!toppsi = toppsi + psi(itau)*tau(itau)
!botpsi = botpsi + psi(itau)
!toppsiz = toppsiz + psiz(itau)*tau(itau)
!botpsiz = botpsiz + psiz(itau)
!enddo
!psimean = toppsi/botpsi
!psizmean = toppsiz/botpsiz
!
!
!!call pgline(Ntau,tau,psi)
!x2(1:2) = psimean
!y2(1) = 0.0
!y2(2) = 1.0
!!call pgline(2,x2,y2)
!
!!call pgsci(2)
!x2(1:2) = psizmean
!!call pgline(2,x2,y2)
!!call pgline(Ntau,tau,psiz)
!
!
!!
!!call pgend
!
!
!write(*,*) 'reshift',z
!write(*,*) 'taurest',psimean
!write(*,*) 'tauobs',psizmean
!write(*,*) 'tauobs/taurest',psizmean/psimean
!write(*,*) '(1+z)^[1+1/alpha]',(1+z)**(1.-1./slopevisc)
!write(*,*) '(mdotobs/mdotrest)**1/3', (uemdot/uemdotref)**(-1./3)
!write(*,*) '(mdotobs/mdotrest)**1/(4*alpha)', (uemdot/uemdotref)**(-1./(4*slopevisc))
!write(*,*) ''
!
!
!end program



!fortran subroutine to calculate the response function for an accretion disc
! assuming SEPARATE	viscous and irradiating temperature laws
!should also be computationally faster and neater

!gfortran tfb_x.f90 ../myrs2ld.f90 ../ran_new2.f90 ../mytr.f90


!code currently assumes tiny lamppost height that doesn't affect time lags
!can only calculate response functions at reference radius of 1 ld
!there is a multiplicative correction of 1.07 between T0v here and T0v in tfb_redshift
!(this one is 1.07 x less than tfb_redshift but I think this one is correct)


subroutine tfbx(tauz,Ntau,wav_prez,T1vin,T1iin,slope_v,slope_i,embh,deginc,rinsch,psi,z)
integer,intent(in):: Ntau
real, intent(in):: tauz(Ntau), T1iin, T1vin, slope_v, slope_i, embh,wav_prez,deginc,z
real, intent(out):: psi(Ntau)
real tau(Ntau)
logical info_test /.False./
!h=6.62607004e-34
!k=1.38064852e-23
!c = 2.99792458e8
!hc/k = 0.014387773538277204 !if wav in angstroms then 1.4387773538277204e8
!how to deal with any negative time lag bins


!organize redshift dependent quantities
zp1 = 1. + z
wav = wav_prez/zp1
tau(1:Ntau) = tauz(1:Ntau)/zp1


ibzero = 0
itau = 1
taunow = tau(itau)
do while (taunow .lt. 0)
taunow = tau(itau)
ibzero = ibzero + 1
itau = itau + 1
enddo

!write(*,*) ibzero,taunow,'checking tau ib0'

x_height = -3.0



if (info_test) then
write(*,*) 'input parms'
write(*,*) 'wav,T1vin,T1iin,slope_v,slope_i,embh,deginc,rinsch'
write(*,*) wav,T1vin,T1iin,slope_v,slope_i,embh,deginc,rinsch
write(*,*)
endif
pi = 3.141592653589793238462643383

psi(1:Ntau) = 0
twopi = 2*pi
xcon = 1.4387773538277204e8/wav
radinc = pi/180*deginc
sininc = sin(radinc)
cosinc = cos(radinc)

rinld = 3*f_rs2ld(rinsch,embh)


r0 = 1.0
if (rinld .gt. r0) then
write(*,*) 'inner ISCO bigger than reference radius!'
write(*,*) 'raising reference radius',rinld, r0
r0 = 10*rinld
endif

!convert mdot to temperature if needed
if (T1vin .lt. 0) then
eta = 0.1
newmmdot = 1
call tr4visc(r0,embh,T1vin,rinld, eta, newmmdot, tv40)
T1v = sqrt(sqrt(tv40))!*1.07
tv40 = tv40!*1.07**4
if (info_test) write(*,*) 'tfx convert mdot to tv',T1vin,T1v
else
T1v = T1vin
TV40 = T1v**4
endif

if (T1iin .lt. 0) then
eta = 0.1
newmmdot = 1
cosang = -666.0
call tr4irad(r0,embh,T1iin, x_height, cosang, eta, newmmdot, ti40)
T1i = sqrt(sqrt(ti40))
if (info_test) write(*,*) 'tfx convert mdot to ti',T1iin,T1i
else
T1i = T1iin
Ti40 = T1i**4
endif


!define a constant to prevent floating point over/underflows
overflow = TV40
TV40 = TV40/overflow
Ti40 = Ti40/overflow
overflow_fourth = overflow**0.25
!
xhld = abs(x_height)/rinsch*rinld!rinld



dtau = (tau(Ntau) - tau(1))/(Ntau-1)



r = rinld
idxr = 1
idxtot = 0
naz_base = 200
!daz = twopi/(naz-1)
slope_i_4 = 4*slope_i
slope_v_4 = 4*slope_v
dlnr = 0.004


!bin smoothing parameters
sig_gaus = dtau!/2!0.1!2.5*dtau
sign = 5.0!3.0
tau_smooth = sig_gaus*sign
bin_smooth = int(floor(tau_smooth/dtau))

X = 0



taumax = tau(Ntau)
!!save output to a text file
if (info_test) then
open (unit = 1, file = 'tfx_testtr.dat')
write(1,*) 'temperatures below divided by TV04**0.25 need to multiply by...'
write(1,*) overflow_fourth
write(*,*) 'doing tfx now'
call system_clock(istart)
endif

rmax = min(taumax/(1.-sininc),taumax*10)
r0rin_const = (1. - sqrt(rinld/max(r0,rinld)))
do while ((X .lt. 18 .and. r .lt. rmax) .or. rprev .lt. 1 )

!update radius logarithmically up to dr = dtau then linearly
rprev = r

if (dr .lt. dtau) then
r     = rinld*exp(idxr*dlnr)
dr    = r - rprev
else
r     = rprev + dtau
dr = dtau
endif


!solid angle
naz_new = min(naz_base,int(floor(twopi*r/dr) + 1))
daz = twopi/(naz_new-1)
dsa = r*dr*daz

!azimuth grid

do ia = 1,naz_new-1
azlo = (ia-1)*daz
azhi = azlo + daz
aznow = ranu(azlo,azhi,iseed)
cosaz = cos(aznow)
rnow = rprev + ranu(0.0,dr,iseed)

!time lag
taulag = rnow*(1.+ sininc*cosaz)

!taulag = (
!dd=(xhld*xhld+rnow*rnow)
!d=sqrt(dd)
!!write(*,*) 'xhld,x_hieght',xhld,x_height
!taulag=( -xhld*cosinc - r*cosaz*sininc + d )

!temperature (can repeat this for every azimuth but random jitter in r only improves
!response function smoothness when recomputing the time lag)
if (ia == 1) then

rnr0   = r0/rnow
rnr0_sv4 = rnr0**slope_v_4
rnr0_si4 = rnr0**slope_i_4
T4visc = Tv40*rnr0_sv4!*(1. - sqrt(rinld/max(rnow,rinld)))/r0rin_const
T4irad = Ti40*rnr0_si4
T4_tot = T4visc + T4irad
T2_tot = sqrt(T4_tot)
T_tot  = sqrt(T2_tot)
X = xcon/T_tot/overflow_fourth !overflow stops floating point overflows and nans for large TV40 temperatures above are divided by overflow_fourth


!response
X2 = X*X
X3 = X2*X
X5 = X3*X2
eX = exp(X)
eX_1 = eX-1.

!floating overflow checks
if (X .gt. 85) then
resp = 0
elseif (X .gt. 50) then
resp = X5/eX * rnr0_si4 * dsa
else
resp = X5*eX/(eX_1*eX_1) * rnr0_si4 * dsa
endif

if (resp .ne. resp) then
write(*,*) 'resp nan',resp,X**5,eX,X,rnr0_si4, dsa
!read(*,*)
exit
endif

!if diagnosing problem save info to output text file
if (info_test) then
write(1,*) rnow, T4visc**0.25, T4irad**0.25
endif


endif


!smooth over several time bins
sum  = 0.d0
iblo = int(max(1.0,ibzero + floor(taulag/dtau)-bin_smooth))
ibhi = int(min(Ntau*1.,ibzero + floor(taulag/dtau)+bin_smooth))
do itau = iblo,ibhi
taubin = tau(itau)
a = (taulag - taubin)/sig_gaus
a2 = a*a/2
weight = exp(-a2)
sum = sum + weight
psi(itau) = psi(itau) + weight*resp
enddo

!write(*,*) tau(iblo:ibhi)
!write(*,*) psi(iblo:ibhi)
!write(*,*) X,weight,resp,X,X5,eX_1,dsa,rnr0_si4!a,-a2,sig_gaus,taulag,taubin
!write(*,*)

!smoothing normalisation doesn't seem to do anything... ignoring it.
!do itau = iblo,ibhi
!psi(itau) = psi(itau)/sum
!enddo



idxtot = idxtot + 1

enddo



idxr = idxr + 1


if (info_test) then
call system_clock(iend)
itime = iend - istart
if (itime .gt. 1.e4) then
write(*,*) 'tf routine taking too long'
write(*,*) 'X,r,rprev,rnow,T_tot',X,r,rprev,rnow,T_tot
endif
endif

enddo



if (info_test) then
write(*,*) 'endof tfx',rmax,X,rprev
write(*,*) psi(1),psi(10),psi(Ntau-10:Ntau)
close(1)
endif





psi(1) = 0.0
psimax = maxval(psi)
do it = 1,Ntau
psi(it) = psi(it)/psimax
enddo

!write(*,*) 'tfbxop'
!write(*,*) tau(1:30)
!write(*,*)
!write(*,*) psi(1:30)
!do the interpolation if necessary
!if (z > 0.001) then
! call itp1d(psiop,ichoplo,ichophi,tau,tauz,psi,Ntau,Ntau,0)
!else
! psiop(1:Ntau) = psi(1:Ntau)
!endif
! how to deal with any negative time lag bins
if (ibzero .gt. 1) psi(1:ibzero-1) = 0

end subroutine





!end program
!



!function to convert from schwarschild radii to light days and back

function f_rs2ld(rin,embh)
real,intent(in):: rin, embh
double precision const, rs1, dl1, a
real myrs2ld
!km
const = 2.9533622197667164
rs1 = const*embh
dl1 = 2.590206837e10

a = rin*rs1/dl1


f_rs2ld = real(a)


end function



function f_ld2rs(rin,embh)
real,intent(in):: rin, embh
double precision const, rs1, dl1, a

!km
const = 2.9533622197667164
rs1 = const*embh
dl1 = 2.590206837e10

a = rin*dl1/rs1
f_ld2rs = real(a)
end function



!----------------------------------------
function dl_zoo( z, om, ol )
! luminosity diameter distance in units of c/H_0
! 2004 Mar Keith Horne @ St-And
dl_zoo = dp_zoo( z, om, ol ) * ( 1. + z )
return
end
!----------------------------------------
function da_zoo( z, om, ol )
! angular diameter distance in units of c/H_0
! 2004 Mar Keith Horne @ St-And
da_zoo = dp_zoo( z, om, ol ) / ( 1. + z )
return
end
!---------------------------------------
function dp_zoo( z, om, ol )
! proper distance (with curvature) in units of c/H_0
! 2004 Mar Keith Horne @ St-And
! 2005 May KDH @ St-And --
d = d_zoo( z, om, ol )
o = om + ol
if( o .gt. 1. ) then
r0 = 1. / sqrt( o - 1. )
dp = r0 * sin( d / r0 )
else if( o .lt. 1. ) then
r0 = 1. / sqrt( 1. - o )
dp = r0 * sinh( d / r0 )
else
dp = d
end if
dp_zoo = dp
return
end
!---------------------------------------
function d_zoo( z, om, ol )
! co-moving coordinate distance in units of c/H_0
! 2005 Dec KDH @ St-And
d = cosmic( 0., z, 0., -1., 0., om, ol )
d_zoo = d
return
end
!---------------------------------------
function t_oo( om, ol )
! age in units of 1/H_0
! 2004 Feb Keith Horne @ St-And
zinf = 1e3
t_oo = t_zzooo( 0., zinf, 0., om, ol )
return
end
!---------------------------------------
function t_zoo( z, om, ol )
! look-back time in units of 1/H_0
! 2004 Feb Keith Horne @ St-And
t_zoo = t_zzooo( 0., z, 0., om, ol )
return
end
!---------------------------------------
function t_zzooo( z1, z2, or, om, ol )
! elapsed time in units of 1/H_0
! 2004 Jul Keith Horne @ St-And
t_zzooo = cosmic( z1, z2, -1., -1., or, om, ol )
return
end
!---------------------------------------
function h_zooo( z, or, om, ol )
! Hubble parameter in units of H_0
! 2004 Jul Keith Horne @ St-And
x = 1. + z
x2 = x * x
ok = 1. - (or + om + ol)
hh = ol + x2 * ( ok + x * (om + x * or ) )
h_zooo = sqrt( hh )
return
end
!---------------------------------------
function cosmic( zlo, zhi, powx, powh, or, om, ol )
! integrate (1+z) ** powx * H(z) ** powh * dz
! 2004 Feb Keith Horne @ St-And
! 2004 Jul KDH @ St-And - implicit real*8
implicit real*8 (a-h,o-z)
real*4 zlo, zhi, powx, powh, or, om, ol,cosmic
logical logint/.false./

xlo = 1. + zlo
xhi = 1. + zhi

if( logint ) then
xloglo = dlog( xlo )
xloghi = dlog( xhi )
powxp = powx + 1.
end if

powh2 = powh / 2.
ok = 1.d0 - ( or + om + ol )

n = abs( zhi - zlo ) * 1000
n = max( 10, min( 10000, n ) )

sum = 0.d0
do i = 1, n
part = ( i - 0.5 ) / n
if( logint ) then
xlog = xloglo * ( 1. - part ) + xloghi * part
x = exp( xlog )
else
x = xlo * ( 1. - part ) + xhi * part
end if

h2 = ol + x*x * ( ok + x * ( om + x * or)  )

add = 1.d0
if( logint ) then
if( powxp .ne. 0. ) add = x ** powxp
if( powh2 .ne. 0. ) add = add * h2 ** powh2
else
if( powx .ne. 0. ) add = x ** powx
if( powh2 .ne. 0. ) add = add * h2 ** powh2
end if
sum = sum + add
end do

if( logint ) then
dlogx = ( xloghi - xloglo ) / n
cosmic = dlogx * sum
else
dx = ( xhi - xlo ) / n
cosmic = dx * sum
end if
return
end



!code eto calculate viscous and irradiating tr law for a give M and Mdot



!! calculates viscous contribution to tr law TO 4TH POWER
!input r, rin, embh, edrat (ld,ld,emsun,dimensionless)
!output tempvisc (T^4) viscous temperature to 4th power
subroutine tr4visc(r_inp,embh,edratin,rin_inp, eta, newmmdot, tv4)
integer, intent(in):: newmmdot
real, intent(in):: r_inp, embh, edratin, eta,rin_inp
real, intent(out):: tv4
logical fcv/.true./
double precision,save:: constvisc,rs2dl
double precision gnewt, sigma, c, pi, rldmet, unitconv,dlkm,rskm
real:: emdot,emsunyr

if (fcv .or. newmmdot == 1) then
gnewt  = 6.67408e-11
sigma  = 5.67036713e-8
pi     = 3.14159265358979323846264338
c      = 2.99792458e8
emsun   = 1.98855e30
emsunyr = 6.30321217e22
rldmet = 2.59020683712e13
unitconv = emsun/rldmet/rldmet*emsunyr/rldmet
if (edratin .gt. 0) then
edrat = edratin
call myedrat2mdot(embh,edrat,eta,emdot)
else
emdot = abs(edratin)
endif
constvisc = 3*gnewt/sigma/8/pi*embh*emdot * unitconv
dlkm = 2.5902068e10
rskm = 3.*embh
rs2dl = rskm/dlkm
!write(*,*) 'mytr.f90 viscous firstcall constant',rin, embh,emdot
!stop
end if

if (r_inp .lt. 0) then
r = abs(r_inp) * rs2dl
else
r = r_inp
endif

if (rin_inp .lt. 0) then
rin = abs(rin_inp) * rs2dl
else
rin = rin_inp
endif


r2 = r*r
r3 = r2*r
tv4_a = (1. - sqrt(rin/max(r,rin))) / r3
tv4 = constvisc*tv4_a

fcv = .false.

end subroutine



!lamppost contribution to temperature radius law (does not require flat accretion disc
!need to work out incidence cosine angle between current surface element and lamppost (cosang)
!used as input to subroutine
!input  r,xh, radius and lamppost height (light days)
!       embh, edrat black hole mass(M_0) and eddington ratio
!       cosang incidence cosine angle between surface element and lamppost
!       if cosang = -666.0 assume flat disk and calculate based on lamppost height
!       eta black hole efficiency parameter
!       newmmdot: if 0 then dont recalculate the eddington ratio or constirad. Do if set to 1
!output ti4 irradiating contribution to tr law raised to 4th power

subroutine tr4irad(r_inp,embh,edratin, xh_inp, cosang, eta, newmmdot, ti4)
integer,intent(in):: newmmdot
real, intent(in):: r_inp,embh,edratin, xh_inp, eta, cosang
real, intent(out):: ti4
logical fci/.true./
double precision,save:: constirad
real,save:: emdot,rs2dl
double precision gnewt, sigma, c, pi, unitconv, emsunyr, rldmet, ti4_a,&
d2,d,dlkm,rskm


if (fci .or. newmmdot == 1) then
sigma  = 5.67036713e-8
pi     = 3.14159265358979323846264338
c      = 2.99792458e8
c2     = c*c
emsunyr = 6.30321217e22
rldmet = 2.59020683712e13
if (edratin .gt. 0) then
edrat = edratin
call myedrat2mdot(embh,edrat,eta,emdot)
else
emdot = abs(edratin)
endif
unitconv  = emsunyr/rldmet/rldmet
constirad = 1./8/pi*c2/sigma*eta*unitconv
dlkm = 2.5902068e10
rskm = 3.*embh
rs2dl = rskm/dlkm
!write(*,*) 'mytr.f90 iradiation firstcall constant', constirad,emdot
end if

if (r_inp .lt. 0) then
r = abs(r_inp) * rs2dl
else
r = r_inp
endif

if (xh_inp .lt. 0) then
xh = abs(xh_inp) * rs2dl
else
xh = xh_inp
endif

if (r < xh) then
ti4 = 0
else
d2 = r*r + xh*xh
!not a flat disk
if (cosang .ne. -666.0) then
ti4_a = emdot/d2 * cosang
else
d = sqrt(d2)
ti4_a = emdot/d2/d*xh
endif

ti4 = ti4_a * constirad
endif

fci = .false.
end subroutine






!subroutine to calculate emdot for given eddington ratio
subroutine myedrat2mdot(um,er,eta,umdot)
real,intent(in):: um,er,eta
real,intent(out):: umdot
umdot = er/460.7e6 / eta * um
end subroutine




subroutine mmdot2edrat(um,umdot,eta,er)
real,intent(in):: um,umdot,eta
real,intent(out):: er

er = umdot*460.7e6*eta/um


end subroutine













!!! Program to take in an input time series, identify the brightewst data point and calculate a luminosity distance
!! according to equation 15 of cackett et al 2007.


!! inputs
!taucent: centroid delay (days)
!wav: (central wavelength (ANGSTROMS))
!cosi: cosine of inclination
!fvb: fnu bright state (mJy)
!fvd: fnu faint state (I assume this and the above mean the max and min data point in mJy)

!! OUTPUT luminosity distance (Mpc)
function dlflat(taucent,wav,cosi,fvb,fvf)

real taucent,wav,cosi,fvb,fvf,dlflat

epsilon=fvf/fvb
b=(1.-epsilon)/(1.-epsilon**1.5)

dlflat = 6.3 * taucent * (wav/1.e4)**(-1.5) * 1/sqrt(fvb/cosi) * b*b

return
end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Same as above but accepts mmdot rather than taucent as input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dlflatm(ummdot,wav,cosi,fvb,fvf)

real taucent,wav,cosi,fvb,fvf,dlflatm,ummdot

epsilon=fvf/fvb
b=(1.-epsilon)/(1.-epsilon**1.5)
ulam0=1000.0


call taumean(tm,ummdot,ulam0,wav)


dlflatm = 6.3 * tm * (wav/1.e4)**(-1.5) * 1/sqrt(fvb/cosi) * b*b

return
end function



!program trtest
!real a


!do i=1,2

!call t_r(t_rout,0.1,1.0e7,1.0,1.*i)
!write(*,*) i,t_rout

!enddo



!! work out mean delay
!write(*,*) 'time delay mean info'

!ummdot=1.e7
!ulam0=1000.0
!ulam=1367.2


!call taumean(tm,ummdot,ulam0,ulam)

!write(*,*) tm,ummdot,ulam0,ulam

!end program


!!! temperature radius profile of flart accretion disk
!!! ur,um,umdot, input radius in rs, um, um dot (mass and accretion rate [M0 and M0/yr])
!!! eta, efficiency parameter 0 to 1 (includes information due to disk albedo)

subroutine t_r(t_rout,eta,um,umdot,ur,uz)
real eta,um,umdot,ur
real t_rout
logical output
output = .False.

! physics constants
pi = 4. * atan2( 1., 1. )
twopi = 2. * pi
fourpi = 4. * pi
c = cgs( 'C' )
bk = cgs( 'K' )
gnewt = cgs( 'G' )
h = cgs( 'H' )
s = cgs( 'SIGMA' )
! astro constants
emsun = cgs( 'MSUN' )
day = 24. * 3600.
yr = cgs( 'YR' )
angst = cgs( 'ANGSTROM' )
urin=1.0
rg=gnewt*um*emsun/c/c

uring=urin/2 !in terms of grav rad
uzg=uz/2
urg=ur/2
dotm=umdot*emsun/yr

bot=10**(alog10(8*pi)+alog10(s)+2*alog10(rg))

T04=3*dotm/bot*c*c


T4=T04/urg**3 * (1.-sqrt(uring/urg) + (eta*uzg)/(sqrt(1.+(uzg/urg)**2)))

t_rout=sqrt(sqrt(T4))

if (output) then
write(*,*) 'uring',uring
write(*,*) 'rg',rg
write(*,*) 'urg',urg
write(*,*) 'ur',ur
write(*,*) 'T04',T04
endif


return
end subroutine




!!! function tau0 work out the reference tau according to horne et al in prep (tau0)
!! input ummdot mass x acc rate (M0**2 yr**-1)
!! input wav0 wavelength (Angstroms)
!! Output reference delay at wavelength wav0

function tau0(ummdot,uwav0)
real tau0,ummdot,wav0,uwav0

emsun=cgs( 'MSUN' )
year = cgs( 'YR' )
hplanck=cgs( 'H' )
c=cgs( 'C' )
G=cgs('G')
pi= 4. * atan2( 1., 1. )

wav0=uwav0*1e-8
dotmm=ummdot!*emsun/year*emsun
tau0=((45*G*dotmm*wav0**4)/(16*pi**6*hplanck))**(1./3)/c**(5./3)
tau0=tau0*emsun**(1./3)/year**(1./3)
tau0=tau0*emsun**(1./3)

!! put on the emsun**2/year at end to prevent overflow

return
end function








!!! routine to calculate the mean delay for face on inclination
!input : ummdot massxacc rate (Msol**2yr)
! ulam0: reference wavelength (ANG)
! ulam: wavelength (ANG)
!!output tm: mean delay (days)

subroutine taumean(tm,ummdot,ulam0,ulam)
real tau0
double precision x,dx
logical output /.True./
!! b is 3/4 for flat disk
b=3./4

!physics constants
pi = 4. * atan2( 1., 1. )
twopi = 2. * pi
fourpi = 4. * pi
c = cgs( 'C' )
bk = cgs( 'K' )
gnewt = cgs( 'G' )
h = cgs( 'H' )
s = cgs( 'SIGMA' )
!astro constants
emsun = cgs( 'MSUN' )
day = 24. * 3600.
yr = cgs( 'YR' )
angst = cgs( 'ANGSTROM' )
day=24.*3600


urin=1.0


tauref=tau0(ummdot,ulam0)



dx=0.1
xlo=0.0
xhi=15.0

nx=ceiling((xhi-xlo)/dx)

top=0.d0
bot=0.d0

!! compute integral
do i=1,nx

x=xlo+i*dx
x2=x*x
x4=x2*x2
X4b=X4/x**(1/b)
csh=cosh(x)-1.0

top=top+x4 /csh !*dx
bot=bot+x4b/csh !/dx

!write(*,*) i,x,xlo,xhi,dx,nx,top,bot,csh,x4b,x4,dx/csh,cosh(x)
!read(*,*)
enddo



tm=tauref*(ulam/ulam0)**(1/b) *top/bot /day
!write(*,*) tauref,'tauref',ummdot,ulam0!top/bot*(ulam/ulam0)**(1/b)
if (output) then
write(*,*) 'taumean routine summary'
write(*,*) 'ulam (Ang):', ulam
write(*,*) 'ummdot (M0^2 /yr):, ummdot'
write(*,*) 'tau mean (days)', tm
endif



return


end subroutine




subroutine plotxy1(x,y,sig,NT)

integer ier,pgbeg,NT
real x(NT),y(NT),sig(NT)

ymax=maxval(y)
ymin=minval(y)
yextra=(ymax-ymin)/10

xmax=maxval(x)
xmin=minval(x)
xextra=(xmax-xmin)/10


!ier = pgbeg(0,'/Xserve',1,1)
!call pgenv(xmin,xmax,ymin-yextra,ymax+yextra,0,1)

!!call pgerrb(6,x,y,sig,NT)
!call pgline(NT,x,y)

!call pgend

end subroutine






subroutine plotxy2(x,y1,y2,sig1,sig2,NT)

integer ier,pgbeg,NT
real x(NT),y1(NT),y2(NT),sig1(NT),sig2(NT)

ymax=max(maxval(y1),maxval(y2))
ymin=min(minval(y1),minval(y2))
yextra=(ymax-ymin)/10

xmax=maxval(x)
xmin=minval(x)
xextra=(xmax-xmin)/10

!ier = pgbeg(0,'/Xserve',1,1)
!call pgenv(xmin,xmax,ymin-yextra,ymax+yextra,0,1)

!call pgerrb(6,x,y1,sig1,NT)
!!call pgline(NT,x,y1)
!call pgpt(NT,x,y1,9)

!call pgsci(2)
!call pgerrb(6,x,y2,sig2,NT)
!call pgline(NT,x,y2)
!call pgsci(1)
!call pgend

end subroutine



subroutine plotxy3(x,y1,y2,y3,sig1,sig2,sig3,NT)

integer ier,pgbeg,NT
real x(NT),y1(NT),y2(NT),y3(NT),sig1(NT),sig2(NT),sig3(NT)

ymax=max(maxval(y1),maxval(y2))
ymin=min(minval(y1),minval(y2))
yextra=(ymax-ymin)/10

xmax=maxval(x)
xmin=minval(x)
xextra=(xmax-xmin)/10

!ier = pgbeg(0,'/Xserve',1,1)
!call pgenv(xmin,xmax,ymin-yextra,ymax+yextra,0,1)

!call pgerrb(6,x,y1,sig1,NT)
!call pgline(NT,x,y1)

!call pgsci(2)
!call pgerrb(6,x,y2,sig2,NT)
!call pgline(NT,x,y2)

!call pgsci(3)
!call pgerrb(6,x,y3,sig3,NT)
!call pgline(NT,x,y3)

!call pgsci(1)
!call pgend

end subroutine











!!! plot a graph with a set of points and a set of lines
! all arrays real
!! input xpt(NTpt),ypt(NTpt),sigpt(NTpt), xline(NTline),yline(NTline)

subroutine plotxyptline(xpt,ypt,xline,yline,sigpt,NTpt,NTline)

integer ier,pgbeg,NTpt,NTline
real xpt(NTpt),ypt(NTpt),sigpt(NTpt),xline(NTline),yline(NTline)

ymax=max(maxval(ypt),maxval(yline))
ymin=min(minval(ypt),minval(yline))
yextra=(ymax-ymin)/10

xmax=maxval(xpt)
xmin=minval(xpt)
xextra=(xmax-xmin)/10


!ier = pgbeg(0,'/Xserve',1,1)
!call pgenv(xmin,xmax,ymin-yextra,ymax+yextra,0,1)

!call pgline(NTline,xline,yline)
!call pgpt(NTpt,xpt,ypt,9)

!call pgend

end subroutine



!! Need some subroutines to go from model to dust corrected result
!0) Bluen the filter centroid wavelengths to the emitted spectrum
!1) convert model f_nu to AGN dust corrected f_nu
!2) Redden the new f_nu
!3) Apply MW dust to reddened f_nu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Calculate the central wavelengths in the emitted fram
!!INPUTS wavobs(ANG), H0(kms-1MPC-1), DL(Mpc)
!! OUTPUTS wavemit (Ang)

function wavemit(wavobs,H0,Dl)
real wavemit,wavobs,H0,DL,cgs

c=cgs('C')


ratio=1.+H0*Dl/c

wavemit=wavobs/ratio

end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!! FUNCTION to apply dust to light curve in emitted frame
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUTS f_mod,wav,ebmv unaltered model flux (mJy),wavelength

! OUTPUT dustAGN (remaining flux after dust (mJy))

function dustAGN(f_mod,wavemit,ebmv)
real f_mod,wav,ebmv,dustAGN,extmag_agn

A=extmag_agn( wavemit, ebmv )

dustAGN=f_mod*10**(-0.4*A)

!write(*,*) 'end of agn dust',f_mod,A,f_mod*10**(-0.4*A),ebmv,wavemit

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!! NOTE 1st May wave and wav0 seem to be useless inputs so long as one remembers that fnumod
!! must be at the emitted wavelength, and the outputted flux is at the observed wavelength

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Redden the AGN dust corrected flux

!!!! INPUT wave (emitted wavelength: Angstroms)
!!!!       fe   (emitted f_nu (mJy)
!!!!       wav0 (observed wavelength (ANG))

!!!! OUTPUT fnuetofnu0 (observed flux mJy)

function fnuetofnu0(wave,fe,wav0,H0,DL)
real fnuetofnu0,fe,wave,wav0,H0

c=2.9979e5 !! c in kms-1

zplus1=1.+c*DL/H0


fnuetofnu0=fe/zplus1*zplus1

end function




!! APPLY MILKY WAY DUST TO THE LIGHTCURVE IN THE OBSERVED FRAME (AT THE OBSERVED WAV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUTS f_obs(f_obs before reddening mJy),wav (ANG),ebmv unaltered model flux (mJy),wavelength

! OUTPUT dustAGN (f_obs ofter reddening mJY)

function dustMW(f_obs,wavobs,ebmv)
real f_obs,wavobs,ebmv,dustMW,extmag_mw

A=extmag_mw( wavobs, ebmv )

dustMW=f_obs*10**(-0.4*A)

!write(*,*) 'end of mw dust',f_obs,A,f_obs*10**(-0.4*A),ebmv,wavobs

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!! APPLY MILKY WAY DUST TO THE LIGHTCURVE IN THE OBSERVED FRAME (AT THE OBSERVED WAV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUTS f_obs, sig_obs (f_obs before reddening mJy),wav (ANG),ebmv unaltered model flux (mJy),wavelength

! OUTPUT fop,sigop (f_obs ofter reddening mJY)

subroutine dustMW_sig(f_obs,sig_obs,wavobs,ebmv,fop,sigop)
real f_obs,wavobs,ebmv,dustMW,extmag_mw, sig_obs

A=extmag_mw( wavobs, ebmv )

b = 10**(-0.4*A)
fop   =f_obs*b
sigop =sig_obs*b


!write(*,*) 'end of mw dust',f_obs,A,f_obs*10**(-0.4*A),ebmv,wavobs

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! i.e take dust away so output is brigther than input
!! DE-APPLY MILKY WAY DUST TO THE LIGHTCURVE IN THE OBSERVED FRAME (AT THE OBSERVED WAV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUTS f_obs, sig_obs (f_obs before reddening mJy),wav (ANG),ebmv unaltered model flux (mJy),wavelength

! OUTPUT fop,sigop (f_obs ofter reddening mJY)

subroutine de_dustMW_sig(f_obs,sig_obs,wavobs,ebmv,fop,sigop)
real f_obs,wavobs,ebmv,dustMW,extmag_mw, sig_obs

A=extmag_mw( wavobs, ebmv )

b = 10**(0.4*A)
fop   =f_obs*b
sigop =sig_obs*b


!write(*,*) 'end of mw dust',f_obs,A,f_obs*10**(-0.4*A),ebmv,wavobs

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


















!!! program to read the number of items in a file testread
!program testread2nline
!integer read2nline,read2nstring
!character(100) fname
!character(1) typ
!
!fname='testread.txt'
!typ='s'
!
!
!n=read2nline(fname,'i')
!
!write(*,*) n
!
!end program


!! function to read the number of lines in a file ignores blank spaces
!! inps character(100) fname
!		character(1) typ 's' for string or not 's' if reading a data file (i.e not using 's', the code will discount all lines beginning with a character i.e can comment things out by putting e.g !)
function read2nline(fname,typ)
character(*):: fname
character(100)::string
character(1):: typ
integer nline,read2nline,irecount

n=0
iend=0
!fname=adjustl(fname)


ifin=0


open(unit=666,file=trim(fname))
do while (ifin .lt.  1)

if (typ =='s') then
read(666,*, iostat=iend) string
string=adjustl(string)
lenstring=len(trim(string))
!write(*,*) lenstring, trim(string),iend

if (lenstring .gt. 0 .and. iend .eq. 0) then
n=n+1  !! this should only count non blank lines
else
ifin=1
cycle
endif

else if (typ .eq. 'r') then
read(666,*, iostat=iend) recount

!write(*,*) recount,iend,'sds',ifin,n
!read(*,*)
if (iend .eq. 0) then
n=n+1
else
ifin=1
cycle
endif


else if (typ .eq. 'i') then
read(666,*, iostat=iend) irecount

!write(*,*) recount,iend,'sds',ifin,n
!read(*,*)
if (iend .eq. 0) then
n=n+1
else
ifin=1
cycle
endif



end if
enddo

close(666)

read2nline=n



return
end function





!! shuffle an array fortran
!! input a(n) input array
!! output a(n) shuffled array

subroutine shuffle(a,n)

integer n
real a(N), ran3

call system_clock(iseed)


do i = n, 1, -1
j = i * ran3(iseed) + 1
t = a(i)
a(i) = a(j)
a(j) = t
end do

end subroutine



! tested 24/10/2014

!program testshuffle
!
!real x(1000), neach(1000)
!
!do i =1,1000
!x(i) = i
!neach(i) = 0
!enddo
!
!
!call shuffle(x,1000)
!nunshift = 0
!
!do i =1,1000
!neach(int(x(i))) = neach(int(x(i)))+1
!enddo
!
!
!do i =1,1000
!!write(*,*) i, x(i), neach(i)
!if (real(i) .eq. x(i)) then
!!write(*,*) 'bollocks'
!nunshuft = nunshift + 1
!endif
!
!if (neach(i) .gt. 1) then
!write(*,*) 'crap'
!stop
!endif
!
!enddo
!
!write(*,*) 'number of unshifted points', nunshift
!end program
!




!!! routine to ingest an input transfer function with an updated either mmdot or wavelength
!.. then scale the x and y axis according to the scaling laws in the book with the red strap

!INPUTS:
!		tau(Ntau),psiold(Ntau), parold, parnew, ipar (this is 1 for mmdot, or 2 for wavelength)
!		iap = -6 then we scale for t0v in T = T0 (r/r0)^-3/4  for flat disk T0 propto mmdot^-1/4 and so we modify mmdot scaling relation for t0 e.g 1/3 -->-4/3, 11/12 --> -11/3, 1/4--> -1


!OUTPUT:
!		parnew(Ntau) the new psi

subroutine tfupdate(Ntau, tau, psiold, psinew, parold, parnew, ipar)

integer,intent(in):: Ntau, ipar
real,intent(in):: psiold(Ntau), tau(Ntau)
real,intent(out):: psinew(Ntau)
integer choplo, chophi
real psimid(Ntau), taunew(Ntau)



!!!!!!!! Different x and y axis scalings for wavelength and mmdot changes
cosh1sub1 = cosh(1.) - 1.
parrat = parnew/parold

if (ipar .eq. 1) then !mmdot !0.91667 = 11/12, 0.333 = 1./3

yscale = parrat**(0.91667) * cosh1sub1 / (cosh(parrat**0.25) - 1.)
xscale = parrat**(0.333)


else if (ipar .eq. -6) then !all powers for -6 should be -ve but for some reason, +ve more accurate

yscale = parrat**(3.66666666) * cosh1sub1 / (cosh(parrat**(1.0)) - 1.)
xscale = parrat**(1.3333333333)


else				  !wavelength

yscale = parrat**(-4.)* cosh1sub1/(cosh(1./parrat) - 1.)
xscale = parrat**(1.33333)

endif
!!!!!!!!

!!!!!!! create the new psi and tau grids
taulo =  tau(1)
dtau  = ( tau(Ntau) - tau(1) )/(Ntau - 1)

xa = dtau*xscale
do itau = 1,Ntau
taunew(itau) = taulo + (itau - 1)*xa
psimid(itau) = psiold(itau) * yscale
enddo



! interpolate the new tau grid onto the old one
call itp1d(psinew,choplo,chophi,taunew,tau,psimid,Ntau,Ntau,0)

!psi new contains the rescaled transfer function


end subroutine



































!! New tf update routine to deal with a reference response function that is sampled with
! a different time resolution to the desired output response function. - 14/6/2017

!!! routine to ingest an input transfer function with an updated either mmdot or wavelength
!.. then scale the x and y axis according to the scaling laws in the book with the red strap

!INPUTS:
!		tau_hr(Ntau_hr), tau(Ntau), psiold(Ntau_hr), parold, parnew, ipar (this is 1 for mmdot, or 2 for wavelength)
!		iap = -6 then we scale for t0v in T = T0 (r/r0)^-3/4  for flat disk T0 propto mmdot^-1/4 and so we modify mmdot scaling relation for t0 e.g 1/3 -->-4/3, 11/12 --> -11/3, 1/4--> -1


!OUTPUT:
!		parnew(Ntau) the new psi


!call tfupdate_hr(Ntau_hr - idxtau0_hr+1, Ntau_hr-idxtau0_hr+1,&
!taugrid_hr(idxtau0_hr:Ntau_hr), taugrid_hr(idxtau0_hr:Ntau_hr),&
!psi_hr(idxtau0_hr:Ntau_hr), psi_hr(idxtau0_hr:Ntau_hr), umdotref, p(NPumbhidx), 1)

subroutine tfupdate_hr(Ntau_hr, Ntau, tau_hr, tau, psiold, psinew, parold, parnew, ipar)

integer,intent(in):: Ntau, ipar, Ntau_hr
real,intent(in):: psiold(Ntau_hr), tau(Ntau), tau_hr(Ntau_hr)
real,intent(out):: psinew(Ntau)
integer choplo, chophi
real psimid(Ntau_hr), taunew_hr(Ntau_hr)



!!!!!!!! Different x and y axis scalings for wavelength and mmdot changes
cosh1sub1 = cosh(1.) - 1.
parrat = parnew/parold

if (ipar .eq. 1) then !mmdot !0.91667 = 11/12, 0.333 = 1./3

yscale = parrat**(0.91667) * cosh1sub1 / (cosh(parrat**0.25) - 1.)
xscale = parrat**(0.333)


else if (ipar .eq. -6) then !all powers for -6 should be -ve but for some reason, +ve more accurate

yscale = parrat**(3.66666666) * cosh1sub1 / (cosh(parrat**(1.0)) - 1.)
xscale = parrat**(1.3333333333)


else				  !wavelength

yscale = parrat**(-4.)* cosh1sub1/(cosh(1./parrat) - 1.)
xscale = parrat**(1.33333)

endif
!!!!!!!!

!!!!!!! create the new psi and tau grids
write(*,*) 'here in tfupdate2'
taulo =  tau(1)
dtau  = ( tau(Ntau) - tau(1) )/(Ntau - 1)
dtau_hr  = ( tau_hr(Ntau_hr) - tau_hr(1) )/(Ntau_hr - 1)


xa_hr = dtau_hr*xscale
do itau = 1,Ntau_hr
taunew_hr(itau) = taulo + (itau - 1)*xa_hr
psimid(itau) = psiold(itau) * yscale
enddo


write(*,*) 'here 2 in tfupdate2',dtau,dtau_hr,Ntau,Ntau_hr,'before itp1d'
! interpolate the new tau grid onto the old one

!write(*,*) psinew(1), psinew(Ntau), Ntau, taunew(1), taunew(Ntau)
!write(*,*) psimid(1), psimid(Ntau_hr),Ntau_hr, taunew_hr(1), taunew_hr(Ntau_hr)
call itp1d(psinew,choplo,chophi,taunew_hr,tau,psimid,Ntau_hr,Ntau,0)

write(*,*) 'here 3 in tfupdate2'
!psi new contains the rescaled transfer function


end subroutine












!! program to test the rescaling routine gflags tfupdate.f90 ran_new2.for delaydist/tfb.f90 cgs.for itp2.for
! mmdot and wavelength scalings working. Now test the scaling interpolation scaling
!
!program tfupdatetest
!
!real, allocatable:: psi(:), tau(:), psinew(:), psilo(:), psihi(:)
!character(10):: tinc, tincp1
!character(100):: pwd
!
!
!
!apower = 3./4
!bpower = 3./4
!taulo = 0.0
!tauhi = 30.0
!dtau = 0.1
!urin = 3.0
!ur0  = 3.0
!embh = 1.e8
!emdot = 1.
!deginc = 0.0
!wav = 3000.0
!uzx = 3.0
!ak = 0.0
!alpha = 0.0
!iversion = 1
!eta  =0.1
!
!
!Ntau = (tauhi - taulo)/dtau + 1
!allocate(tau(Ntau), psi(Ntau))
!do itau = 1,Ntau
!tau(itau) = taulo + (itau - 1)*dtau
!enddo
!
!!! make the orriginal transfer function
!call tfb(tau,psi,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!
!
!open(unit = 1,file = 'testtfupdate0.op')
!do itau = 1,Ntau
!write(1,*) tau(itau), psi(itau)
!enddo
!close(1)
!
!!!
!
!
!
!!! now do the rescaling... tes for mmdot 2 an 3 times orriginal
!
!allocate(psinew(Ntau))
!
!embhnew = 3.*embh
!call tfupdate(Ntau, tau, psi, psinew, embh, embhnew, 1)
!call tfb(tau,psi,Ntau,wav,embhnew,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!
!a = maxval(psi)
!psi=psi/a
!a = maxval(psinew)
!psinew=psinew/a
!open(unit = 1,file = 'testtfupdate1.op')
!do itau = 1,Ntau
!write(1,*) tau(itau), psinew(itau), psi(itau)
!enddo
!close(1)
!
!
!embhnew = 2.*embh
!psinew(:) = 0
!psi(:) = 0
!call tfb(tau,psi,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!call tfupdate(Ntau, tau, psi, psinew, embh, embhnew, 1)
!call tfb(tau,psi,Ntau,wav,embhnew,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!open(unit = 1,file = 'testtfupdate2.op')
!do itau = 1,Ntau
!write(1,*) tau(itau), psinew(itau), psi(itau)
!enddo
!close(1)
!
!
!
!
!
!
!!!! now test for the wavelength dependence
!wavnew = 3*wav
!call tfb(tau,psi,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!call tfupdate(Ntau, tau, psi, psinew, wav, wavnew, 2)
!call tfb(tau,psi,Ntau,wavnew,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!
!a = maxval(psi)
!psi = psi/a
!a = maxval(psinew)
!psinew=psinew/a
!open(unit = 1,file = 'testtfupdate_w1.op')
!do itau = 1,Ntau
!write(1,*) tau(itau), psinew(itau), psi(itau)
!enddo
!close(1)
!
!
!
!
!wavnew = 2*wav
!!call system_clock(istart)
!call tfb(tau,psi,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!!call system_clock(iend)
!
!!call system_clock(istart2)
!call tfupdate(Ntau, tau, psi, psinew, wav, wavnew, 2)
!!call system_clock(iend2)
!
!!call system_clock(istart)
!call tfb(tau,psi,Ntau,wavnew,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!!call system_clock(iend)
!open(unit = 1,file = 'testtfupdate_w2.op')
!do itau = 1,Ntau
!write(1,*) tau(itau), psinew(itau), psi(itau)
!enddo
!close(1)
!
!
!!write(*,*) 'Time check tf, ssclae:', iend-istart, iend2-istart2
!
!
!
!
!
!
!
!!!1!!! test inclination scaling
!!deginc = 31.5
!!dstore = 5.0
!!
!!
!!
!!call tfb(tau,psilo,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!!call system_clock(istart)
!!call getcwd(pwd)
!!call chdir('./tfstore')
!!write(tinc, '(F10.2)') nint(deginc/dstore)*dstore
!!tinc = adjustl(tinc)
!!
!!write(tincp1, '(F10.2)') (nint(deginc/dstore) + 1.)*dstore
!!tincp1 = adjustl(tincp1)
!!
!!allocate(psilo(Ntau), psihi(Ntau))
!!
!!open(unit = 1, file = 'tftemp_'//adjustl(trim(tinc))//'.dat')
!!do itau = 1,Ntau
!! read(1,*) tau(itau), psilo(itau)
!!enddo
!!close(1)
!!
!!open(unit = 1, file = 'tftemp_'//adjustl(trim(tincp1))//'.dat')
!!do itau = 1,Ntau
!! read(1,*) tau(itau), psihi(itau)
!!enddo
!!close(1)
!!call system_clock(iend)
!!
!!yscale = (deginc - nint(deginc/dstore)*dstore)/dstore
!!write(*,*) 'yscale',yscale
!!do itau = 1,Ntau
!!psinew(itau) = psilo(itau) + (psihi(itau) - psilo(itau))*yscale
!!enddo
!!
!!call system_clock(istart1)
!!call tfb(tau,psi,Ntau,wav,embh,emdot,uzx,eta,urin,deginc,bpower,apower,ur0,ak,alpha,iversion)
!!call system_clock(iend1)
!!
!!
!!write(*,*) iend1 - istart1, 'evaluate again'
!!write(*,*) iend - istart, 'interpolate'
!!
!!call chdir(adjustl(trim(pwd)))
!!open(unit=1, file = 'testtfincscale.op')
!!do itau = 1,Ntau
!!write(1,*) tau(itau), psilo(itau), psihi(itau), psi(itau), psinew(itau)
!!enddo
!!close(1)
!!
!!
!!
!!
!!write(*,*) nint(deginc/dstore)
!
!
!
!
!
!
!
!
!
!end program
!
!!
!!
!!
!!
!!
!!
!!
!!




!! triangular shaped set of scatter plots showing the model in parameter space after Nits itertions
! input
!    n parms
!    Nits number of iterations
!    p(Nits,n) stored parameters


subroutine scatterparmsb(n,Nits,p,axlab)

integer n,Nits,pgopen,ier
real p(Nits,n),mean(n),sd(n),sdev,avg,xw0(2),yw0(2)
character*1000 plab(n),axlab(n),information,header
character*10 tinits,tsd(n),tmean(n)
integer ipcol(Nits,n)


ddeg=10.0
Ndeg=9
pi=atan2(1.,1.)*4
deg2rad=pi/180



!!! calculate the standard deviation of points
do i=1,n
mean(i)=avg(p(1:nits,i),nits)
write(tmean(i),*) mean(i)
tmean(i)=adjustl(tmean(i))


sd(i)=sdev(p(1:nits,i),nits)
write(tsd(i),*) sd(i)
tsd(i)=adjustl(tsd(i))

!sd4(i)=rms(p(1:nits,i),nits)
sd1=sd(i)
sd2=2*sd(i)
sd3=3*sd(i)
sd4=4*sd(i)
do it=1,Nits
ipcol(it,i)=1
if (abs(p(it,i)-mean(i)) .lt. sd4) ipcol(it,i)=2
if (abs(p(it,i)-mean(i)) .lt. sd3) ipcol(it,i)=3
if (abs(p(it,i)-mean(i)) .lt. sd2) ipcol(it,i)=4
if (abs(p(it,i)-mean(i)) .lt. sd1) ipcol(it,i)=5
enddo

enddo


!! trim the labels
do i=1,n
plab(i)=axlab(i)
plab(i)=adjustl(plab(i))
enddo



!ier=pgopen('parmparmplot.PS/CPS')!


!call pgbgw
svpup=0.95
svpdown=0.1
svpleft=0.1
svpright=0.95
!!call pgsvp(svpleft,svpright,svpdown,svpup)


svpdploty=(svpup-svpdown)/n
svpdplotx=(svpright-svpleft)/n




do ix=1,n

svpxright=svpleft+ix*svpdplotx
svpxleft=svpxright-svpdplotx




do iy=1,n
if (iy .lt. ix) cycle ! Do not mirror the plot




xmin=mean(ix)-4*sd(ix)!minval(p(1:nits,ix))
xmax=mean(ix)+4*sd(ix)!maxval(p(1:Nits,ix))
xrange=(xmax-xmin)/10
xminplot=xmin-xrange
xmaxplot=xmax+xrange

svpybot=svpup-svpdploty*iy
svpytop=svpybot+svpdploty

!!!! Include histograms
if (iy .eq. ix) then
Nbin=100
histmin=0.0
histmax=0.2*nits


!!!! Add the mean and standard deviations
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgbox( '', 0., 0, '', 0., 0 )
!call pgswin(xminplot,xmaxplot,histmin,histmax)
xw0(1:2)=mean(ix)
yw0(1)=histmin
yw0(2)=histmax
!call pgline(2,xw0,yw0)


!call pgsci(4)
!call pgsls(3)
!call pgsch(1.3)
xw0(1:2)=mean(ix)+sd(ix)
!call pgline(2,xw0,yw0)
xw0(1:2)=mean(ix)-sd(ix)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!call pgsch(1.0)

!!!!! Plot Histograms
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgswin(xminplot,xmaxplot,histmin,histmax)
!call PGHIST(nits, p(1:nits,ix), xminplot,xmaxplot , Nbin, 3)



header='Mean:'//trim(tmean(ix))//'    sd:'//trim(tsd(ix))
header=adjustl(header)
!call pgmtxt('T',-1.0,0.5,0.5,trim(header))


if (ix .eq. n) then
!call pglab(trim(plab(ix)),'','')
!call pgbox( 'bcstn', 0., 10, 'bcst', 0., 10 )
else
!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )
endif

write(*,*) ix,iy,svpxleft,svpxright,svpybot,svpytop,xminplot,xmaxplot

cycle
endif
!!!


!!!! make the scatter plot
ymin=mean(iy)-4.*sd(iy)!minval(p(1:Nits,iy))
ymax=mean(iy)+4.*sd(iy)!maxval(p(1:Nits,iy))
yrange=(ymax-ymin)/10
yminplot=ymin-yrange
ymaxplot=ymax+yrange

!write(*,*) 'here',ix,iy,svpxleft,svpxright,svpybot,svpytop
!write(*,*) 'here too',xminplot,xmaxplot,yminplot,ymaxplot
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgswin(xminplot,xmaxplot,yminplot,ymaxplot)
!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )

!!call pgenv(xminplot,xmaxplot,yminplot,ymaxplot,0,1)
do it=1,Nits
!!call pgsci(ipcol(it,
!call pgsci(min(ipcol(it,ix),ipcol(it,iy))) !! only if both parameters are in the confidence region does it pass to the higher confidece region
!call pgpt(1,p(it,ix),p(it,iy),-1)
enddo
!call pgsci(1)
!do it=1,n
!write(*,*) trim(plab(it))
!enddo
!!! decide whether to plot the axis labels

if ((ix.eq.1) .AND. (iy .ne. n)) then !! plot the axis labels on the left hand side
!call pglab('',trim(plab(iy)),'')
!call pgbox( 'bcst', 0., 10, 'bcstn', 0., 10 )

else if ((ix .eq. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),trim(plab(iy)),'')
!call pgbox( 'bcstn', 0., 10, 'bcstn', 0., 10 )

else if ((ix .ne. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),'','')
!call pgbox( 'bcstn', 0., 10, 'bcst', 0., 10 )
endif

read(*,*)



enddo ! end ix loop

enddo ! end iy loop


!!!!!! add on relevant information in the top right hand corner

svpxleft=svpxright-svpdplotx
svptop=svpup
svpbot=svpup-svpdploty

!call pgsvp(svpxleft,svpxright,svpbot,svptop)
!call pgswin(0.0,1.0,0.0,1.0)
write(tinits,'(I7.3)') nits
tinits=adjustl(tinits)

information='N='//trim(tinits)
information=adjustl(information)

!call pgtext(0.5,0.5,trim(information))

!call pgend

end subroutine













































































! scatter parms 2 (incorporates badness of fit information to produce confidence regions)

!! triangular shaped set of scatter plots showing the model in parameter space after Nits itertions
! input
!    n parms
!    isteplog(n) integer array stating whether parameter should be plotted as log (i.e for logarithmic steps) 0 for linear, 1 for log10
!    Nits number of iterations
!    p(Nits,n) stored parameters
!    bof(Nits)
!    preal(n) real array of actual input parameter values (for fake data) if set to p(:) = 0.0 do not plot



subroutine scatterparms(n,Nits,p,preal,axlab,bof,iburnin,isteplog)

integer n,Nits,pgopen,ier,isteplog(n)
real p(Nits,n),mean(n),sd(n),sdev,avg,bof(nits),xw0(2),yw0(2),preal(n)
character*1000 plab(n),axlab(n),information,header
integer ipcol(Nits),iburnin
logical info /.false./, posterplot /.true./
character*20 tinits,tsd(n),tmean(n)
character*10 tempincsd,tempincmean,tempinclo,tempinchi,tdegax
character*1000 tdeginc

pi=4.*atan2(1.,1.)
ddeg=10.0
Ndeg=6
deg2rad=pi/180



!!! character information tnits to add as header to plots
write(tinits,'(I7.3)') nits
tinits=adjustl(tinits)

!!! color of sd and mean lines
icmean=4
icsd=4
icreal=2
icoliburnin=2
!!! annotating character size
sizeann=1.65
sizedef=1.65

if (info) then
write(*,*) 'Scatterparms input parms...'
do i=1,Nits
write(*,*)i, bof(i), p(i,:)
enddo
endif


!!! calculate the mean of the points
do i=1,n
sum=0.d0
sd(i)=sdev(p(iburnin:nits,i),nits-iburnin+1)
write(tsd(i),'(F10.1)') sd(i)
tsd(i)=adjustl(tsd(i))

do it=iburnin,nits
sum=sum+p(it,i)
write(*,*) 'calculating mean',i,it,sum

enddo
mean(i)=sum/(nits-iburnin+1)
write(tmean(i),'(F20.1)') mean(i)
tmean(i)=adjustl(tmean(i))
write(*,*) trim(tmean(i)),mean(i)
read(*,*)
enddo

!!! calculate the standard deviation of points
bofmin=minval(bof)
dcs1=2.3
dcs2=6.17
dcs3=11.8
do it=1,Nits
ipcol(it)=1
if (abs(bof(it)-bofmin) .lt. dcs3) ipcol(it)=2
if (abs(bof(it)-bofmin) .lt. dcs2) ipcol(it)=3
if (abs(bof(it)-bofmin) .lt. dcs1) ipcol(it)=4
enddo
!do it=1,NIts
!write(*,*) bofmin,bof(it)
!enddo
!stop



!! trim the labels
do i=1,n
plab(i)=axlab(i)
plab(i)=adjustl(plab(i))

enddo


!ier=pgopen('parmparmplot.PS/CPS')!


!call pgslw(6)
!call pgsch(sizedef)
if (posterplot) then
!call pgsch(1.4)
!else
!!call pgsch(sizedef)
endif

!call pgsfs(2)

!call pgbgw
svpup=0.9
svpdown=0.15
svpleft=0.1
svpright=0.90
!!call pgsvp(svpleft,svpright,svpdown,svpup)


svpdploty=(svpup-svpdown)/n
svpdplotx=(svpright-svpleft)/n



do ix=1,n

svpxright=svpleft+ix*svpdplotx
svpxleft=svpxright-svpdplotx




do iy=1,n
if (iy .lt. ix) cycle ! Do not mirror the plot










xmin=mean(ix)-4*sd(ix)!minval(p(1:nits,ix))
xmax=mean(ix)+4*sd(ix)!maxval(p(1:Nits,ix))
xrange=(xmax-xmin)/10
xminplot=xmin-xrange
xmaxplot=xmax+xrange

svpybot=svpup-svpdploty*iy
svpytop=svpybot+svpdploty

!!!! Include histograms
if (iy .eq. ix) then
histmax=1.*(nits-iburnin+1)!/Nbin!0.5*iburnin
histmin=0.0
Nbin=100

!!!! Add the mean and standard deviations
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgbox( '', 0., 0, '', 0., 0 )
!call pgswin(xminplot,xmaxplot,histmin,histmax)



!!!! Add the histogram
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)

! The next if statements make it so the last histogram on the x axis is plotted rotated by 90 degrees)
if (ix .ne. n .or. n .eq. 1) then   !! if only one parameter varied or not on last histogram along
!call pgswin(xminplot,xmaxplot,histmin,histmax)
write(*,*) ix,'here'
!call PGHIST(nits-iburnin+1, p(iburnin:nits,ix), xminplot,xmaxplot , Nbin, 3)

if (n .ne. 1) then !! if more than one parameter varied and not on last one along

if (isteplog(ix) .eq. 0) then  !linear
!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstl', 0., 10, 'bcst', 0., 10 )
endif

else  !!if only one parameter varied

!! decide whether to plot log histograms (inputs must already be logged)
if (isteplog(ix) .eq. 0) then  !linear
!call pgbox( 'bcstn', 0., 10, 'bcst', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstnl', 0., 10, 'bcst', 0., 10 )
endif


!call pglab(trim(axlab(1)),'Number','')
endif

!! add the lines for mean and sd
xw0(1:2)=mean(ix)
yw0(1)=histmin
yw0(2)=histmax
!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
xw0(1:2)=mean(ix)+sd(ix)
!call pgline(2,xw0,yw0)
xw0(1:2)=mean(ix)-sd(ix)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)


















else if (ix .eq. n) then

if (isteplog(ix) == 3) then  !! force axis on histogram plot to be 0 to 1 if stepping in cos theta
iplotlimforce = 1


else
iplotlimforce = 0
endif

plotmin = 0.0
plotmax = 0.0
!call pghistalong(nits-iburnin+1, p(iburnin:nits,ix), xminplot, xmaxplot,plotmin,plotmax,dbin, Nbin, 2, 2.0,iplotlimforce) ! pgswin called inside this routine

!!call pglab(trim(plab(ix)),'','')

if (isteplog(ix) .ne. 1) then  !linear
!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcst', 0., 10, 'bcstl', 0., 10 )
endif


!! If we are plotting cos(inc) isteplog(i)=3. Then we add on the plot annotations to show the mean and lower and upper confidence limits in degrees
if (isteplog(ix) .eq. 3) then


!write(tempinclo,'(F10,2)') 180/pi*acos(mean(ix)-sd(ix))
write(tempincmean,'(F10.0)') 180/pi*acos(mean(ix))
tempincmean=adjustl(tempincmean)

!write(tempinchi,'(F10,2)') 180/pi*acos(mean(ix)+sd(ix))
write(tempincsd,'(F10.0)') abs( 180/pi*( acos(mean(ix)+sd(ix)) - acos(mean(ix)-sd(ix)) )/2 )
tempincsd=adjustl(tempincsd)

tdeginc=trim('i\uo\d = '//trim(tempincmean)//' \(2233) '//trim(tempincsd))
tdeginc=adjustl(tdeginc)
!call pgsch(sizeann)
!call pgmtxt('T',2.0,0.5,0.5,trim(tdeginc))
!call pgsch(sizedef)
endif
!!



!! add on an extra axis if cosine parameter (so that you can see the cosine in degrees)
do i=1,Ndeg
deg=ddeg*i
yw0(1:2)=cos(deg*deg2rad)

if ((yw0(1) .gt. xminplot) .and. (yw0(1) .lt. xmaxplot)) then
write(tdegax,'(I0)') int(deg)
tdegax=trim(trim(tdegax)//'\uo')
tdegax=adjustl(tdegax)
xw0(1)=plotmin+0.9*(plotmax-plotmin)
xw0(2)=plotmax
!call pgline(2,xw0,yw0)
!call pgmtxt('RV',2.0,(yw0(1)-xmaxplot)/(xminplot-xmaxplot),1.0,trim(tdegax))
endif

enddo
!!





!! add the lines for mean and sd
yw0(1:2)=mean(ix)
xw0(1)=plotmin
xw0(2)=plotmax
!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
yw0(1:2)=mean(ix)+sd(ix)
!call pgline(2,xw0,yw0)
yw0(1:2)=mean(ix)-sd(ix)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)








endif




header=trim(tmean(ix))//'\(2233)'//trim(tsd(ix))

if (isteplog(ix) .eq. 1) header='10\u'//trim(header) !! make the header into 10^ form

if (n .eq. 1) header=trim(header)//'    N:'//trim(tinits)
header=adjustl(header)
!call pgmtxt('T',0.2,0.5,0.5,trim(header))


write(*,*) ix,iy,svpxleft,svpxright,svpybot,svpytop,xminplot,xmaxplot,nits

cycle
endif
!!!


!!!! make the scatter plot
ymin=mean(iy)-4.*sd(iy)!minval(p(1:Nits,iy))
ymax=mean(iy)+4.*sd(iy)!maxval(p(1:Nits,iy))
yrange=(ymax-ymin)/10
yminplot=ymin-yrange
ymaxplot=ymax+yrange

if (isteplog(iy) .eq. 3) then
yminplot = -0.1
ymaxplot = 1.1
endif

!write(*,*) 'here',ix,iy,svpxleft,svpxright,svpybot,svpytop
write(*,*) 'here too',ix,iy,xminplot,xmaxplot,yminplot,ymaxplot
stop
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgswin(xminplot,xmaxplot,ymaxplot,yminplot)
!!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )

!!call pgenv(xminplot,xmaxplot,yminplot,ymaxplot,0,1)

do it=1,iburnin-1
!call pgsci(icoliburnin)
!call pgpt(1,p(it,ix),p(it,iy),-1)
enddo

do it=iburnin,Nits
!call pgsci(ipcol(it))
!call pgpt(1,p(it,ix),p(it,iy),-1)
enddo

!! plot the actual parameter values if any (if not 0)
if ((preal(ix) .ne. 0.0) .and. (preal(iy) .ne. 0.0)) then
!call pgsci(icreal)
!!call pg
!call pgpt(1,preal(ix),preal(iy),8)

endif

!call pgsci(1)



!! add the lines for mean and sd  (vertical lines for all scatter plots)
xw0(1:2)=mean(ix)
yw0(1)=yminplot
yw0(2)=ymaxplot
!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
xw0(1:2)=mean(ix)+sd(ix)
!call pgline(2,xw0,yw0)
xw0(1:2)=mean(ix)-sd(ix)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)


if (iy .eq. n) then
!! add the lines for mean and sd  (vertical lines for all scatter plots)
yw0(1:2)=mean(iy)
xw0(1)=xminplot
xw0(2)=xmaxplot
!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
yw0(1:2)=mean(iy)+sd(iy)
!call pgline(2,xw0,yw0)
yw0(1:2)=mean(iy)-sd(iy)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)
endif







!do it=1,n
!write(*,*) trim(plab(it))
!enddo
!!! decide whether to plot the axis labels

if ((ix.eq.1) .AND. (iy .ne. n)) then !! plot the axis labels on the left hand side
!call pglab('',trim(plab(iy)),'')
!call pgbox( 'bcst', 0., 10, 'bcstn', 0., 10 )

else if ((ix .eq. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),trim(plab(iy)),'')

!! decide whether to plot log histograms (inputs must already be logged)
if (isteplog(ix) .eq. 0) then  !linear
!call pgbox( 'bcstn', 0., 10, 'bcstn', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstnl', 0., 10, 'bcstn', 0., 10 )
endif

else if ((ix .ne. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),'','')
!call pgbox( 'bcstn', 0., 10, 'bcst', 0., 10 )
endif



read(*,*)



enddo ! end ix loop

enddo ! end iy loop




!!!!!! add on relevant information in the top right hand corner

svpxleft=svpxright-svpdplotx
svptop=svpup
svpbot=svpup-svpdploty

!call pgsvp(svpxleft,svpxright,svpbot,svptop)
!call pgswin(0.0,1.0,0.0,1.0)


information='N='//trim(tinits)
information=adjustl(information)

!if (n .ne. 1) call pgtext(0.5,0.5,trim(information))









read(*,*)
!call pgend

end subroutine






































































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






! scatter parms 2 (incorporates badness of fit information to produce confidence regions)

!! triangular shaped set of scatter plots showing the model in parameter space after Nits itertions
! input
!    n parms
!    isteplog(n) integer array stating whether parameter should be plotted as log (i.e for logarithmic steps) 0 for linear, 1 for log10
!    Nits number of iterations
!    p(Nits,n) stored parameters
!    bof(Nits)
!    preal(n) real array of actual input parameter values (for fake data) if set to p(:) = 0.0 do not plot



subroutine scatterparmsmulti(n,Nits,p,preal,axlab,bof,iburnin,isteplog,Ngroups)

integer n,Nits,pgopen,ier,isteplog(n),icgroup(Ngroups)
real p(Nits,n),mean(n),sd(n),sdev,avg,bof(nits),xw0(2),yw0(2),preal(n),meangroup(Ngroups,n),sdgroup(ngroups,n)
character*1000 plab(n),information,header
character*100 axlab(n)
integer ipcol(Nits),iburnin
logical info /.false./, posterplot /.true./
character*20 tinits,tsd(n),tmean(n)
character*10 tempincsd,tempincmean,tempinclo,tempinchi,tdegax,tglo,tghi
character*1000 tdeginc,tgroupinf

plotmin = 0.0
plotmax = 0.0
!write(*,*) isteplog, 'scatterparms'
cosincforcelo = 0.0  !0.866 = 30 deg
cosincforcehi = 1.0

pi=4.*atan2(1.,1.)
ddeg=10.0
deglo = 30.0
Ndeg=7
deg2rad=pi/180
igroupsize=Nits/Ngroups

if (Ngroups .ge. Nits) then
write(*,*) 'scatterparms.f90... cant have more groups than iterations (Nits)'
stop
endif

ic=0
do i=1,Ngroups !! skip out green histograms as some people are red green colour blind
ic=ic+1
icgroup(i)=ic
if (ic .eq. 2) ic = ic+1
enddo

!!! character information tnits to add as header to plots
write(tinits,'(I7.3)') nits
tinits=adjustl(tinits)

!!! color of sd and mean lines
icmean=4
icsd=4
icreal=2
icoliburnin=2
!!! annotating character size
sizeann=1.65
sizedef=1.65

if (info) then
write(*,*) 'Scatterparms input parms...',Nits
do i=1,Nits
write(*,*)i, bof(i), p(i,:)
enddo
endif



!! calculate mean of points in group
do i=1,Ngroups

idxmin=(i-1)*igroupsize +1
idxmax=idxmin+igroupsize-1
do ip=1,n
sum=0.d0
do idx=idxmin,idxmax
sum=sum+p(idx,ip)
enddo
meangroup(i,ip)=sum/(idxmax-idxmin+1)
!write(*,*) 'here multiscatterparms', i, ip, idxmin, idxmax, igroupsize, 'fff',p
sdgroup(i,ip)=sdev(p(idxmin:idxmax,ip),idxmax-idxmin+1)
enddo
enddo !! end i


!!! calculate the mean of the points
do i=1,n
sum=0.d0
sd(i)=sdev(p(iburnin:nits,i),nits-iburnin+1)
write(tsd(i),'(F10.1)') sd(i)
tsd(i)=adjustl(tsd(i))

do it=iburnin,nits
sum=sum+p(it,i)
!write(*,*) 'calculating mean',i,it,sum

enddo
mean(i)=sum/(nits-iburnin+1)
write(tmean(i),'(F10.1)') mean(i)
tmean(i)=adjustl(tmean(i))
!write(*,*) trim(tmean(i)),mean(i),sd(i),trim(tsd(i))
!read(*,*)
enddo



!write(*,*) 'hereeerrere in scatter parms...'
!!! calculate the standard deviation of points
bofmin=minval(bof)
dcs1=2.3
dcs2=6.17
dcs3=11.8
do it=1,Nits
ipcol(it)=1
if (abs(bof(it)-bofmin) .lt. dcs3) ipcol(it)=2
if (abs(bof(it)-bofmin) .lt. dcs2) ipcol(it)=3
if (abs(bof(it)-bofmin) .lt. dcs1) ipcol(it)=4
enddo
!do it=1,NIts
!write(*,*) bofmin,bof(it)
!enddo
!stop



!! trim the labels
do i=1,n
plab(i)=axlab(i)
plab(i)=adjustl(plab(i))

enddo


!ier=pgopen('parmparmplot.PS/CPS')!

!call pgslw(6)
!call pgsch(sizedef)
if (posterplot) then
!call pgsch(1.4)
!else
!!call pgsch(sizedef)
endif

!call pgsfs(2)

!call pgbgw
svpup=0.9
svpdown=0.15
svpleft=0.1
svpright=0.90
!!call pgsvp(svpleft,svpright,svpdown,svpup)


svpdploty=(svpup-svpdown)/n
svpdplotx=(svpright-svpleft)/n

!write(*,*) 'now here'

do ix=1,n

svpxright=svpleft+ix*svpdplotx
svpxleft=svpxright-svpdplotx




do iy=1,n
if (iy .lt. ix) cycle ! Do not mirror the plot









!if (isteplog(ix) .ne. 3) then
xmin=mean(ix)-4*sd(ix)!minval(p(1:nits,ix))
xmax=mean(ix)+4*sd(ix)!maxval(p(1:Nits,ix))
xrange=(xmax-xmin)/10
!else
!xmin=0.0
!xmax=1.0
!xrange=0.0
!endif

xminplot=xmin-xrange
xmaxplot=xmax+xrange

if (isteplog(ix) == 3) then
xminplot = -0.1
xmaxplot = 1.1
endif


svpybot=svpup-svpdploty*iy
svpytop=svpybot+svpdploty

!!!! Include histograms
if (iy .eq. ix) then
Nbin=100
histmax=real(Nits)/Nbin*5.
histmin=0.0



!write(*,*) 'hereeerrere in scatter parms... 222222'

!!!! Add the mean and standard deviations
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgbox( '', 0., 0, '', 0., 0 )
!call pgswin(xminplot,xmaxplot,histmin,histmax)

!write(*,*) 'argh!!!!',ix,iy,xminplot,xmaxplot, svpxleft,svpxright,svpybot,svpytop, histmin, histmax
!read(*,*)

!!!! Add the histogram


!call pgsvp(svpxleft,svpxright,svpybot,svpytop)

! The next if statements make it so the last histogram on the x axis is plotted rotated by 90 degrees)
if (ix .ne. n .or. n .eq. 1) then   !! if only one parameter varied or not on last histogram along
!call pgswin(xminplot,xmaxplot,histmin,histmax)



igroupsize=Nits/Ngroups
do igroup=1,Ngroups  !! have a do loop here to deal with different histograms for different groups
idxmin=(igroup-1)*igroupsize +1
idxmax=idxmin+igroupsize-1
!call pgsci(icgroup(igroup))
!call PGHIST(idxmax-idxmin+1, p(idxmin:idxmax,ix), xminplot,xmaxplot , Nbin, 5)
!call pgsci(1)

!write(*,*) xminplot,xmaxplot,histmin,histmax,'histogram', n, ix, isteplog(ix)
!write(*,*) isteplog
!read(*,*)
enddo



!write(*,*) 'hereeerrere in scatter parms... 3'



if (n .ne. 1) then !! if more than one parameter varied and not on last one along


!write(*,*) 'hereeerrere in scatter parms... 3b'
if (isteplog(ix) .ne. 1) then  !linear
!call pgbox( 'bcst', 0., 10, 'bc', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstl', 0., 10, 'bc', 0., 10 )
endif

!write(*,*) 'hereeerrere in scatter parms... 4'

else  !!if only one parameter varied


!write(*,*) 'hereeerrere in scatter parms... 4b'
!! decide whether to plot log histograms (inputs must already be logged)
if (isteplog(ix) .ne. 1) then  !linear
!call pgbox( 'bcstn', 0., 10, 'bc', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstnl', 0., 10, 'bc', 0., 10 )
endif

!write(*,*) 'hereeerrere in scatter parms... 4c',axlab

!call pglab(trim(axlab(1)),'Number','')
endif


!write(*,*) 'hereeerrere in scatter parms... 5'
do igroup=1,Ngroups
!call pgsci(icgroup(igroup))!! add the lines for mean and sd
xw0(1:2)=meangroup(igroup,ix)
yw0(1)=histmin
yw0(2)=histmax
!!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!!call pgsci(icsd)
!!call pgsls(3)
!!call pgsch(1.3)
!xw0(1:2)=meangroup(igroup,ix)+sdgroup(igroup,ix)
!!call pgline(2,xw0,yw0)
!xw0(1:2)=meangroup(igroup,ix)-sdgroup(igroup,ix)
!!call pgline(2,xw0,yw0)
!call pgsci(1)
!!call pgsls(1)
!!call pgsch(1.0)
enddo















!write(*,*) 'i realy hate science',ix,n

else if (ix .eq. n) then


do igroup=1,Ngroups  !! have a do loop here to deal with different histograms for different groups
idxmin=(igroup-1)*igroupsize +1
idxmax=idxmin+igroupsize-1
nidx = idxmax - idxmin +1
!call pgsci(icgroup(igroup))
if (isteplog(ix) == 3) then  !! force axis on histogram plot to be 0 to 1 if stepping in cos theta
iplotlimforce = 1
if (cosincforcelo .gt. 0) then
xminplot = cosincforcelo
xmaxplot = cosincforcehi
else
xminplot = -0.1
xmaxplot = 1.1
endif
!write(*,*) 'pghistalong plotlims before routine', plotmin, plotmax,xminplot,xmaxplot
!write(*,*) ix,'here'
!write(*,*) 'argh!!!!',xminplot,xmaxplot,plotmin,plotmax! svpxleft,svpxright,svpybot,svpytop
!read(*,*)

else
iplotlimforce = 0
endif

!write(*,*) 'before histlong'
!call pghistalong(Nidx, p(idxmin:idxmax,ix), xminplot, xmaxplot,plotmin,plotmax,dbin, Nbin, 2, 2.0,iplotlimforce) ! pgswin called inside this routine
!write(*,*) 'after hist along'


if (isteplog(ix) == 3) then
!call pgslw(8)
!call pgsci(1)
xw0(1) = plotmin
xw0(2) = plotmax
yw0(1:2) = 0.0
!call pgline(2, xw0,yw0)
yw0(1:2) = 1.0
!call pgline(2,xw0, yw0)
!call pgslw(6)
endif

!write(*,*) 'hereeerrere in scatter parms... 6'

!call pgsci(1)
enddo
!write(*,*) 'horhist',xminplot,xmaxplot
!!call pglab(trim(plab(ix)),'','')

if (isteplog(ix) .ne. 1) then  !linear
!call pgbox( 'bc', 0., 10, 'bcstn', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bc', 0., 10, 'bcstl', 0., 10 ) !! on histograms tipped over, x axis is yaxis so do not put ticks on the xaxis which indicates (number) or posterior probability distribution
endif


!! If we are plotting cos(inc) isteplog(i)=3. Then we add on the plot annotations to show the mean and lower and upper confidence limits in degrees
if (isteplog(ix) .eq. 3) then


!write(tempinclo,'(F10,2)') 180/pi*acos(mean(ix)-sd(ix))
write(tempincmean,'(I10)') int(180/pi*acos(mean(ix)))
tempincmean=adjustl(tempincmean)
!write(tempinchi,'(F10,2)') 180/pi*acos(mean(ix)+sd(ix))
if (sd(ix) .gt. 0) then
write(tempincsd,'(I10)') int(abs(sd(ix)/sin(acos(mean(ix))))*180/pi)!! by propagation of uncertainty      !abs( 180/pi*( acos(mean(ix)+sd(ix)) - acos(mean(ix)-sd(ix)) )/2 )
else
tempincsd='NOO!'
endif

tempincsd=adjustl(tempincsd)
!write(*,*)   acos(mean(ix)+sd(ix)) ,acos(mean(ix)-sd(ix)),'incs'
!read(*,*)
tdeginc=trim('i\uo\d = '//trim(tempincmean)//' \(2233) '//trim(tempincsd))
tdeginc=adjustl(tdeginc)
!call pgsch(sizeann)
!call pgmtxt('T',2.0,0.5,0.5,trim(tdeginc))
!call pgsch(sizedef)
endif
!!



!! add on an extra axis if cosine parameter (so that you can see the cosine in degrees)
do i=1,Ndeg
deg=deglo + ddeg*(i-1)
yw0(1:2)=cos(deg*deg2rad)

if ((yw0(1) .gt. xminplot) .and. (yw0(1) .lt. xmaxplot)) then
write(tdegax,'(I0)') int(deg)
tdegax=trim(trim(tdegax)//'\uo')
tdegax=adjustl(tdegax)
xw0(1)=plotmin+0.9*(plotmax-plotmin)
xw0(2)=plotmax
!call pgline(2,xw0,yw0)
!call pgmtxt('RV',2.0,(yw0(1)-xmaxplot)/(xminplot-xmaxplot),1.0,trim(tdegax))
endif

enddo
!!


!write(*,*) 'hereeerrere in scatter parms... 7'


!!! add the lines for mean and sd
!yw0(1:2)=mean(ix)
!xw0(1)=plotmin
!xw0(2)=plotmax
!!call pgsci(icmean)
!!call pgline(2,xw0,yw0)
!!call pgsci(icsd)
!!call pgsls(3)
!!!call pgsch(1.3)
!yw0(1:2)=mean(ix)+sd(ix)
!!call pgline(2,xw0,yw0)
!yw0(1:2)=mean(ix)-sd(ix)
!!call pgline(2,xw0,yw0)
!!call pgsci(1)
!!call pgsls(1)
!!!call pgsch(1.0)


do igroup=1,Ngroups
!call pgsci(icgroup(igroup))!! add the lines for mean and sd
yw0(1:2)=meangroup(igroup,ix)
xw0(1)=plotmin
xw0(2)=plotmax
!!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
yw0(1:2)=meangroup(igroup,ix)+sdgroup(igroup,ix)
!!call pgline(2,xw0,yw0)
yw0(1:2)=meangroup(igroup,ix)-sdgroup(igroup,ix)
!!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)
enddo







endif




header=trim(tmean(ix))//'\(2233)'//trim(tsd(ix))

if (isteplog(ix) .eq. 1) header='10\u'//trim(header) !! make the header into 10^ form

if (n .eq. 1) header=trim(header)//'    N:'//trim(tinits)
header=adjustl(header)
!call pgmtxt('T',0.2,0.5,0.5,trim(header))


!write(*,*) ix,iy,svpxleft,svpxright,svpybot,svpytop,xminplot,xmaxplot,nits

cycle
endif
!!!


!!!! make the scatter plot
if (isteplog(iy) .ne. 3) then
ymin=mean(iy)-4.*sd(iy)!minval(p(1:Nits,iy))
ymax=mean(iy)+4.*sd(iy)!maxval(p(1:Nits,iy))
yrange=(ymax-ymin)/10
else

ymin =0.0! mean(iy)
ymax =1.0! mean(iy)

endif

!else
!ymin=0.0
!ymax=1.0
!yrange=0.0
!endif



yminplot=ymin-yrange
ymaxplot=ymax+yrange

if (isteplog(iy) .eq. 3) then
if (cosincforcelo .gt. 0) then
yminplot = cosincforcelo
ymaxplot = cosincforcehi
else
yminplot = -0.1
ymaxplot = 1.1
endif
endif

!write(*,*) 'here',ix,iy,svpxleft,svpxright,svpybot,svpytop
!write(*,*) 'here too',ix,iy,xminplot,xmaxplot,yminplot,ymaxplot
!call pgsvp(svpxleft,svpxright,svpybot,svpytop)
!call pgswin(xminplot,xmaxplot,ymaxplot,yminplot)
!call pgbox( 'bcst', 0., 10, 'bcst', 0., 10 )

!!call pgenv(xminplot,xmaxplot,yminplot,ymaxplot,0,1)


do igroup=1,Ngroups  !! have a do loop here to deal with different histograms for different groups
idxmin=(igroup-1)*igroupsize +1
idxmax=idxmin+igroupsize-1
!call pgsci(icgroup(igroup))
do it=idxmin,idxmax
!call pgpt(1,p(it,ix),p(it,iy),-1)
enddo
enddo
!write(*,*) xminplot,xmaxplot,yminplot,ymaxplot,'scatter'
!write(*,*) svpxleft, svpxright



!! plot the actual parameter values if any (if not 0)
if ((preal(ix) .ne. 0.0) .and. (preal(iy) .ne. 0.0)) then
!call pgsci(icreal)
!!call pg
!call pgpt(1,preal(ix),preal(iy),8)

endif

!call pgsci(1)






do igroup=1,Ngroups
!call pgsci(icgroup(igroup))!! add the lines for mean and sd
xw0(1:2)=meangroup(igroup,ix)
yw0(1)=yminplot
yw0(2)=ymaxplot
!!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
xw0(1:2)=meangroup(igroup,ix)+sdgroup(igroup,ix)
!call pgline(2,xw0,yw0)
xw0(1:2)=meangroup(igroup,ix)-sdgroup(igroup,ix)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)
enddo





if (iy.eq. n) then
do igroup=1,Ngroups
!call pgsci(icgroup(igroup))!! add the lines for mean and sd
yw0(1:2)=meangroup(igroup,iy)
xw0(1)=xmaxplot
xw0(2)=xminplot
!!call pgsci(icmean)
!call pgline(2,xw0,yw0)
!!call pgsci(icsd)
!call pgsls(3)
!!call pgsch(1.3)
yw0(1:2)=meangroup(igroup,iy)+sdgroup(igroup,iy)
!call pgline(2,xw0,yw0)
yw0(1:2)=meangroup(igroup,iy)-sdgroup(igroup,iy)
!call pgline(2,xw0,yw0)
!call pgsci(1)
!call pgsls(1)
!!call pgsch(1.0)
enddo
end if






!do it=1,n
!write(*,*) trim(plab(it))
!enddo
!!! decide whether to plot the axis labels

if ((ix.eq.1) .AND. (iy .ne. n)) then !! plot the axis labels on the left hand side
!call pglab('',trim(plab(iy)),'')
!call pgbox( 'bcst', 0., 10, 'bcstn', 0., 10 )

else if ((ix .eq. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),trim(plab(iy)),'')

!! decide whether to plot log histograms (inputs must already be logged)
if (isteplog(ix) .eq. 0) then  !linear
!call pgbox( 'bcstn', 0., 10, 'bcstn', 0., 10 )
else if (isteplog(ix) .eq. 1) then !log10
!call pgbox( 'bcstnl', 0., 10, 'bcstn', 0., 10 )
endif

else if ((ix .ne. 1) .AND. (iy .eq. n)) then
!call pglab(trim(plab(ix)),'','')
!call pgbox( 'bcstn', 0., 10, 'bcst', 0., 10 )
endif



!read(*,*)



enddo ! end ix loop

enddo ! end iy loop




!!!!!! add on relevant information in the top right hand corner

svpxleft=svpxright-svpdplotx
svptop=svpup
svpbot=svpup-svpdploty

!call pgsvp(svpxleft,svpxright,svpbot,svptop)
!call pgswin(0.0,1.0,0.0,1.0)




if (n .ne. 1) then

do igroup=1,Ngroups
idxmin=(igroup-1)*igroupsize +1
idxmax=idxmin+igroupsize-1



tgroupinf=' '
write(tglo,'(I10.5)') idxmin
write(tghi,'(I10.5)') idxmax
tglo=adjustl(tglo)
tghi=adjustl(tghi)


tgroupinf=trim(tgroupinf)//' '//trim(tglo)//' - '//trim(tghi)
tgroupinf=adjustl(tgroupinf)
!call pgsci(icgroup(igroup))
!call pgmtxt('T',-2.0-1.*igroup,0.5,0.5,trim(tgroupinf))




if (igroup .eq. Ngroups) then
information='N='//trim(tghi)
information=adjustl(information)
!call pgsci(1)
!call pgmtxt('T',-1.0,0.5,0.5,trim(information))
endif

enddo
endif !! end the additional annotations









!read(*,*)
!call pgend

!write(*,*) 'end of scatter routine'
end subroutine





















!!!! test scatter parms

!program testscatter
!
!integer n
!real,allocatable::p(:,:)
!character*1000,allocatable::axlab(:)
!character*3 tr
!rms=4.0
!n=4
!Nits=5000
!allocate(p(Nits,n),axlab(n))
!
!
!iseed=2314344
!
!pmin=10.0
!pmax=80.0
!
!do i=1,n
!
!write(tr,'(I1.1)') i
!tr=adjustl(tr)
!axlab(i)='parm '//tr
!
!pmean=ranu(pmin,pmax,iseed)
!do it=1,Nits
!
!p(it,i)=rang(pmean,rms,iseed)
!enddo
!
!
!enddo
!
!do it=1,Nits
!write(*,*) p(it,:)
!enddo
!call scatterparms(n,Nits,p,axlab)
!
!end program
!
!!! gfortran scatterparms.f90 -fbounds-check -fno-backtrace -Wuninitialized -L/Users/ds207/software/pgplot -lpgplot -L/usr/X11/lib -lX11 -lpng ran_new2.for pgbgw.for avg.f90
!
!
!
!
!
!






!!gfortran pghistalong.f90 -fbounds-check -fno-backtrace -Wuninitialized -L/Users/ds207/software/pgplot -lpgplot -L/usr/X11/lib -lX11 -lpng ran_new2.for
!!! subroutine to plot a horizontal histogram
!! view port must already be set

!!! INPUTS x(nx), xmin, xmax the array of data, its minimum and maximum respectively
!			nbin the desired number of bins
!			pgflag (if 2, will plot the histogram with x and y axis reversed i.e bins start from left of plot and go out)
!			pgoutline if 2, draws only outline of histogram new 20/10/2014
!			plotlimforce: if 1, plotmin and plotmax are inputs rather than outputs
subroutine pghistalong(nx, x, xmin, xmax, plotmin, plotmax, dx, Nbin, pgflag,pgoutline,iplotlimforce)

integer pgflag, nx, nbin, iplotlimforce
real x(nx), bin(nbin), x2(2),y2(2)


bin(:)=0.0
!! bin spacing
dx = (xmax-xmin)/Nbin


!! calculate the number of elements in each bin
do i=1,nx
!if (x(i) .lt. 0) then
!write(*,*) x(i), i,'less than 0!!'
!stop
!endif

do i2 =1,nbin

binlow=xmin+dx*(i2-1)
binhi=binlow+dx

!write(*,*) i2,
if ( (x(i) .ge. binlow) .and. (x(i) .lt. binhi) ) then
bin(i2) = bin(i2) + 1.0
endif


!if (i .eq. nx) write(*,*) i2,binlow,binhi,bin(i2)

enddo
enddo



!! end bin assignments
ymax=maxval(bin)
ymin=0.0
yextra=ymax/2
!write(*,*) xmin, xmax, ymin, ymax,dx


!if (iplotlimforce == 1) then

!write(*,*) 'pghistalong plotlims', plotmin, plotmax,xmin,xmax
!read(*,*)
!endif
if (pgflag .eq. 2) then

plotmin=ymin
plotmax=ymax+yextra



!call pgswin(plotmin, plotmax, xmax, xmin)
endif
!! plot the rectangles
do in=1,nbin

if (pgflag .eq. 2) then
ybot=xmin+dx*(in-1)
ytop=ybot+dx
xleft=0.0
xright=bin(in)
endif


!!!!!!!!!!!!!!!!!!!!!!!!!! The next few lines plot only the outline
if (pgoutline ==2) then



if (in .eq. 1) then !! if first bin draw left line
x2(1) = xleft
x2(2) = xright
y2(1:2) = ybot
!call pgline(2,x2,y2)
endif

x2(1:2)=xright !! always draw top vertical bar
y2(1) = ybot
y2(2) = ytop
!call pgline(2,x2,y2)

if (in .lt. nbin) then
xrightnext = bin(in+1)
!if rightcondition) then  !!draw right bar if last bin
x2(1) = xright
x2(2) = xrightnext
else
x2(1) = xleft
x2(2) = xright
endif
y2(1:2) = ytop
!call pgline(2,x2,y2)


else


!call pgrect(xleft,xright,ybot,ytop)

!!!!!!!!!!!!!!!!!!!!!!
endif !end if pgoutline

enddo


end subroutine








!!!!! test pghistalong
!program testpghistalong
!integer ier,pgopen
!real dat(1000),rang
!
!
!n=1000
!sd=10.0
!av=320.0
!iseed=3124523
!
!do i=1,n
!dat(i)=rang(av,sd,iseed)
!enddo
!
!xright=0.95
!xleft=0.05
!ybot=0.05
!ytop=0.95
!Nbin=10
!xmin=minval(dat)
!xmax=maxval(dat)
!dx=(xmax-xmin)/Nbin
!do i=1,n
!write(*,*)i, dat(i)
!enddo
!
!
!!ier=pgopen('/XServe')
!
!!call pgsvp(xleft,xright,ybot,ytop)
!!!call pgswin(xleft,xright,ybot,ytop)
!
!!call pghistalong(n,dat,xmin,xmax,Nbin, 2)
!!call pgbox( 'bcstn', 0., 0, 'bcstm', 0, 0 )
!
!read(*,*)
!!call pgend
!
!end program
!
!
!
!
!
!



!Update 22jun 15 p0 output altered so that it yields rms = 1

!! double precision for frequency and time better for higher values


!  input tdat(Ndat), dat(Ndat), sigdat(Ndat), w(NW)
!  deginc, ememdot, wav
!  input / output sk(NW), ck(NW), sigsk(NW), sigck(NW) estimates for the parameters and uncertainties
!  .......    p0, w0 not yet decided if these should be inputs or just clculated in routine
! presently w0 set to w(2) (w(1) is 0), p0 set to 12.5*amp**2xs

! consider inputting logged data to prevent model from going below 0 (unphysical for fluxes)


! NOTE does not incorporate convolution with prior at this point, just -2 power law slope
! ememdot, deginc, wav useless

!NEW 21st april 15: if w0 and p0 -ve they will be changed by the routine

subroutine iosecho(Ndat, tdat, dat, sigdat,&
NW, w, w0, p0, sk, ck, sigsk, sigck)

integer Ndat, Nw,pgbeg

real ememdot, deginc, wav, tdat(Ndat), sigdat(Ndat), dat(ndat), sw(Ndat,NW),&
cw(Ndat,NW), sk(NW), sigsk(NW), ck(NW), sigck(NW), medspac, w(NW), datcp(Ndat), pk2(NW)

real,allocatable:: xmod(:), tmod(:), xmod_inc(:)

logical diag /.True./

!sigdat(:) = 1.0e-2
pi 			= 3.14159265
twopi 		= 2*pi
Nits        = 1000
!! default parameters
dw = w(2) - w(1)
convergence = 1.e-5

!estimate w0 and p0
if (w0 .le. 0) w0 = w(2)


datmax = maxval(dat)
datmin = minval(dat)
sd = (datmax - datmin)/5


if (p0 .lt. 0) p0 = 12.5*sd*sd !see derivation in small black book




!!!! perform a large loop over all the frequencies using iterated optimal scaling as we go
datcp(1:Ndat) = dat(1:Ndat)

do iw = 1,Nw !pre calculate sin and cosine terms
wnow = w(iw)
do it = 1,Ndat
cw(it,iw) = cos(wnow*tdat(it))
sw(it,iw) = sin(wnow*tdat(it))
enddo
enddo

a = 0.5*p0*dw*w0*w0
do iw = 1,NW ! pre calculate the prior terms
wnow = w(iw)
pk2(iw) = a/(wnow*wnow)
enddo


do iw = 1,NW
cnow = ck(iw)
snow = sk(iw)
wnow = w(iw)
pk2now = pk2(iw)



do iteration = 1,nits

cnowold = cnow
snowold = snow
!if (iteration .eq. 10) write(*,*) w(iw), iteration, cnow


!do S_k first
top = 0.d0
bot = 0.d0

!write(*,*) snow, iw, iteration
do it = 1,Ndat
fnow    = datcp(it)
signow  = sigdat(it)
sn2     = signow*signow
tnow    = tdat(it)
cwnow   = cw(it,iw)
swnow   = sw(it,iw)


top = top + (fnow - cnow*cwnow)*swnow/sn2
bot = bot + swnow*swnow/sn2
enddo
!apply prior
if (wnow .ne. 0) then
a = sqrt(pk2(iw))
bot = bot + 1./(a*a)
endif


snow = top/bot
sigsnow = sqrt(1./bot)





if (wnow .eq. 0) then
snowold = 1.0
snow    = 1.0
sigsnow = 1.0
pk2(iw) = 1.0
endif


!write(*,*) snow,a, 'snow after'
!read(*,*)






!Now C_k
top = 0.d0
bot = 0.d0

!write(*,*) cnow, iw, iteration

do it = 1,Ndat
fnow    = datcp(it)
signow  = sigdat(it)
sn2     = signow*signow
!tnow    = tdat(it)
cwnow   = cw(it,iw)
swnow   = sw(it,iw)


top = top + (fnow - snow*swnow)*cwnow/sn2
bot = bot + cwnow*cwnow/sn2
enddo
!apply prior
if (wnow .ne. 0) then
a = sqrt(pk2(iw))
bot = bot + 1/(a*a)
endif

cnow = top/bot
sigcnow = sqrt(1./bot)



!snow = sqrt(pk2(iw))

ck(iw)    = cnow
sigck(iw) = sigcnow
sk(iw)    = snow
sigsk(iw) = sigsnow





! convergence test here
a = abs((snow - snowold)/sigsnow)
b = abs((cnow - cnowold)/sigcnow)
if ((a .lt. convergence) .and. (b .lt. convergence) .and. (iteration .gt. 10)) then
!write(*,*) 'Convergence after',iteration, 'iterations', snow, snowold, a
exit
endif

!write(*,*) iteration, iw, snow,cnow,a, b!, snow, snowold, sigsnow

enddo !end iteration

!read(*,*)

datcp(1:Ndat) = datcp(1:Ndat) - snow*sw(1:Ndat,iw) - cnow*cw(1:Ndat,iw)

if (wnow .eq. 0)then !!0 freq check
sum = 0.d0
do it = 1,Ndat
sum = sum + dat(it)
enddo
dmean = sum/Ndat
write(*,*) dmean, ck(iw), sk(iw), 'mean, 0freq ck and sk'
endif


enddo ! end iw







if (diag) then

tlo = tdat(1)
thi = tdat(ndat)
xlo = minval(dat)
xhi = maxval(dat)
sigxlo = minval(sigdat)
sigxhi = maxval(sigdat)

!!call pgenv(tlo,thi,xlo,xhi,0,0)
!!call pgpt(Ndat, tdat, dat, 6)
NTmod = Ndat*10
dtmod = (thi-tlo)/(NTmod-1)
allocate(tmod(NTmod), xmod(NTmod), xmod_inc(Ntmod))

do it = 1,Ntmod
tmod(it) = tlo + dtmod*(it-1)
enddo
xmod_inc(:) = 0.d0


!do it = 1,NTmod
!sum = 0.d0
!do iw = 1,NW
!sum = sum + sk(iw)*sin(w(iw)*tmod(it)) + ck(iw)*cos(w(iw)*tmod(it))
!
!!! make incremental plpot
!do it1 = 1,Ntmod
!xmod_inc(it1) = xmod_inc(it1) + sk(iw)*sin(w(iw)*tmod(it1)) + ck(iw)*cos(w(iw)*tmod(it1))
!!!write(*,*) tmod(it1), xmod_inc(it1)
!enddo
!!ier = pgbeg(0,'/Xserve',1,2)
!!call pgenv(tlo,thi,xlo,xhi,0,0)
!!call pgline(Ntmod,tmod,xmod_inc)
!!call pgpt(Ndat, tdat, dat, 6)
!!call pgenv(tlo,thi,sigxlo,sigxhi,0,0)
!!call pgpt(Ndat, tdat, sigdat, 6)
!!
!
!!do it2 = 1,Ntmod
!!write(*,*) w(iw),tmod(it2), sk(iw)*sin(w(iw)*tmod(it2)) + ck(iw)*cos(w(iw)*tmod(it2))
!!enddo
!
!!read(*,*)
!!call pgend
!enddo
!xmod(it) = sum
!write(*,*) tmod(it), xmod(it)
!if (xmod(it) .ne. xmod(it)) then
! do iw  =1,NW
! write(*,*) sk(iw), ck(iw)
! enddo
! write(*,*) 'nan problem model ios_echo2.f90'
! stop
!endif


!enddo


!!call pgline(Ntmod,tmod,xmod)

!xmodlo = minval(xmod)
!xmodhi = maxval(xmod)
wplo = w(1)
wphi = w(NW)
pslo = -11.
pshi = 2.0
!!ier = pgbeg(0,'pspec_spike.ps/CPS',1,1)
!!call pgenv(alog10(w(2)/twopi),alog10(wphi/twopi),pslo,pshi,0,0)
!twopi = 2*3.14159265
!do iw = 1,NW
!! call pgpt(1, alog10(w(iw)/twopi), alog10(sk(iw)**2 + ck(iw)**2), 6)
!enddo
!!call pgmtxt('B',2.0,0.5,0.5,'log\d10\u(frequency / [cycles/day])')
!!call pgmtxt('L',2.0,0.5,0.5,'log\d10\u(P(f))')
!!call pgend
!
endif






!!!!! New 22nd June 2015!!!!!!!!!!!
!sum = 0.d0
!do iw = 1,NW
! sum = sum + 1./w(iw)**2
!enddo
!p0 =  2./dw/w0**2 / sum
!sk(1:NW) = sk(1:NW)/sqrt(p0/2)
!sigsk(1:NW) = sigsk(1:NW)/sqrt(p0/2)
!ck(1:NW) = ck(1:NW)/sqrt(p0/2)
!sigck(1:NW) = sigck(1:NW)/sqrt(p0/2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





end subroutine







!!! testing program gflags iosecho3.f90 ../dft.for ../delaydist/tfb.f90 ../medspac.for ../cgs.for ../ran_new2.for ../numline.for
!program testiosecho
!
!real, allocatable,dimension(:):: tdat, dat, sigdat, sk, ck, sigck, sigsk, w,&
!tmod,xmod,pspec
!integer pgbeg, numline
!character(100) fname,dir
!
!dir   = '../mcmcmulti3/fake_snr_300.0/'
!dir   = adjustl(dir)
!call chdir(trim(dir))
!
!
!fname =  '300.0_resamp_fake_3556.51.dat'
!fname =  adjustl(fname)
!Ndat  = numline(trim(fname)) - 1
!
!write(*,*) Ndat
!allocate(tdat(Ndat), dat(ndat), sigdat(Ndat))
!
!open(unit = 1,file=trim(fname))
!do i = 1,Ndat
!read(1,*) tdat(i), dat(i), sigdat(i)
!!write(*,*)tdat(i), dat(i), sigdat(i)
!enddo
!close(1)
!!stop
!tlo = tdat(1)
!thi = tdat(Ndat)
!dt  = (thi - tlo)/(Ndat - 1)
!tlen = thi - tlo
!
!
!pi = 3.14159265
!twopi = pi*2
!
!wav = 4000.0
!ememdot = 1.e8
!deginc = 30.0
!
!
!! load in a fake dataset
!
!
!!t
!!tlo    = 0.0
!!thi    = 100.0
!!tlen   = thi - tlo
!!dt     = 1.0
!!Nt     = (thi - tlo)/dt + 1
!
!wlo    = twopi/tlen * 0.5
!whi    = twopi/dt * 10.
!dw     = wlo
!NW     = ceiling((whi - wlo)/dw) + 1
!
!!write(*,*) 'NW;,',dw, whi!, wlo
!allocate(w(NW), sk(NW), sigsk(NW), sigck(NW), ck(NW))
!do iw = 1,NW
!w(iw) = (iw-1)*dw
!enddo
!
!
!
!sk(1:NW) = 1
!ck(1:NW) = 1
!
!!call iosecho(Ndat, tdat, alog(dat(1:Ndat)), sigdat/dat(1:Ndat), wav, ememdot, deginc,&
!! NW, w, w0, p0, sk, ck, sigsk, sigck)
!
!call iosecho(Ndat, tdat, dat(1:Ndat), sigdat, wav, ememdot, deginc,&
! NW, w, w0, p0, sk, ck, sigsk, sigck)
!
!
!
!
!
!
!
!
!
!!! plot the results
!Ntmod = Ndat*100
!allocate(tmod(NTmod), xmod(NTmod))
!dtmod = (thi - tlo)/(1.*Ntmod - 1.)
!do it = 1,NTmod
!tmod(it) = (it-1)*dtmod + tlo
!sum = 0.d0
!do iw = 1,NW
!sum = sum + sk(iw)*sin(w(iw) * tmod(it)) + ck(iw)*cos(w(iw)*tmod(it))
!enddo
!!xmod(it) = exp(sum)
!xmod(it) = sum
!!write(*,*) it,tmod(it),xmod(it)
!
!enddo
!
!!do it =1,Ndat
!!write(*,*) tdat(it), dat(it)
!!enddo
!
!xlo = min(minval(xmod), minval(dat))
!xhi = max(maxval(xmod), maxval(dat))
!
!
!!ier = pgbeg(0,'/Xserve',1,2)
!!call pgenv(tlo,thi,xlo,xhi,0,0)
!!call pgpt(Ndat, tdat, dat, 6)
!!call pgline(Ntmod,tmod,xmod)
!
!
!
!allocate(pspec(NW))
!do iw = 1,NW
!pspec(iw) = sk(iw)**2 + ck(iw)**2
!enddo
!
!xlo = alog10(minval(pspec))
!xhi = alog10(maxval(pspec))
!wlo = alog10(w(2)/twopi)
!whi = alog10(w(NW)/twopi)
!
!!call pgenv(wlo,whi,xlo,xhi,0,-2)
!!call pgbox( 'bcnstl', 0., 0.0, 'bcnstl', 0., 0. )
!!call pgpt(NW, alog10(w/twopi), alog10(pspec), 6)
!!call pgend
!
!end program






!! return the minimum spacing between elements in an array


subroutine minspac(N,x,dxmin)

integer N
real x(N),dxmin


do i = 2,N,2

xhi = x(i)
xlo = x(i-1)
dx  = xhi - xlo


if (i.eq. 2) then
dxmin = dx
else
dxmintemp = dx
if (dxmintemp .lt. dxmin) dxmin = dxmintemp
endif

enddo

end subroutine



subroutine corpar(NW,cov,cor)
integer NW
real cov(NW,NW),cor(NW,NW)


do iw2 = 1,NW
var2 = cov(iw2,iw2)
sd2  = sqrt(var2)

do iw1 = 1,NW
var1 = cov(iw1,iw1)
sd1  = sqrt(var1)

cor(iw1,iw2) = cov(iw1,iw2)/(sd1*sd2)

!if(cor
!if (iw1 .eq. iw2 .and. iw1 .eq. 801) then
!write(*,*) sd1, sd2, var1, var2, cov(iw1,iw2), 'corpar', cor(iw1,iw2)
!read(*,*)
!endif

enddo
enddo

end subroutine






subroutine covpar(NW,Nrep,x,cov)

integer NW,Nrep
real x(NW,Nrep), cov(NW,NW),xmean(NW)



!! evaluate the parameter mean
do iw = 1,NW
sum = 0.d0
do irep = 1,Nrep
sum = sum + x(iw,irep)
enddo
xmean(iw) = sum /Nrep
enddo

!

! evaluate covariance

!write(*,*) 'dkjldf',Nrep,NW
do iw2 = 1,NW
xm2 = xmean(iw2)
do iw1 = 1,NW

xm1 = xmean(iw1)

covsum = 0.d0
do irep = 1,Nrep
!write(*,*) 'cov:',iw1,iw2,NW,irep,Nrep
covsum = covsum + (x(iw1,irep) - xm1) * (x(iw2,irep) - xm2)

enddo
cov(iw1,iw2) = covsum / Nrep

enddo
enddo


end subroutine














!like covpar but reverses row and columns of input cov matrix
!... faster due to fortran array indexing
subroutine covpar_2018(NW,Nrep,x,cov)

integer NW,Nrep
real x(Nrep,NW), cov(NW,NW),xmean(NW)



!! evaluate the parameter mean
do iw = 1,NW
sum = 0.d0
do irep = 1,Nrep
sum = sum + x(irep,iw)
enddo
xmean(iw) = sum /Nrep
enddo

!

! evaluate covariance

!write(*,*) 'dkjldf',Nrep,NW
do iw2 = 1,NW
xm2 = xmean(iw2)
do iw1 = 1,NW

xm1 = xmean(iw1)

covsum = 0.d0
do irep = 1,Nrep
!write(*,*) 'cov:',iw1,iw2,NW,irep,Nrep
covsum = covsum + (x(irep,iw1) - xm1) * (x(irep,iw2) - xm2)

enddo
cov(iw1,iw2) = covsum / Nrep

enddo
enddo


end subroutine





!!! calculate index of min and max of 1d array
!! inp x(N)
!! op idxmin,idxmax

subroutine minmax(x,N,idxmin,idxmax)

integer N, idxmin, idxmax
real x(N)

do i = 1,N
if (i .eq. 1) then
xmin = x(1)
xmax = xmin
idxmin = 1
idxmax = 1
else
xmintemp = x(i)
xmaxtemp = xmintemp

if (xmintemp .lt. xmin) then
xmin = xmintemp
idxmin = i
endif

if (xmaxtemp .gt. xmax) then
xmax = xmaxtemp
idxmax = i
endif

endif

enddo

end subroutine
!!!!!!!!!!!!!!!!!!!



!! calculate the min and max value of 2d array fortran

!! INPUTS x(N1,N2)

!! OP idxlo1,idxlo2, idxhi1, idxhi2 the indicees of the min and max value
subroutine minmax2d(x,N1,N2,idxlo1,idxlo2,idxhi1,idxhi2)

integer N1, N2, idxlo1, idxlo2, idxhi1, idxhi2
real xlo, xhi, x(N1,N2)

!! minval


do i2 = 1,N2
do i1 =1,N1


if ((i1 .eq. 1) .and. (i2 .eq. 1)) then

xmin = x(1,1)
xmax = xmin
idxlo1 = 1
idxlo2 = 1
idxhi1 = 1
idxhi2 = 1

else

xmintemp = x(i1,i2)
xmaxtemp = xmintemp

if (xmintemp .lt. xmin) then
xmin = xmintemp
idxlo1 = i1
idxlo2 = i2
endif

if (xmaxtemp .gt. xmax) then
xmax = xmaxtemp
idxhi1 = i1
idxhi2 = i2
endif

endif

enddo !end i2
enddo ! end i1
!write(*,*) 'minmax2d report...', idxhi1,idxhi2,xmax, idxlo1,idxlo2,xmin,x(1,1)
end subroutine

















!! calculate the min and max value of 3d array fortran

!! INPUTS x(N1,N2,N3)

!! OP idxlo1,idxlo2, idxhi1, idxhi2 the indicees of the min and max value
subroutine minmax3d(x,N1,N2,N3,idxlo1,idxlo2,idxlo3,idxhi1,idxhi2,idxhi3)

integer N1, N2,N3, idxlo1, idxlo2, idxhi1, idxhi2
real xlo, xhi, x(N1,N2,N3)

!! minval

do i3 = 1,N3
do i2 = 1,N2
do i1 =1,N1


if ((i1 .eq. 1) .and. (i2 .eq. 1) .and. (i3.eq.1)) then

xmin = x(1,1,1)
xmax = xmin
idxlo1 = 1
idxlo2 = 1
idxhi1 = 1
idxhi2 = 1

else

xmintemp = x(i1,i2,i3)
xmaxtemp = xmintemp

if (xmintemp .lt. xmin) then
xmin = xmintemp
idxlo1 = i1
idxlo2 = i2
idxlo3 = i3
endif

if (xmaxtemp .gt. xmax) then
xmax = xmaxtemp
idxhi1 = i1
idxhi2 = i2
idxhi3 = i3
endif

endif

enddo ! end i3
enddo !end i2
enddo ! end i1
!write(*,*) 'minmax2d report...', idxhi1,idxhi2,xmax, idxlo1,idxlo2,xmin,x(1,1)
end subroutine






















!! calculate the min and max value of 4d array fortran

!! INPUTS x(N1,N2,N3,N4)

!! OP idxlo1,idxlo2, idxhi1, idxhi2 the indicees of the min and max value
subroutine minmax4d(x,N1,N2,N3,N4,idxlo1,idxlo2,idxlo3,idxlo4,&
idxhi1,idxhi2,idxhi3,idxhi4)

integer N1, N2,N3, idxlo1, idxlo2, idxlo3, idxlo4, idxhi1, idxhi2,idxhi3, idxhi4
real xlo, xhi, x(N1,N2,N3,N4)

!! minval

do i4 = 1,N4
do i3 = 1,N3
do i2 = 1,N2
do i1 = 1,N1


if ((i1 .eq. 1) .and. (i2 .eq. 1) .and. (i3.eq.1) .and. (i4.eq.1)) then

xmin = x(1,1,1,1)
xmax = xmin
idxlo1 = 1
idxlo2 = 1
idxlo3 = 1
idxlo4 = 1
idxhi1 = 1
idxhi2 = 1
idxhi3 = 1
idxhi4 = 1

else

xmintemp = x(i1,i2,i3,i4)
xmaxtemp = xmintemp

if (xmintemp .lt. xmin) then
xmin = xmintemp
idxlo1 = i1
idxlo2 = i2
idxlo3 = i3
idxlo4 = i4
endif

if (xmaxtemp .gt. xmax) then
xmax = xmaxtemp
idxhi1 = i1
idxhi2 = i2
idxhi3 = i3
idxhi4 = i4
endif

endif

enddo ! end i4
enddo ! end i3
enddo ! end i2
enddo ! end i1




end subroutine



!! code to make parameter vs iteration plots
!! makes a series of vertically stacked line plots with par(1:nits,ipar) on each row
!! fname what to name output ps file
!! xlab (character len=100) what to put on y axis
!! par(nits,np) total pparameter array
!! idxpar(npar) idx of total par array you want to make into line plots
!! ymin_in(npar), ymax_in(npar) array of minimum and maximium values on y axis you
! want to plot (set to 666 to let the program chose)





subroutine mypgline(parmain,nits,fname,idxpar,npar,nparmain,axlab,ymin_in,ymax_in)

character(len=100) axlab(npar),fname
integer npar,idxpar(npar),nparmain
real parmain(nparmain,nits),x_its(nits),ymin_in(npar),ymax_in(npar),&
par(nits,npar)



!write(*,*) 'somewhere u phere'
do ip = 1,npar
idp = idxpar(ip)
do it = 1,nits
par(it,ip) =  parmain(idp,it)
enddo
enddo
!write(*,*) 'somewhere in the middle'



do it = 1,nits
x_its(it) = it
enddo

!write(*,*) 'somewhere down here'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ier = pgbeg(0,trim(adjustl(fname))//'/CPS',1,npar)
!!ier = pgbeg(0,'/XServe',1,npar)

!write(*,*) 'fsdfsdfsfdsfsfsdfs',npar
do ip = 1,npar
if (ymin_in(ip) .eq. 666.) then
ymin = minval(par(1:nits,ip))
else
ymin = ymin_in(ip)
endif

if (ymax_in(ip) .eq. 666.) then
ymax = maxval(par(1:nits,ip))
else if (ymax_in(ip) .eq. 555.) then
a = rms(par(1:nits,ip),nits)
ymax = avg(par(1:nits,ip),nits) + 5*a
else
ymax = ymax_in(ip)
endif

yrange = ymax - ymin

xmin = 1.0
xmax = nits

!call pgenv(xmin,xmax,ymin,ymax,0,0)
!call pgline(nits,x_its,par(1:nits,ip))

if (ip == npar) then
!call pglab('Iteration',trim(adjustl(axlab(ip))),'')
else
!call pglab('',trim(adjustl(axlab(ip))),'')
endif

enddo


!call pgend





end subroutine



!shade region between two curves

subroutine mypgshade(x,y1,y2,N)
integer N
real x(N), y1(N), y2(N),xp(4),yp(4)


!call pgsfs(3)
do i = 1,N-1
xp(1:2) = x(i:i+1)
xp(3:4) = xp(1:2)
yp(1:2) = y1(i:i+1)
yp(3:4) = y2(i:i+1)

!call pgpoly(4, xp, yp)
enddo

end subroutine





!! turn pgplot colours into wavelength colours for pgsci


subroutine wavcol(wav,icol)

!dwav = 800

if (wav .le. 2800) then
icol = 1
else if (wav .le. 3500) then
icol = 6
else if (wav .le. 4200) then
icol = 4
else if (wav .le. 4900) then
icol = 5
else if (wav .le. 5600) then
icol = 3
else if (wav .le. 6300) then
icol = 7
else if (wav .le. 7000) then
icol = 8
else
icol = 2
endif

end subroutine




!!! given a list return another list removing repeated elements e.g 1,2,2,3,4,3,3,5  --> 1,2,4,3,5
! input dat(n)
! op datop(n), nophi[index up to which uniqe elements exist in datop]

subroutine unique(dat,n,datop,nophi)

integer n,nophi
real dat(n), datop(n)

datop(2:n) = -666.0
datop(1) = dat(1)

idx_op = 2
do i = 2,n
datnow = dat(i)
iskip = 0
j     = 1

do while (iskip .eq. 0)
if (datop(j) == datnow) iskip = 1
if (j .eq. n) then
datop(idx_op) = datnow
idx_op = idx_op + 1
exit
endif
j = j + 1
enddo
enddo
nophi = idx_op - 1

end subroutine


!works 12 jun 15
!program testunitque
!
!real dat(10),datop(10)
!
!do i = 1,10
!dat(i) = i
!enddo
!
!
!dat(3) = 1
!dat(7) = 9
!
!
!
!n = 10
!
!call unique(dat,n,datop,nhi)
!
!write(*,*) dat
!write(*,*) datop
!write(*,*) nhi
!
!end program




!! subroutine to plot the fvg diagram within CREAM
!! input a file
subroutine cream_fvgplot(Nwavold, Ndatmaxold, t_dattemp_old, datinold, sigdatinin,&
taumean, sigtaumean, wavold,sigexpand, idxwav0,&
cosinc, sigcosinc, idxlc_old, ebmv_mw, z, omega_m, omega_l)




integer, intent(in):: Nwavold, Ndatmaxold, idxlc_old(Nwavold),idxwav0
real, intent(in):: datinold(Ndatmaxold,Nwavold), sigdatinin(Ndatmaxold,Nwavold), taumean(Nwavold), sigtaumean(Nwavold),&
wavold(Nwavold),cosinc, sigcosinc, t_dattemp_old(Ndatmaxold,Nwavold),sigexpand(Nwavold)
real wavlog(2), fnulog(2)
character(100):: ttemp,txlab,tx2lab
real xw2(2), yw2(2),sigdatin_old(Ndatmaxold,Nwavold)
integer,allocatable:: idxlc(:)
real,allocatable:: wavmod(:),fmodred(:),sigfmodred_log(:),fmod_unred_max(:),fmod_unred_min(:),&
t_datin(:,:), datin(:,:), sigdatin(:,:), wav(:),dat(:,:),sigdat(:,:),&
sigfmin(:), sigfmax(:),alam(:), sigalam(:),&
clam(:), sigclam(:),fmin(:), fmax(:),galmin(:), galmax(:),&
siggalmin(:), siggalmax(:)


!write(*,*) 'starting fvg plot'
wav0 = wavold(idxwav0)
!write(*,*) 'starting fvg plot'
rln10 = alog(10.)



!! expand data by error bars before proceeding new feature of 8 dec 2015
!write(*,*) 'starting fvg plot'
!write(*,*) Nwavold,idxwavnow_old
!write(*,*) 'starting fvg plot'

do i = 1,Nwavold
senow = sigexpand(i)
idxwavnow_old = idxlc_old(i)
coef = max(1.0,senow)
!write(*,*) 'here',i,Nwavold,idxwavnow_old,senow,coef,Ndatmaxold
sigdatin_old(1:idxwavnow_old,i) = coef * sigdatinin(1:idxwavnow_old,i)

enddo
!!
!stop
open(unit = 1,file = 'cream_fvgplot_8decpy.dat')
do i = 1,Nwavold
do it = 1,idxlc_old(i)
write(1,*) wavold(i),t_dattemp_old(it,i), datinold(it,i), sigdatin_old(it,i)
enddo
enddo
!! now call python script to combine data at same wavelength
call system('cp ../../../cream_fvgplot_8decpy.py ./')
call system('python cream_fvgplot_8decpy.py')
call system('rm cream_fvgplot_8decpy.py')


open(unit = 1,file='tempy_fvg.dat')
read(1,*) ndatmax,nwav
allocate(idxlc(nwav), datin(ndatmax,nwav), t_datin(ndatmax,nwav),&
sigdatin(ndatmax,nwav),wav(nwav),sigdat(Ndatmax,nwav), dat(ndatmax,nwav),&
sigfmin(nwav), sigfmax(nwav),alam(nwav), sigalam(nwav),&
clam(nwav), sigclam(nwav),fmin(nwav), fmax(nwav),galmin(nwav), galmax(nwav),&
siggalmin(nwav), siggalmax(nwav))



read(1,*) wav(1:nwav)
read(1,*) idxlc(1:nwav)
do ilc =1,nwav
do it = 1,idxlc(ilc)
read(1,*) crap, t_datin(it,ilc), datin(it,ilc), sigdatin(it,ilc)
enddo
enddo
close(1)







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !! deredden the light curves for MW extinction !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sigdat(1:Ndatmax,1:Nwav) = 0
dat(1:Ndatmax,1:Nwav) = 0

if (ebmv_mw .gt. 0) then
do ilc = 1, Nwav
do it = 1, idxlc(ilc)
call dustMW_sig(datin(it,ilc),sigdatin(it,ilc),wav(ilc),ebmv_mw,&
dat(it,ilc),sigdat(it,ilc))
enddo
enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!  calculate constant and variable components

call fvg_aug17(dat, sigdat, t_datin, Ndatmax, wav, Nwav, idxlc, alam, clam, sigalam, sigclam,fmin,fmax,&
sigfmin,sigfmax,xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0,galmin,galmax,&
siggalmin,siggalmax)

!call fvg_aug28(dat, sigdat, t_dattemp, Ndatmax, wav, Nwav, idxlc, alam, clam, sigalam, sigclam,fmin,fmax,&
!sigfmin,sigfmax,xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0,galmin,galmax,&
!siggalmin,siggalmax)
!!



!!

!! de-redden the variable disc spectrum
call dered_aug28(Nwav, wav, wav0, alam, sigalam, f0, sigf0, ebmv, sigebmv,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! New 19th August - to get luminosity distance, de-redden both the faint and bright spectrum!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call dered_aug28(Nwav, wav, wav0, fmin, sigfmin, f0_min, sigf0_min, ebmv_min, sigebmv_min,1)
call dered_aug28(Nwav, wav, wav0, fmax, sigfmax, f0_max, sigf0_max, ebmv_max, sigebmv_max,1)







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*) 'after dered'

wavlo = max(10.,minval(wav) - 100.)
wavhi = maxval(wav) + 1000.
wavloglo = alog10(wavlo)
wavloghi = alog10(wavhi)

wavlog(1) = wavloglo
wavlog(2) = wavloghi
fnulog(1) = alog10(f0) - 0.3333*alog10(wavlo/wav0)
fnulog(2) = alog10(f0) - 0.3333*alog10(wavhi/wav0)
f0log = alog10(f0)
sigf0log = 1./(alog(10.)*f0)*sigf0

Nmod = 100
dwmod = (wavhi-wavlo)/(Nmod - 1)


allocate(wavmod(Nmod), fmodred(Nmod), sigfmodred_log(Nmod),fmod_unred_max(Nmod), fmod_unred_min(Nmod))
do it = 1,Nmod

wm = wavlo + dwmod*(it-1)
ext = extmag_agn(wm,1.0)
wavmod(it) = wm
a0 = (wm/wav0)**(-1./3)
a = f0_max*a0
siga = a0*sigf0_max
b = 10**(-0.4*ext*ebmv_max)
sigb = -0.4*ext*rln10*b*sigebmv_max
fm = a*b
fmodred(it) = fm
sigfmodred_log(it) = 1./(rln10*fm)* fm*sqrt((siga/a)**2 + (sigb/b)**2)
enddo
!!

!do ilc = 1,Nwav
!write(*,*) wav(ilc),clam(ilc),sigclam(ilc),alam(ilc),sigalam(ilc), ebmv
!enddo


!stop

!!!! write fmin, sigfmin and fmax, sigfmax to separate text files to be fed into diffdisk_mcmc_2.f90 in fvg_nov14!!
open (unit = 1,file = 'diffdisk_mcmc_faint.dat')
do ilc = 1,Nwav
write(1,*) wav(ilc), fmin(ilc), sigfmin(ilc)
enddo
close(1)

open (unit = 1,file = 'diffdisk_mcmc_bright.dat')
do ilc = 1,Nwav
write(1,*) wav(ilc), fmax(ilc), sigfmax(ilc)
enddo
close(1)

open (unit = 1,file = 'diffdisk_mcmc_galmin.dat')
do ilc = 1,Nwav
write(1,*) wav(ilc), galmin(ilc), siggalmin(ilc)
enddo
close(1)
!!!!!




!ier = pgopen('cream_fvgplot_29aug.ps/CPS')
!call pgslw(4)
!call pgsch(1.6)
xmin = wavloglo
xmax = wavloghi
floor = 1.e-5
floorlog = alog10(floor)
ymin = alog10(min(minval(alam), minval(clam),minval(fmin))/1.1)
ymin = max(ymin,alog10(floor))
ymax = alog10(max(maxval(alam), maxval(clam),maxval(galmax),maxval(fmax))*1.1)

do ilc = 1,Nwav
if (galmin(ilc) .le. 0) galmin(ilc) = floor
enddo
!call pgsvp(0.15,0.9,0.15,0.9)
!call pgswin(xmin,xmax,ymin,ymax)


!call pgbox( 'bcnstl', 0., 0.0, 'bcnstl', 0., 0. )


do ilc = 1,Nwav
!call pgsci(4)
!call pgpt(1,alog10(wav(ilc)), max(alog10(galmin(ilc)),floorlog),6)					!! plot galin and galmax

if (galmin(ilc) .gt. 1.e-5) then
galpt = galmin(ilc)
else
galpt=galmin(ilc)+siggalmin(ilc)
endif
a1 = siggalmin(ilc)/(alog(10.)*galpt)

!write(*,*) ilc,galpt, galmin(ilc), siggalmin(ilc),a1,'sfsdfdsfds'

!if (galmin(ilc) .gt. 0) then
!call pgerrb(6, 1, alog10(wav(ilc)), alog10(galpt), a1, 2.0)
!endif
!!call pgpt(1,alog10(wav(ilc)), alog10(galmax(ilc)),6)
!a2 = siggalmax(ilc)/(alog(10.)*galmax(ilc))
!!call pgerrb(6, 1, alog10(wav(ilc)), alog10(galmax(ilc)), a2, 2.0)


!call pgsci(2)

!!call pgpt(1, alog10(wav(ilc)), alog10(alam(ilc)), 6)
!call pgpt(1, alog10(wav(ilc)), alog10(fmin(ilc)), 6)
!call pgpt(1, alog10(wav(ilc)), alog10(fmax(ilc)), 6)
!b = sigalam(ilc)/(alog(10.)*alam(ilc))
b1 = sigfmin(ilc)/(alog(10.)*fmin(ilc))
b2 = sigfmax(ilc)/(alog(10.)*fmax(ilc))
!!call pgerrb(6, 1, alog10(wav(ilc)), alog10(alam(ilc)), b, 2.0)
!call pgerrb(6, 1, alog10(wav(ilc)), alog10(fmin(ilc)), b1, 2.0)
!call pgerrb(6, 1, alog10(wav(ilc)), alog10(fmax(ilc)), b2, 2.0)



if (ilc .gt. 1) then  !! plot lines connecting points for clam and alam
xw2(1) = alog10(wav(ilc-1))
xw2(2) = alog10(wav(ilc))
!  yw2(1) = alog10(alam(ilc-1)) - b
!  yw2(2) = alog10(alam(ilc)) - b
!
!!  if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)
!  yw2(1) = alog10(alam(ilc-1)) + b
!  yw2(2) = alog10(alam(ilc)) + b
!!  if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)


yw2(1) = alog10(fmin(ilc-1)+sigfmin(ilc-1))
yw2(2) = alog10(fmin(ilc)+sigfmin(ilc))
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)



yw2(1) = alog10(fmin(ilc-1)-sigfmin(ilc-1))
yw2(2) = alog10(fmin(ilc)-sigfmin(ilc))
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)


yw2(1) = alog10(fmax(ilc-1)+sigfmax(ilc-1))
yw2(2) = alog10(fmax(ilc)+sigfmax(ilc))
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)



yw2(1) = alog10(fmax(ilc-1)-sigfmax(ilc-1))
yw2(2) = alog10(fmax(ilc)-sigfmax(ilc))
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)



!call pgsci(4)
yw2(1) = alog10(galmin(ilc-1)-siggalmin(ilc-1))! - a1
yw2(2) = alog10(galmin(ilc)-siggalmin(ilc))! - a1
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)



yw2(1) = alog10(galmin(ilc-1)+siggalmin(ilc-1))
yw2(2) = alog10(galmin(ilc)+siggalmin(ilc))
!if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)


!  yw2(1) = alog10(galmax(ilc-1)) - a2
!  yw2(2) = alog10(galmax(ilc)) - a2
!!  if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)
!
!  yw2(1) = alog10(galmax(ilc-1)) + a2
!  yw2(2) = alog10(galmax(ilc)) + a2
!!  if ((yw2(1) .gt. floorlog) .and. (yw2(2) .gt. floorlog)) call pgline(2,xw2,yw2)

endif
enddo


!! plot de-reddened disc slope

!!call pgsci(1)
!!call pgline(2,wavlog,fnulog)
!!call pgline(2,wavlog,fnulog+sigf0log)
!!call pgline(2,wavlog,fnulog-sigf0log)

!call pgsls(3)
!call pgsci(2)
!call pgline(Nmod,alog10(wavmod),alog10(fmodred))
!call pgline(Nmod,alog10(wavmod),alog10(fmodred)+sigfmodred_log)
!call pgline(Nmod,alog10(wavmod),alog10(fmodred)-sigfmodred_log)
!call pgsls(1)






!!!!!!!!!!! calculate faint de-reddened spectrum !!!!!!!!!!
do it = 1,Nmod

wm = wavmod(it)
ext = extmag_agn(wm,1.0)

a0 = (wm/wav0)**(-1./3)
a = f0_min*a0

fmod_unred_min(it) = a
fmod_unred_max(it) = f0_max*a0

siga = a0*sigf0_min
b = 10**(-0.4*ext*ebmv_min)
sigb = -0.4*ext*rln10*b*sigebmv_min
fm = a*b
fmodred(it) = fm
sigfmodred_log(it) = 1./(rln10*fm)* fm*sqrt((siga/a)**2 + (sigb/b)**2)
!write(*,*) wavmod(it), fmodred(it), alog10(fmodred(it)), sigfmodred_log(it)
enddo

!call pgsls(3)
!call pgsci(2)
!call pgline(Nmod,alog10(wavmod),alog10(fmodred))
!call pgline(Nmod,alog10(wavmod),alog10(fmodred)+sigfmodred_log)
!call pgline(Nmod,alog10(wavmod),alog10(fmodred)-sigfmodred_log)
!call pgsls(1)

!call pgline(Nmod,alog10(wavmod), alog10(fmod_unred_min))
!call pgline(Nmod,alog10(wavmod), alog10(fmod_unred_max))

!call pgsci(1)







a = 0.4*extmag_agn(wav0,1.0)
dered = 10**(a*ebmv)
!! calculate the Luminosity distance
disk_min = f0_min         !fmin(idxwav0)
disk_max_red = f0_max/dered     !fmax(idxwav0)
sig_diskmin = sigf0_min    !sigfmin(idxwav0)
sig_diskmax_red = sigf0_max/dered         !sigfmax(idxwav0)
tauave = taumean(idxwav0)
sigtauave = sigtaumean(idxwav0)


sigdered = dered*a*rln10*sigebmv
disk_max = dered*disk_max_red
sig_diskmax = disk_max*sqrt((sigdered/dered)**2 + &
(sig_diskmax_red/disk_max_red)**2)

eps      = disk_min/disk_max
sigeps   = eps* sqrt((sig_diskmin/disk_min)**2 + (sig_diskmax/disk_max)**2)

a_eps    = (1. - eps)/(1.-eps**(1.5))
a = (wav0/1.e4)**(-3./2)! * sqrt(40.1e3)
siga     = a_eps* sqrt( (sigeps/eps)**2 + (1.5*eps**0.5*sigeps)**2 )





c_f        = cosinc / disk_max
c_f_root   = sqrt(c_f)
sig_cf     = c_f* sqrt( (sigcosinc/cosinc)**2 + (sig_diskmax/disk_max)**2 )
sig_cfroot = 0.5/c_f_root*sig_cf


f0max_jy = f0_max*0.001
sig_f0maxjy = sigf0_max*0.001

root_f0max = sqrt(f0max_jy)
sig_root_f0max = 0.5/root_f0max * sig_f0maxjy

dlum = 6.3*a * tauave * c_f_root * a_eps /root_f0max
sigdlum = dlum*sqrt( (sigtauave/tauave)**2 + (sig_cfroot/c_f_root)**2 &
+ (sig_root_f0max/root_f0max)**2 + (siga/a_eps)**2)

!! annotate plot with luminosity distance and E(B-V)

! labels - ebmv
!call pgsci(1)
write(ttemp,'(F100.3)') ebmv ! legend
ttemp = adjustl(ttemp)
txlab = trim(ttemp)
write(ttemp,'(F100.3)') sigebmv
!write(*,*) ebmv(ilc), sigebmv
ttemp = adjustl(ttemp)
!write(*,*) trim(ttemp), ebmv,sigebmv,'sfdjr'
tx2lab = 'E(B-V) = '//trim(txlab)//'\(2233)'//trim(ttemp)
tx2lab = adjustl(tx2lab)
!call pgmtxt('LV', -1.0, 0.96,0.0, trim(tx2lab))


! labels - dl
!call pgsci(1)
write(ttemp,'(F100.2)') dlum ! legend
ttemp = adjustl(ttemp)
txlab = trim(ttemp)
write(ttemp,'(F100.2)') sigdlum
ttemp = adjustl(ttemp)
tx2lab = 'D\dL\u (Mpc) = '//trim(txlab)//'\(2233)'//trim(ttemp)
tx2lab = adjustl(tx2lab)
!call pgmtxt('LV', -1.0, 0.88,0.0, trim(tx2lab))




!! calculate the h0 value
a = zlum_ml( z, omega_m, omega_l )
oH = cgs('C')*a/(dlum * 100000)
sigoh = oH/dlum * sigdlum

write(ttemp,'(F100.2)') oH! legend
ttemp = adjustl(ttemp)
txlab = trim(ttemp)
write(ttemp,'(F100.2)') sigoh
ttemp = adjustl(ttemp)
tx2lab = 'H\d0\u (kms\u-1\dMpc\u-1\d) = '//trim(txlab)//'\(2233)'//trim(ttemp)
tx2lab = adjustl(tx2lab)
!call pgmtxt('LV', -1.0, 0.80,0.0, trim(tx2lab))

!dl=cgs('C')*zlum_ml( z, omega_m, omega_l )/Ho/100000 !! dl in Mpc (100000 is as ho should be in right units)


!call PGMTXT ('B', 3.0, 0.5, 0.5, 'Wavelength \A')
!call PGMTXT ('L', 3.0, 0.5, 0.5, 'F\d\gn\u(\gl) (mJy)')


!!! annotate plot with legend !!
! dlen = 0.1 !length of line to plot in terms of length f x axis
! uxstart = 0.05
! xrange = xmax - xmin
! Nlab = 4
!
!
! xlstart = xmin + uxstart*xrange
! xlen = xlstart + dlen*xrannge
!xw2(1) = xlstart
!xw2(2) = xlen
!
!!call pgsci(1)
!!call pgsls(1)



!call pgend

write(*,*) 'cream_fvgplot sumary:'
write(*,*) 'DL=',dlum, sigdlum
write(*,*) 'E(B-V)=',ebmv,sigebmv, ebmv_min,sigebmv_min,ebmv_max,sigebmv_max
write(*,*) 'wavarb, f0, sigf0=',wav0,f0,sigf0, f0_min, sigf0_min, f0_max, sigf0_max
write(*,*) 'epsilon correction=',eps
write(*,*) 'diskmin &max=',disk_min,disk_max
end subroutine









!!!!!! program to test the cream_fvgplot routine
!gflags cream_fvgplot.f90 ../numline.for ../fvg_nov14/fvg_jul15.f90 ../fvg_nov14/dered_jul15.f90 ../minmax.f90 ../readfromfile.f90 ../extmag_agn.for ../mystats.f90 ../inverse.f90
!program fvgtest
!
!real, allocatable:: taumean(:), sigtaumean(:), wav(:), dat(:,:),time(:,:), sig(:,:)
!integer, allocatable:: idxlc(:)
!character(100):: fname
!Nwav = 19
!Ndatmax = 500
!idxwav0 = 10
!
!allocate(taumean(Nwav), sigtaumean(Nwav), dat(Ndatmax, Nwav), wav(Nwav), idxlc(Nwav),&
!time(Ndatmax,Nwav), sig(Ndatmax,Nwav))
!
!
!wav(1:Nwav) = (/1157.5, 1367.0, 1478.5, 1746.0, 2120.0, 2310.0, 2910.0, 3440.0, 3551.0,&
!4340.0, 4353.0,4686.0, 5430.0, 5477.0, 6231.0, 6349.0, 7625.0, 8797.0, 8931.0/)        !! input wavs
!
!taumean(1:Nwav) = (/0.37, 0.45, 0.50, 0.61, 0.79, 0.88, 1.19, 1.48, 1.54, 2.01, 2.02,&
!2.23, 2.71, 2.74, 3.25, 3.34, 4.24, 5.10, 5.19/)  		 !! tau mean
!
!sigtaumean(1:Nwav) = (/0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,&
!0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1/)		 !! sig tau
!cosinc = 1.0
!sigcosinc = 0.05	    											 !! inc and uncertainty (from mcmc)
!
!
!fname = './stormlc_jun28/cream_fvgtest.dat'
!
!call readfromfile(fname,Ndatmax,Nwav,idxlc,time,dat,sig)
!
!do ilc = 1,Nwav
!do it = 1,idxlc(ilc)
!write(*,*) dat(it,ilc), sig(it,ilc)
!if (sig(it,ilc) .eq. 0) then
!write(*,*) 'zero error', it, ilc
!stop
!endif
!enddo
!enddo
!
!
!
!
!!call cream_fvgplot(Nwav, Ndatmax, dat, sig, taumean, sigtaumean, wav, idxwav0,&
!!cosinc, sigcosinc, idxlc)
!
!write(*,*) 'Finished'
!
!end program
!



!!fit model fnu(ilc,t) = clam(ilc) + xl(it)*slam(ilc) to light curves
!! differs from previous versions. We define a time grid at 1 day spacing going from 1st to final data point.
!! if 2 points in grid, points combined and error barrs added in quadrature.
!! if no points in grid, L(t) set as half way between adjacent positions

!! Update 18/08/15
!! update if you subtract mean and divide by rms inside ios loop, bad things happen.
!! ...you get a spikey L(t) function which leads to a constant galaxy flux brighter than the data points
!! ...fixed by leaving L(t) alone inside loop. Subtract mean and divide by rms after iterations finished.
!!... then also have to modify C(ilc) --> C(ilc) + mean(L(t))*S_old(ilc), S(ilc) --> S_old(ilc)*RMS(L(t))



!!! input fnu(Ndatmax,NLC), tnu(Ndatmax,NLC), sigfnu(Ndatmax,NLC) the light curve fluxes,times and errors at wavlength ilc


!!! output slam(NLC), clam(NLC) the constant and variable back ground spectrum
!!!....xlmin(-clam(1)/slam(1) i.e xl(it) value at 0 flux for shortes wavelength. Use to get host galaxy contribution.



subroutine fvg_aug17(fnu,sigfnu,tnu, Ndatmax, wav,NLC, idxlc, slam, clam, sigslam, sigclam,fmin,fmax,&
sigfmin,sigfmax,xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0,galmin,diskbg,siggalmin,sigdiskbg)



integer,intent(in):: NLC, Ndatmax, idxlc(NLC)
real,intent(in):: fnu(Ndatmax,NLC),sigfnu(Ndatmax,NLC), tnu(Ndatmax,NLC),wav(NLC)
real,intent(out):: slam(NLC), clam(NLC),sigclam(NLC), sigslam(NLC),fmin(NLC),fmax(NLC),&
sigfmin(NLC), sigfmax(NLC), xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0, galmin(NLC),&
diskbg(NLC),siggalmin(NLC), sigdiskbg(NLC)
real xl(Ndatmax, NLC), sigxl(Ndatmax, NLC)
real,allocatable:: fnuref(:,:), sigfnuref(:,:),ymod(:),&
tgrid(:),fnugrid(:,:), sigfnugrid(:,:), xlgrid(:), sigxlgrid(:), xlgridsort(:)
integer,allocatable:: inumgrid(:,:), igriditp_lo(:), igriditp_hi(:),inum_tot(:),&
idxkey(:)
logical testplot/.true./, diagnostic/.false./,custom_lab/.False./,custom_col/.False./
real slamold(NLC), clamold(NLC),xw2(2),yw2(2)
integer pgopen, isciwav(NLC),icuscol(NLC)
character(10) txlab,ttemp, t0, t1, t2, t3, t4, t5, t6, t7, t8
character(10) cuslab(NLC)
character(1000) tat



icount_all = 0
do ilc = 1,NLC
icount_all = idxlc(ilc) + icount_all
enddo


clam(1:NLC) = 1.0
slam(1:NLC) = 1.0



!!!!! set custom annotations !!!!
if (custom_lab) then
cuslab(1)  = 'HST 1158\A'
cuslab(2)  = 'HST 1367\A'
cuslab(3)  = 'HST 1478\A'
cuslab(4)  = 'HST 1746\A'
cuslab(5)  = 'Swift UVW2'
cuslab(6)  = 'Swift UVM2'
cuslab(7)  = 'Swift UVW1'
cuslab(8)  = 'Swift U'
cuslab(9)  = 'u'
cuslab(10) = 'B'
cuslab(11) = 'Swift B'
cuslab(12) = 'g'
cuslab(13) = 'Swift V'
cuslab(14) = 'V'
cuslab(15) = 'r'
cuslab(16) = 'R'
cuslab(17) = 'i'
cuslab(18) = 'I'
cuslab(19) = 'z'
endif

if(custom_col) then
icuscol(1)  = 12
icuscol(2)  = 4
icuscol(3)  = 8
icuscol(4)  = 2
icuscol(5)  = 1
icuscol(6)  = 12
icuscol(7)  = 4
icuscol(8)  = 11
icuscol(9)  = 5
icuscol(10) = 10
icuscol(11) = 8
icuscol(12) = 13
icuscol(13) = 2
icuscol(14) = 6
icuscol(15) = 1
icuscol(16) = 12
icuscol(17) = 4
icuscol(18) = 2
icuscol(19) = 6
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set up time grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dtgrid = 1.0
tlo = 1.e9
thi = -10.0
do ilc = 1,NLC
do it = 1,idxlc(ilc)
tnow = tnu(it,ilc)
if (tnow .lt. tlo) tlo = tnow
if (tnow .gt. thi) thi = tnow
enddo
enddo
thi = thi + dtgrid/5
tlo = tlo - dtgrid/5



Ntgrid = (thi - tlo)/dtgrid + 1
allocate(tgrid(Ntgrid),fnugrid(Ntgrid,NLC),sigfnugrid(Ntgrid,NLC),inumgrid(Ntgrid,NLC),&
xlgrid(NTgrid), igriditp_hi(Ntgrid), igriditp_lo(Ntgrid),sigxlgrid(NTgrid),&
inum_tot(Ntgrid),idxkey(Ntgrid), xlgridsort(Ntgrid))


igriditp_lo(1:Ntgrid) = 0
igriditp_hi(1:ntgrid) = 0
xlgrid(1:Ntgrid) = 1.0


do it = 1,Ntgrid
tgrid(it) = tlo + (it - 1)*dtgrid
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! assign fnu and sigfnu points to correct bins !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do it = 1,Ntgrid
tgnow = tgrid(it)
tglo = max(tgrid(1),tgnow - dtgrid/2)
tghi = min(tgrid(Ntgrid),tgnow + dtgrid/2)

do ilc = 1,NLC
sum_fnugrid     = 0.d0
sum2_sigfnugrid = 0.d0

inbin= 0
do it2 = 1,idxlc(ilc)
tnow = tnu(it2,ilc)

if ((tnow .ge. tglo) .and. (tnow .lt. tghi)) then !! assign point to correct bin
sum_fnugrid = sum_fnugrid + fnu(it2,ilc)
sum2_sigfnugrid = sum2_sigfnugrid + sigfnu(it2,ilc)**2
inbin = inbin + 1
endif
enddo

if (inbin .gt. 0) then
fnugrid(it,ilc)    = sum_fnugrid/inbin
sigfnugrid(it,ilc) = sqrt(sum2_sigfnugrid/inbin)
inumgrid(it,ilc)   = inbin
else
fnugrid(it,ilc)    = 0
sigfnugrid(it,ilc) = 0
inumgrid(it,ilc)   = 0
endif

enddo


enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! grid to interpolate xl grid (at time bins with no data for any wavelengths  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do it = 1,Ntgrid
isum = 0
do ilc = 1,NLC
isum = isum + inumgrid(it,ilc)
enddo
inum_tot(it) = isum
enddo
if (inum_tot(1) .eq. 0) then
write(*,*) 'fvg_aug17.f90, no points in first bin something is wrong'
stop
endif
if (inum_tot(Ntgrid) .eq. 0) then
write(*,*) 'fvg_aug17.f90, no points in last bin something is wrong'
!stop
endif


do it =2,Ntgrid-1
isum = inum_tot(it)
if (isum .eq. 0) then
ifound = 0
do it2 = it-1,1,-1
do ilc = 1,NLC
if(inumgrid(it2,ilc) .gt. 0) then
ititplo = it2
ifound = 1
exit
endif
enddo
if (ifound .eq. 1) exit
enddo
igriditp_lo(it) = ititplo

ifound = 0
do it2 = it+1,NTgrid
do ilc = 1,NLC
if(inumgrid(it2,ilc) .gt. 0) then
ititphi = it2
ifound = 1
exit
endif
enddo
if (ifound .eq. 1) exit
enddo
igriditp_hi(it) = ititphi
endif


enddo








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! perform ios analysis on binned data !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nits = 1000



do iteration = 1,nits
!! begin iterated optimal scaling to find the slam, clam and xl parameters

!! in


do ilc = 1,NLC
slamold(ilc) = slam(ilc)
top = 0.d0
bot = 0.d0
do it = 1,Ntgrid
if (inumgrid(it,ilc) .gt. 0) then
top = top + (fnugrid(it,ilc) - clam(ilc)) * xlgrid(it)/sigfnugrid(it,ilc)**2
bot = bot + (xlgrid(it) / sigfnugrid(it,ilc))**2
endif

enddo

slam(ilc) = abs(top/bot)
sigslam(ilc) = sqrt(sqrt(1./bot))



enddo

!! Now clam

do ilc = 1,NLC
clamold(ilc) = clam(ilc)
top = 0.d0
bot = 0.d0
do it = 1,Ntgrid
!write(*,*) fnugrid(it,ilc), sigfnugrid(it,ilc), it, 'testing fit'
if (inumgrid(it,ilc) .gt. 0) then
top = top + (fnugrid(it,ilc) - slam(ilc) * xlgrid(it))/sigfnugrid(it,ilc)**2
bot = bot + (1. / sigfnugrid(it,ilc))**2
endif
enddo

clam(ilc) = top/bot
sigclam(ilc) = sqrt(sqrt(1./bot))
!write(*,*) ilc, clam(ilc), sigclam(ilc), iteration,'testing fit!!'
!read(*,*)
enddo


!! Now xl
do it = 1,Ntgrid
top = 0.d0
bot = 0.d0


do ilc = 1,NLC
if (inumgrid(it,ilc) .gt. 0) then
top = top + (fnugrid(it,ilc) - clam(ilc)) * slam(ilc)/(sigfnugrid(it,ilc))**2
bot = bot + (slam(ilc)/sigfnugrid(it,ilc))**2
endif
enddo
if (top .gt. 0) then
xlgrid(it) = top/bot
sigxlgrid(it) = 1./sqrt(real(bot))
endif

enddo




!!! Not all L(t) points will have any data in t (some bins may be empty)
! set L(t) = half way between adjacent points in this case!!
do it = 1,Ntgrid


if (inum_tot(it) .eq. 0) then
idxhi = igriditp_hi(it)
idxlo = igriditp_lo(it)

xlhi = xlgrid(idxhi)
xllo = xlgrid(idxlo)
xlgrid(it) = xllo + 0.5*(xlhi - xllo)
endif
enddo






! convergence test
idx_count = 0
frac_cutoff = 1.e-5
do ilc = 1,NLC
a = abs(clam(ilc) - clamold(ilc))/clamold(ilc)
b = abs(slam(ilc) - slamold(ilc))/slamold(ilc)
if ( (a .lt. frac_cutoff) .and. (b .lt. frac_cutoff) ) then
idx_count = idx_count + 1
endif
enddo

if (idx_count == NLC) then
write(*,*) 'fvg_jul15.f90: convergence after iteration',iteration
exit
endif


enddo !end iteration loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! sub mean and rms
sum = 0.d0
do it = 1,Ntgrid
sum = sum + xlgrid(it)
enddo
xlmean = sum/Ntgrid
sum = 0.d0
do it = 1,Ntgrid
a = xlgrid(it) - xlmean
sum = sum + a*a
enddo
xlrms = sqrt(sum/Ntgrid)

!write(*,*) 'xlrms',xlrms, xlmean, slam(1:5)
!stop

xlgrid(1:Ntgrid) = (xlgrid(1:Ntgrid) - xlmean)/xlrms
sigxlgrid(1:Ntgrid) = sigxlgrid(1:Ntgrid)/xlrms



clam(1:NLC) = clam(1:NLC) + xlmean*slam(1:NLC)



slam(1:NLC) = slam(1:NLC)*xlrms
sigslam(1:NLC) = sigslam(1:NLC)*xlrms

!write(*,*) xlrms
!write(*,*) slam(1:5)/xlrms
!stop





!!!!!!!!!!!!!!! save results for python program
!
!do ilc = 1,NLC
!open(unit = 1,file = 'pytest_inp.dat')
!do it = 1,Ntgrid
!write(1,*) xlgrid(it), fnugrid(it,ilc), sigfnugrid(it,ilc)
!enddo
!close(1)
!call system('python pytest_leastsq_fluxflux.py')
!write(*,*) slam(ilc), sigslam(ilc), clam(ilc), sigclam(ilc),'fortran'
!read(*,*)
!enddo
!
!
!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  The rest is identical to the old routine !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!! find the min and max value of xl
xlmin = minval(xlgrid)
xlmax = maxval(xlgrid)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! find the xlgalmin, xldiskbg, xldiskmax !!!!!!!!!!!!!!!
xlgalmin = -4.15!-clam(1)/slam(1)          !! xl at min gal flux


xldiskbg = 10000
do it = 1,Ntgrid
if (inum_tot(it) .gt. 0) then
if (xlgrid(it) .lt. xldiskbg) xldiskbg = xlgrid(it)
endif
enddo
!xldiskbg = minval(xlgrid)                !! xlat disk background

xldiskmax = -10
do it = 1,Ntgrid
if (inum_tot(it) .gt. 0) then
if (xlgrid(it) .gt. xldiskmax) xldiskmax = xlgrid(it)
endif
enddo
!xldiskmax = maxval(xlgrid)                   !! xlat disk background
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!XL0 = -clam(1)/slam(1) i.e when flux goes to zero for shortest wavelenght (gives you L(t) value for the host galaxy contribution).
tempmod = -clam(1)/slam(1)
sigtempmod = tempmod * sqrt( (sigclam(1)/clam(1))**2 + (sigclam(1)/slam(1))**2 )
xl0 = tempmod
sigxl0 = sigtempmod


!f = c +L*A

call sort1d(xlgrid,Ntgrid,xlgridsort,idxkey)

do ilc = 1,NLC

do idx = 1,Ntgrid
it = idxkey(idx)
if (fnugrid(it,ilc) .gt. 0) xldiskmax = xlgrid(it)
enddo

do idx = 1,Ntgrid
it = idxkey(idx)
if (fnugrid(it,ilc) .gt. 0) then
xldiskbg = xlgrid(it)
exit
endif
enddo

!write(*,*) ilc, xldiskbg, xldiskmax
!do it = 1,Ntgrid
! write(*,*) it, xlgrid(it)
!enddo
!read(*,*)

!! The galaxy component
galmin(ilc) = clam(ilc) + xlgalmin*slam(ilc)
!write(*,*) 'werewwt',wav(ilc), clam(ilc), sigclam(ilc), slam(ilc), sigslam(ilc)
!read(*,*)

siggalmin(ilc) = sqrt( sigclam(ilc)**2 + (xlgalmin*sigslam(ilc))**2)

!! The disk back ground
diskbg(ilc) = clam(ilc) + xldiskbg*slam(ilc)
sigdiskbg(ilc) = sqrt( sigclam(ilc)**2 + (xldiskbg*sigslam(ilc))**2)


!! the minimum and maximum accretion disc components
a = abs(xldiskbg - xlgalmin)
fmin(ilc) = a*slam(ilc)
sigfmin(ilc) =a*sigslam(ilc)

a = abs(xldiskmax - xlgalmin)
fmax(ilc) = a*slam(ilc)
sigfmax(ilc) = a * sigslam(ilc)
enddo




if (testplot) then

call minmax2d(fnugrid,Ntgrid,NLC,idxlo1,idxlo2,idxhi1,idxhi2)

xmin = min(xlmin, -clam(1)/slam(1))*1.1
xmax = 4.0
ymin = 0.1
roof = 1.e4
ymax = min(roof,fnugrid(idxhi1,idxhi2))*1.1

if (xmin .ne. xmin) xmin = -10.0


!!ier = pgopen('/Xserve')
!ier = pgopen('test_fluxflux_fvgaug17.ps/VCPS')
!call pgsvp(0.15,0.9,0.15,0.9)

!call pgslw(4)
!call pgsch(1.6)
xmax_for_ann = 6.0
!call pgswin(xmin,xmax_for_ann,ymin,ymax)
!call pgbox( 'bcnst', 0., 0.0, 'bcnst', 0., 0. )
!call pgmtxt('L',2.0,0.5,0.5,'f\d\gn\u(\gl) (mJy)')
!call PGMTXT ('B', 2.0, 0.5, 0.5, 'L(t)')


xpos = (xlgalmin -xmin)/(xmax_for_ann-xmin)
!call pgmtxt('T',0.8,xpos,0.5,'L\d gal \u')

xpos = (xldiskbg -xmin)/(xmax_for_ann-xmin)
!call pgmtxt('T',0.8,xpos,0.5,'L\d F \u')

xpos = (xldiskmax -xmin)/(xmax_for_ann-xmin)
!call pgmtxt('T',0.8,xpos,0.5,'L\d B \u')


!! plot vertical lines to show the xlmin and xlmax positions
yw2(1) = ymin
yw2(2) = ymax
xw2(1:2) = xlgalmin
!call pgline(2,xw2,yw2)
xw2(1:2) = xldiskbg
!call pgline(2,xw2,yw2)
xw2(1:2) = xldiskmax
!call pgline(2,xw2,yw2)


open(unit = 1, file = 'fluxflux_tab_aug17.dat')
do iw = 1,NLC
write(t0,'(I10.3)') int(wav(iw))
write(t1,'(F10.2)') galmin(iw)
write(t2,'(F10.2)') siggalmin(iw)
write(t3,'(F10.2)') fmin(iw)
write(t4,'(F10.2)') sigfmin(iw)
write(t5,'(F10.2)') fmax(iw)
write(t6,'(F10.2)') sigfmax(iw)

write(*,*) galmin(iw), siggalmin(iw), fmin(iw), sigfmin(iw), fmax(iw), sigfmin(iw)
!read(*,*)

b = fmax(iw)/fmin(iw)
sigb = b*sqrt( (sigfmin(iw)/fmin(iw))**2 + (sigfmax(iw)/fmax(iw))**2 )

write(t7,'(F10.2)') b
write(t8,'(F10.2)') sigb

tat = '$'//trim(adjustl(t0))//'$'
tat = trim(adjustl(tat))//' & $'//trim(adjustl(t1))//'\pm'//trim(adjustl(t2))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t3))//'\pm'//trim(adjustl(t4))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t5))//'\pm'//trim(adjustl(t6))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t7))//'\pm'//trim(adjustl(t8))
tat = trim(adjustl(tat))//'$ //'

write(1,*) trim(adjustl(tat))
enddo
close(1)



!! obtain correct colours for lines
call sciwav(wav, isciwav, NLC)


do ilc = 1,NLC

iclcnow = isciwav(ilc)
!call wavcol(wav(ilc),iclcnow)

if (custom_col) then
!call pgsci(icuscol(ilc))
else
!call pgsci(iclcnow)
endif


if (custom_lab) then
txlab=trim(adjustl(cuslab(ilc)))
else
write(ttemp,'(I10.3)') int(wav(ilc)) ! legend
ttemp = adjustl(ttemp)
txlab = trim(ttemp)//'\A'
txlab = adjustl(txlab)
endif
!write(*,*) wav(ilc), trim(adjustl(ttemp)),'fdfdf',iclcnow
!call pgmtxt('RV', -1.0, 0.01 + 1.*ilc/NLC*0.95 ,1.0, trim(txlab))


!do it = 1,Ntgrid
!write(*,*) ilc, it, xlgrid(it), fnugrid(it,ilc), 'arararararara'
!enddo

!call pgpt(Ntgrid,xlgrid,fnugrid(1:Ntgrid,ilc),6)
!call pgerrb(6,Ntgrid,xlgrid,fnugrid(1:ntgrid,ilc),sigfnugrid(1:ntgrid,ilc),2.0)
xw2(1) = xlgalmin
xw2(2) = xldiskmax
yw2(1:2) = xw2(1:2)*slam(ilc) + clam(ilc)



!call pgline(2,xw2,yw2)





enddo
!call pgend
endif



write(*,*) 'slam and uncerts...'
write(*,*) slam
write(*,*) sigslam
write(*,*) 'clam and uncerts...'
write(*,*) clam
write(*,*) sigclam
write(*,*)

open(unit = 1,file = 'fvg_aug17_clamslam.dat')
do ilc = 1,NLC
write(t1,'(F10.2)') clam(ilc)
write(t2,'(F10.2)') sigclam(ilc)
write(t3,'(F10.2)') slam(ilc)
write(t4,'(F10.2)') sigslam(ilc)

write(1,*) '$'//trim(adjustl(t1))//'\pm',trim(adjustl(t2))//'$ & ',&
'$'//trim(adjustl(t3))//'\pm',trim(adjustl(t4))//'$ '
enddo
close(1)


!!! diagnostic plot !!!
if (diagnostic) then
do ilc = 1,NLC
!ier = pgopen('/Xserve')
!call pgsvp(0.1,0.9,0.5,0.9)
xmin = tgrid(1)
xmax = tgrid(Ntgrid)
ymin = minval(xlgrid)
ymax = maxval(xlgrid)
!call pgswin(xmin,xmax,ymin,ymax)
!call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!call pgline(Ntgrid,tgrid,xlgrid)

!call pgsvp(0.1,0.9,0.1,0.4)
xmin = tgrid(1)
xmax = tgrid(Ntgrid)
allocate(ymod(Ntgrid))

ymod(1:Ntgrid) = clam(ilc) + xlgrid(1:Ntgrid)*slam(ilc)
ymin = 0.0!minval(ymod)
ymax = maxval(ymod)

!call pgswin(xmin,xmax,ymin,ymax)
!call pgline(Ntgrid,tgrid,ymod)
!call pgbox('bcnst',0.0,0,'bcnst',0.0,0)

!call pgpt(idxlc(ilc),tnu(1:idxlc(ilc),ilc),fnu(1:idxlc(ilc),ilc),6)
!call pgerrb(6,idxlc(ilc),tnu(1:idxlc(ilc),ilc),fnu(1:idxlc(ilc),ilc),&
!sigfnu(1:idxlc(ilc),ilc), 1.0)

!call pgsci(2)
do it = 1,Ntgrid
if (inum_tot(it) .eq. 0) then
!call pgsci(4)
else if(inumgrid(it,ilc) .eq. 0) then
!call pgsci(3)
else
!call pgsci(2)
endif
!call pgpt(1,tgrid(it),fnugrid(it,ilc),6)
!call pgerrb(6,1,tgrid(it),fnugrid(it,ilc),&
!sigfnugrid(it,ilc), 1.0)
enddo

!call pgsci(1)
deallocate(ymod)
read(*,*)
!call pgend
enddo
endif








write(*,*) 'End of routine fvg_aug17.f90'

end subroutine













!!fit model fnu(ilc,t) = clam(ilc) + xl(it)*slam(ilc) to light curves
!! differs from previous versions. We define a time grid at 1 day spacing going from 1st to final data point.
!! if 2 points in grid, points combined and error barrs added in quadrature.
!! if no points in grid, L(t) set as half way between adjacent positions

!! Update 18/08/15
!! update if you subtract mean and divide by rms inside ios loop, bad things happen.
!! ...you get a spikey L(t) function which leads to a constant galaxy flux brighter than the data points
!! ...fixed by leaving L(t) alone inside loop. Subtract mean and divide by rms after iterations finished.
!!... then also have to modify C(ilc) --> C(ilc) + mean(L(t))*S_old(ilc), S(ilc) --> S_old(ilc)*RMS(L(t))



!!! input fnu(Ndatmax,NLC), tnu(Ndatmax,NLC), sigfnu(Ndatmax,NLC) the light curve fluxes,times and errors at wavlength ilc


!!! output slam(NLC), clam(NLC) the constant and variable back ground spectrum
!!!....xlmin(-clam(1)/slam(1) i.e xl(it) value at 0 flux for shortes wavelength. Use to get host galaxy contribution.



subroutine fvg_aug28(fnu,sigfnu,tnu, Ndatmax, wav,NLC, idxlc, slam, clam, sigslam, sigclam,fmin,fmax,&
sigfmin,sigfmax,xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0,galmin,diskbg,siggalmin,sigdiskbg)



integer,intent(in):: NLC, Ndatmax, idxlc(NLC)
real,intent(in):: fnu(Ndatmax,NLC),sigfnu(Ndatmax,NLC), tnu(Ndatmax,NLC),wav(NLC)
real,intent(out):: slam(NLC), clam(NLC),sigclam(NLC), sigslam(NLC),fmin(NLC),fmax(NLC),&
sigfmin(NLC), sigfmax(NLC), xlmin,xlmax,xl0,sigxlmin,sigxlmax,sigxl0, galmin(NLC),&
diskbg(NLC),siggalmin(NLC), sigdiskbg(NLC)
real xl(Ndatmax, NLC), sigxl(Ndatmax, NLC)
real,allocatable:: fnuref(:,:), sigfnuref(:,:),ymod(:),&
tgrid(:),fnugrid(:,:), sigfnugrid(:,:), xlgrid(:), sigxlgrid(:),pat(:,:),pout(:),sigpout(:)
integer,allocatable:: inumgrid(:,:), igriditp_lo(:), igriditp_hi(:),inum_tot(:)
logical testplot/.true./, diagnostic/.false./
real slamold(NLC), clamold(NLC),xw2(2),yw2(2)
integer pgopen, isciwav(NLC)
character(10) txlab,ttemp,t0,t1,t2,t3,t4,t5,t6,t7,t8
character(1000) tat



icount_all = 0
do ilc = 1,NLC
icount_all = idxlc(ilc) + icount_all
enddo


clam(1:NLC) = 1.0
slam(1:NLC) = 1.0




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set up time grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dtgrid = 1.0
tlo = 1.e9
thi = -10.0
do ilc = 1,NLC
do it = 1,idxlc(ilc)
tnow = tnu(it,ilc)
if (tnow .lt. tlo) tlo = tnow
if (tnow .gt. thi) thi = tnow
enddo
enddo
thi = thi + dtgrid/5
tlo = tlo - dtgrid/5




Ntgrid = (thi - tlo)/dtgrid + 1
allocate(tgrid(Ntgrid),fnugrid(Ntgrid,NLC),sigfnugrid(Ntgrid,NLC),inumgrid(Ntgrid,NLC),&
xlgrid(NTgrid), igriditp_hi(Ntgrid), igriditp_lo(Ntgrid),sigxlgrid(NTgrid),inum_tot(Ntgrid))
igriditp_lo(1:Ntgrid) = 0
igriditp_hi(1:ntgrid) = 0
xlgrid(1:Ntgrid) = 1.0


do it = 1,Ntgrid
tgrid(it) = tlo + (it - 1)*dtgrid
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! assign fnu and sigfnu points to correct bins !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do it = 1,Ntgrid
!write(*,*) it,tgnow,ilc,Ntgrid,NLC
tgnow = tgrid(it)
tglo = max(tgrid(1),tgnow - dtgrid/2)
tghi = min(tgrid(Ntgrid),tgnow + dtgrid/2)

do ilc = 1,NLC
sum_fnugrid     = 0.d0
sum2_sigfnugrid = 0.d0


inbin= 0
do it2 = 1,idxlc(ilc)
tnow = tnu(it2,ilc)

if ((tnow .ge. tglo) .and. (tnow .lt. tghi)) then !! assign point to correct bin
sum_fnugrid = sum_fnugrid + fnu(it2,ilc)
sum2_sigfnugrid = sum2_sigfnugrid + sigfnu(it2,ilc)**2
inbin = inbin + 1
!write(*,*) ilc, tnow, tglo, tghi,'dfdfdfd'
endif
enddo

if (inbin .gt. 0) then
fnugrid(it,ilc)    = sum_fnugrid/inbin
sigfnugrid(it,ilc) = sqrt(sum2_sigfnugrid/inbin)
inumgrid(it,ilc)   = inbin
else
fnugrid(it,ilc)    = 0
sigfnugrid(it,ilc) = 0
inumgrid(it,ilc)   = 0
endif

enddo


enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! grid to interpolate xl grid (at time bins with no data for any wavelengths  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do it = 1,Ntgrid
isum = 0
do ilc = 1,NLC
isum = isum + inumgrid(it,ilc)
enddo
inum_tot(it) = isum
enddo
if (inum_tot(1) .eq. 0) then
write(*,*) 'fvg_aug28.f90, no points in first bin something is wrong'
stop
endif
if (inum_tot(Ntgrid) .eq. 0) then
write(*,*) 'fvg_aug28.f90, no points in last bin something is wrong'
!stop
endif


do it =2,Ntgrid-1
isum = inum_tot(it)
if (isum .eq. 0) then
ifound = 0
do it2 = it-1,1,-1
do ilc = 1,NLC
if(inumgrid(it2,ilc) .gt. 0) then
ititplo = it2
ifound = 1
exit
endif
enddo
if (ifound .eq. 1) exit
enddo
igriditp_lo(it) = ititplo

ifound = 0
do it2 = it+1,NTgrid
do ilc = 1,NLC
if(inumgrid(it2,ilc) .gt. 0) then
ititphi = it2
ifound = 1
exit
endif
enddo
if (ifound .eq. 1) exit
enddo
igriditp_hi(it) = ititphi
endif


enddo

!call chdir('../../../fvg_nov14/')
!open(unit=1,file='testdump.dat')
!write(1,*) Ndatmax, NLC
!write(1,*) idxlc
!do it = 1,Ndatmax
!write(1,*) tnu(it,1:NLC), fnu(it,1:NLC), sigfnu(it,1:NLC)
!enddo
!close(1)
!stop
Npat = 2
allocate(pat(Ntgrid,Npat), pout(Npat), sigpout(Npat))
pat(1:Ntgrid,2) = 1.0
xlglo = -1.0!-5.0
xlghi = 1.0!5.0
dxlg  = (xlghi - xlglo)/(Ntgrid-1)
do it = 1,Ntgrid
xlgrid(it) = xlglo + (it-1)*dxlg
enddo

nits = 1000
do iteration = 1, nits


!!!! use hes fit to fit the slam and c lam components of the model

pat(1:Ntgrid,1) = xlgrid

!open(unit=1,file='testdump.dat')
!write(1,*) Ntgrid,NLC
!do it = 1,Ntgrid
!write(1,*) xlgrid(it) , fnugrid(it,ilc), sigfnugrid(it,ilc)
!enddo
!close(1)
!stop

do  ilc = 1,NLC
!write(*,*) 'before hesfit fvg_aug28'
call hesfit2(xlgrid,fnugrid(1:Ntgrid,ilc),sigfnugrid(1:Ntgrid,ilc),pat,Ntgrid,NPat,pout,sigpout,0.0) !! fit use lleast squares and ignore grid elements with no points in the bin (i.e fnugrid = 0)
!write(*,*) 'after hesfit',iteration,slam,clam
slam(ilc)    = pout(1)
sigslam(ilc) = sigpout(1)
clam(ilc)    = pout(2)
sigclam(ilc) = sigpout(2)
enddo


!!! Now use an ios routine to optimise xlgrid

!! Now xl

!write(*,*) xlgrid(1), xlgrid(10), xlgrid(Ntgrid), iteration
do it = 1,Ntgrid
top = 0.d0
bot = 0.d0


do ilc = 1,NLC
if (inumgrid(it,ilc) .gt. 0) then
top = top + (fnugrid(it,ilc) - clam(ilc)) * slam(ilc)/(sigfnugrid(it,ilc))**2
bot = bot + (slam(ilc)/sigfnugrid(it,ilc))**2
endif
enddo
if (top .gt. 0) then
xlgrid(it) = top/bot
sigxlgrid(it) = 1./sqrt(real(bot))
endif

enddo





!!! Not all L(t) points will have any data in t (some bins may be empty)
!! set L(t) = half way between adjacent points in this case!!
do it = 1,Ntgrid


if (inum_tot(it) .eq. 0) then
idxhi = igriditp_hi(it)
idxlo = igriditp_lo(it)

xlhi = xlgrid(idxhi)
xllo = xlgrid(idxlo)
xlgrid(it) = xllo + 0.5*(xlhi - xllo)
endif
enddo




! convergence test
idx_count = 0
frac_cutoff = 1.e-4
do ilc = 1,NLC
a = abs(clam(ilc) - clamold(ilc))/clamold(ilc)
b = abs(slam(ilc) - slamold(ilc))/slamold(ilc)
if ( (a .lt. frac_cutoff) .and. (b .lt. frac_cutoff) ) then
idx_count = idx_count + 1
endif
enddo

if (idx_count == NLC) then
write(*,*) 'fvg_jul15.f90: convergence after iteration',iteration
exit
endif



!do it = 1,Ntgrid
!write(*,*) xlgrid(it), xlmean, xlrms, 'dffd', iteration
!enddo


enddo !end iteration loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! sub mean and rms
sum = 0.d0
do it = 1,Ntgrid
sum = sum + xlgrid(it)
enddo
xlmean = sum/Ntgrid
sum = 0.d0
do it = 1,Ntgrid
a = xlgrid(it) - xlmean
sum = sum + a*a
enddo
xlrms = sqrt(sum/Ntgrid)

xlgrid(1:Ntgrid) = (xlgrid(1:Ntgrid) - xlmean)/xlrms


sigxlgrid(1:Ntgrid) = sigxlgrid(1:Ntgrid)/xlrms

clam(1:NLC) = clam(1:NLC) + xlmean*slam(1:NLC)



slam(1:NLC) = slam(1:NLC)*xlrms
sigslam(1:NLC) = sigslam(1:NLC)*xlrms









open(unit = 1,file = 'fvg_aug28_clamslam.dat')
do ilc = 1,NLC
write(t1,'(F10.2)') clam(ilc)
write(t2,'(F10.2)') sigclam(ilc)
write(t3,'(F10.2)') slam(ilc)
write(t4,'(F10.2)') sigslam(ilc)

write(1,*) '$'//trim(adjustl(t1))//'\pm',trim(adjustl(t2))//'$ & ',&
'$'//trim(adjustl(t3))//'\pm',trim(adjustl(t4))//'$ '
enddo
close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  The rest is identical to the old routine !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!! find the min and max value of xl
xlmin = minval(xlgrid)
xlmax = maxval(xlgrid)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! find the xlgalmin, xldiskbg, xldiskmax !!!!!!!!!!!!!!!

fmin_spec = 1.e6
do ilc = 1,NLC
a = clam(ilc) + xlmin*slam(ilc)
if (a .lt. fmin_spec) then
fmin_spec = a
ilclo = ilc
endif
enddo


xlgalmin = -clam(ilclo)/slam(ilclo)          !! xl at min gal flux

!call minmax2d(fnugrid,Ntgrid,NLC,idxlo1,idxlo2,idxhi1,idxhi2)

xldiskbg = 10000
do it = 1,Ntgrid
if (inum_tot(it) .gt. 0) then
if (xlgrid(it) .lt. xldiskbg) xldiskbg = xlgrid(it)
endif
enddo
!xldiskbg = minval(xlgrid)                !! xlat disk background

xldiskmax = -10
do it = 1,Ntgrid
if (inum_tot(it) .gt. 0) then
if (xlgrid(it) .gt. xldiskmax) xldiskmax = xlgrid(it)
endif
enddo
!xldiskmax = maxval(xlgrid)                   !! xlat disk background
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!XL0 = -clam(1)/slam(1) i.e when flux goes to zero for shortest wavelenght (gives you L(t) value for the host galaxy contribution).
tempmod = -clam(ilclo)/slam(ilclo)
sigtempmod = tempmod * sqrt( (sigclam(ilclo)/clam(ilclo))**2 + (sigclam(ilclo)/slam(ilclo))**2 )
xl0 = tempmod
sigxl0 = sigtempmod


!f = c +L*A

do ilc = 1,NLC

do it = 1,Ntgrid
if (fnugrid(it,ilc) .gt. 0) xldiskmax = xlgrid(it)
enddo

do it = 1,Ntgrid
if (fnugrid(it,ilc) .gt. 0) then
xldiskbg = xlgrid(it)
exit
endif
enddo



!! The galaxy component
galmin(ilc) = clam(ilc) + xlgalmin*slam(ilc)
siggalmin(ilc) = sqrt( sigclam(ilc)**2 + (xlgalmin*sigslam(ilc))**2)
!! The disk back ground
diskbg(ilc) = clam(ilc) + xldiskbg*slam(ilc)
sigdiskbg(ilc) = sqrt( sigclam(ilc)**2 + (xldiskbg*sigslam(ilc))**2)



!! the minimum and maximum accretion disc components
a = abs(xldiskbg - xlgalmin)


fmin(ilc) = a*slam(ilc)
sigfmin(ilc) =a*sigslam(ilc)

a = abs(xldiskmax - xlgalmin)
fmax(ilc) = a*slam(ilc)
sigfmax(ilc) = a * sigslam(ilc)



!adim = maxval(fnugrid(1:Ntgrid,ilc))
!abri = adim
!
!do it = 1,Ntgrid
!a = fnugrid(it,ilc)
!if ((a .lt. adim) .and. (a .gt. 0)) adim = a
!enddo
!
!
!fmax(ilc) = abri
!fmin(ilc) = adim






enddo




if (testplot) then

call minmax2d(fnugrid,Ntgrid,NLC,idxlo1,idxlo2,idxhi1,idxhi2)

xmin = min(xlmin, -clam(1)/slam(1))*1.1
xmax = 4.0
ymin = 0.1
roof = 1.e4
ymax = min(roof,fnugrid(idxhi1,idxhi2))*1.1

write(*,*) 'making testplot fvg_jul24.f90'
!!ier = pgopen('/Xserve')
!ier = pgopen('test_fluxflux_fvgaug28.ps/CPS')
!call pgsvp(0.15,0.9,0.15,0.9)

!call pgslw(4)
!call pgsch(1.6)
!call pgswin(xmin,xmax,ymin,ymax)
!call pgbox( 'bcnst', 0., 0.0, 'bcnst', 0., 0. )
!call pgmtxt('L',2.0,0.5,0.5,'F\d\gn\u(\gl) (mJy)')
!call PGMTXT ('B', 2.0, 0.5, 0.5, 'L(t)')



!! plot vertical lines to show the xlmin and xlmax positions
yw2(1) = ymin
yw2(2) = ymax
xw2(1:2) = xlgalmin
!call pgline(2,xw2,yw2)
xw2(1:2) = xldiskbg
!call pgline(2,xw2,yw2)
xw2(1:2) = xldiskmax
!call pgline(2,xw2,yw2)


!write(*,*) xlgalmin,xldiskbg,xldiskmax
!write(*,*) xlmin, xlmax, ymin, ymax, 'fvg_aug xlmin/max'



!! obtain correct colours for lines
call sciwav(wav, isciwav, NLC)


do ilc = 1,NLC

iclcnow = isciwav(ilc)
!call wavcol(wav(ilc),iclcnow)
!call pgsci(iclcnow)



write(ttemp,'(I10.3)') int(wav(ilc)) ! legend
ttemp = adjustl(ttemp)
txlab = trim(ttemp)//'\A'
txlab = adjustl(txlab)

!write(*,*) wav(ilc), trim(adjustl(ttemp)),'fdfdf',iclcnow
!call pgmtxt('RV', -1.0, 0.01 + 1.*ilc/NLC*0.95 ,1.0, trim(txlab))



!call pgpt(Ntgrid,xlgrid,fnugrid(1:Ntgrid,ilc),6)
!call pgerrb(6,Ntgrid,xlgrid,fnugrid(1:ntgrid,ilc),sigfnugrid(1:ntgrid,ilc),2.0)
xw2(1) = xlgalmin
xw2(2) = xldiskmax
yw2(1:2) = xw2(1:2)*slam(ilc) + clam(ilc)



!call pgline(2,xw2,yw2)





enddo
!call pgend
endif


write(*,*) 'fvg_aug28.f90'
write(*,*) 'slam and uncerts...'
write(*,*) slam
write(*,*) sigslam
write(*,*) 'clam and uncerts...'
write(*,*) clam
write(*,*) sigclam
write(*,*) ''



open(unit = 1, file = 'fluxflux_tab_aug28.dat')
do iw = 1,NLC
write(t0,'(I10.3)') int(wav(iw))
write(t1,'(F10.2)') galmin(iw)
write(t2,'(F10.2)') siggalmin(iw)
write(t3,'(F10.2)') fmin(iw)
write(t4,'(F10.2)') sigfmin(iw)
write(t5,'(F10.2)') fmax(iw)
write(t6,'(F10.2)') sigfmax(iw)



b = fmax(iw)/fmin(iw)
sigb = b*sqrt( (sigfmin(iw)/fmin(iw))**2 + (sigfmax(iw)/fmax(iw))**2 )

write(*,*) wav(iw),fmin(iw), fmax(iw), &
b,sigb,fmin(iw)/fmax(iw)
read(*,*)

write(t7,'(F10.2)') b
write(t8,'(F10.2)') sigb

tat = '$'//trim(adjustl(t0))//'$'
tat = trim(adjustl(tat))//' & $'//trim(adjustl(t1))//'\pm'//trim(adjustl(t2))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t3))//'\pm'//trim(adjustl(t4))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t5))//'\pm'//trim(adjustl(t6))
tat = trim(adjustl(tat))//'$ & $'//trim(adjustl(t7))//'\pm'//trim(adjustl(t8))
tat = trim(adjustl(tat))//'$ //'

write(1,*) trim(adjustl(tat))
enddo
close(1)


!!! diagnostic plot !!!
if (diagnostic) then
do ilc = 1,NLC
!ier = pgopen('/Xserve')
!call pgsvp(0.1,0.9,0.5,0.9)
xmin = tgrid(1)
xmax = tgrid(Ntgrid)
ymin = minval(xlgrid)
ymax = maxval(xlgrid)
!call pgswin(xmin,xmax,ymin,ymax)
!call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!call pgline(Ntgrid,tgrid,xlgrid)

!call pgsvp(0.1,0.9,0.1,0.4)
xmin = tgrid(1)
xmax = tgrid(Ntgrid)
allocate(ymod(Ntgrid))

ymod(1:Ntgrid) = clam(ilc) + xlgrid(1:Ntgrid)*slam(ilc)
ymin = 0.0!minval(ymod)
ymax = maxval(ymod)

!call pgswin(xmin,xmax,ymin,ymax)
!call pgline(Ntgrid,tgrid,ymod)
!call pgbox('bcnst',0.0,0,'bcnst',0.0,0)

!call pgpt(idxlc(ilc),tnu(1:idxlc(ilc),ilc),fnu(1:idxlc(ilc),ilc),6)
!call pgerrb(6,idxlc(ilc),tnu(1:idxlc(ilc),ilc),fnu(1:idxlc(ilc),ilc),&
!sigfnu(1:idxlc(ilc),ilc), 1.0)

!call pgsci(2)
do it = 1,Ntgrid
if (inum_tot(it) .eq. 0) then
!call pgsci(4)
else if(inumgrid(it,ilc) .eq. 0) then
!call pgsci(3)
else
!call pgsci(2)
endif
!call pgpt(1,tgrid(it),fnugrid(it,ilc),6)
!call pgerrb(6,1,tgrid(it),fnugrid(it,ilc),&
!sigfnugrid(it,ilc), 1.0)
enddo

!call pgsci(1)
deallocate(ymod)
read(*,*)
!call pgend
enddo
endif








write(*,*) 'End of routine fvg_aug28.f90'

end subroutine













!! subroutine to calculate de-reddened disc spectrum
!! from variable spectrum and errors,


!!input wav(Nwav) wavlengths, wav0 ref wavelength, alam(Nwav), sigalam(Nwav)
! wavelength and uncertainties.

!! output f0, sigf0, ebmv, sigebmv (extinction, ref flux and uncertainties)

subroutine dered_jul15(Nwav, wav, wav0, alam, sigalam, f0, sigf0, ebmv, sigebmv)

integer, intent(in):: Nwav
real, intent(in):: wav(Nwav), wav0, alam(Nwav), sigalam(Nwav)
real, intent(out):: f0, sigf0, ebmv, sigebmv


real wav_wavref_log(Nwav), ext_lam(Nwav),alamlog(Nwav), sigalamlog(Nwav)
double precision alam_hes(2,2),alam_hes_inv(2,2)


!! convert flux and wavelength to log units
do ilc = 1, Nwav
wav_wavref_log(ilc) = alog10(wav(ilc)/wav0)
ext_lam(ilc) = 0.4*extmag_AGN(wav(ilc), 1.0)
call mystats(-1,alam(ilc),sigalam(ilc),alamlog(ilc),sigalamlog(ilc))
!write(*,*) alam(ilc), sigalam(ilc)
!write(*,*) alamlog(ilc), sigalamlog(ilc)
enddo







!! set up Hessian matrix
sum0 = 0.d0
sum1 = 0.d0
sum2 = 0.d0

do ilc = 1,Nwav
a = 1./sigalamlog(ilc)**2
b = ext_lam(ilc)
c = a*b
sum0 = sum0 + a
sum1 = sum1 + c
sum2 = sum2 + c*b

enddo
alam_hes(1,1) = sum0
alam_hes(1,2) = -1.*sum1
alam_hes(2,1) = -1.*sum1
alam_hes(2,2) = sum2
!! calculate inverse hessian matrix
call inverse(alam_hes,alam_hes_inv,2)


!! calculate par array 9array to time by inverse hessian to get parameters
sum0 = 0.d0
sum1 = 0.d0
do ilc = 1,Nwav
a = 1./sigalamlog(ilc)**2 * (alamlog(ilc) + 1./3*wav_wavref_log(ilc))
sum0 = sum0 + a
sum1 = sum1 - a*ext_lam(ilc)


enddo

f0log = alam_hes_inv(1,1)*sum0 + alam_hes_inv(1,2)*sum1
sigf0log = sqrt( alam_hes_inv(1,1) )

!! calculate error in f0 (from log f0 above)
call mystats(1,f0log,sigf0log,f0,sigf0)



alamwavarb = extmag_AGN(wav_arb, 1.0)



ebmv = alam_hes_inv(2,1)*sum0 + alam_hes_inv(2,2)*sum1
sigebmv = sqrt( alam_hes_inv(2,2) )

a = 0.d0
a =  a+ alam_hes_inv(2,1) + alam_hes_inv(2,2)
write(*,*) ebmv, sigebmv, f0log, sigf0log, f0, sigf0, sum1,sum0
write(*,*) 'his'
write(*,*) alam_hes(1,1),alam_hes(1,2)
write(*,*) alam_hes(2,1),alam_hes(2,2)

write(*,*) 'his_inv'
write(*,*) alam_hes_inv(1,1),alam_hes_inv(1,2)
write(*,*) alam_hes_inv(2,1),alam_hes_inv(2,2)



end subroutine



!! subroutine to calculate de-reddened disc spectrum
!! from variable spectrum and errors,


!!input wav(Nwav) wavlengths, wav0 ref wavelength, alam(Nwav), sigalam(Nwav)
! wavelength and uncertainties.
! i_whichlaw : Which reddening law to use?
!			 . 1 = gaskel agn 2004 law
!			 . 2 = smc gordon 03 law
!			 . 3 = lmc gordon 03 law
!			 . 4 = mw law



!! output f0, sigf0, ebmv, sigebmv (extinction, ref flux and uncertainties)

subroutine dered_aug28(Nwav, wav, wav0, alam, sigalam, f0, sigf0, ebmv, sigebmv,i_whichlaw)

integer, intent(in):: Nwav
real, intent(in):: wav(Nwav), wav0, alam(Nwav), sigalam(Nwav)
real, intent(out):: f0, sigf0, ebmv, sigebmv


real wav_wavref_log(Nwav), ext_lam(Nwav),alamlog(Nwav), sigalamlog(Nwav)
double precision alam_hes(2,2),alam_hes_inv(2,2)


!! convert flux and wavelength to log units
do ilc = 1, Nwav
wav_wavref_log(ilc) = alog10(wav(ilc)/wav0)

if (i_whichlaw .eq. 1) then
ext_lam(ilc) = 0.4*extmag_AGN(wav(ilc), 1.0)
else if (i_whichlaw .eq. 2) then
ext_lam(ilc) = 0.4*extmag_smc(wav(ilc), 1.0)
else if (i_whichlaw .eq. 3) then
ext_lam(ilc) = 0.4*extmag_lmc(wav(ilc), 1.0)
else if (i_whichlaw .eq. 4) then
ext_lam(ilc) = 0.4*extmag_mw(wav(ilc), 1.0)
endif


call mystats(-1,alam(ilc),sigalam(ilc),alamlog(ilc),sigalamlog(ilc))
!write(*,*) alam(ilc), sigalam(ilc)
!write(*,*) alamlog(ilc), sigalamlog(ilc)
enddo







!! set up Hessian matrix
sum0 = 0.d0
sum1 = 0.d0
sum2 = 0.d0

do ilc = 1,Nwav
a = 1./sigalamlog(ilc)**2
b = ext_lam(ilc)
c = a*b
sum0 = sum0 + a
sum1 = sum1 + c
sum2 = sum2 + c*b

enddo
alam_hes(1,1) = sum0
alam_hes(1,2) = -1.*sum1
alam_hes(2,1) = -1.*sum1
alam_hes(2,2) = sum2
!! calculate inverse hessian matrix
call inverse(alam_hes,alam_hes_inv,2)


!! calculate par array 9array to time by inverse hessian to get parameters
sum0 = 0.d0
sum1 = 0.d0
do ilc = 1,Nwav
a = 1./sigalamlog(ilc)**2 * (alamlog(ilc) + 1./3*wav_wavref_log(ilc))
sum0 = sum0 + a
sum1 = sum1 - a*ext_lam(ilc)


enddo

f0log = alam_hes_inv(1,1)*sum0 + alam_hes_inv(1,2)*sum1
sigf0log = sqrt( alam_hes_inv(1,1) )

!! calculate error in f0 (from log f0 above)
call mystats(1,f0log,sigf0log,f0,sigf0)





if (i_whichlaw .eq. 1) then
alamwavarb = 0.4*extmag_AGN(wav_arb, 1.0)
else if (i_whichlaw .eq. 2) then
alamwavarb = 0.4*extmag_smc(wav_arb, 1.0)
else if (i_whichlaw .eq. 3) then
alamwavarb = 0.4*extmag_lmc(wav_arb, 1.0)
else if (i_whichlaw .eq. 4) then
alamwavarb = 0.4*extmag(wav_arb, 1.0)
endif





ebmv = alam_hes_inv(2,1)*sum0 + alam_hes_inv(2,2)*sum1
sigebmv = sqrt( alam_hes_inv(2,2) )

a = 0.d0
a =  a+ alam_hes_inv(2,1) + alam_hes_inv(2,2)
!write(*,*) ebmv, sigebmv, f0log, sigf0log, f0, sigf0, sum1,sum0
!write(*,*) 'his'
!write(*,*) alam_hes(1,1),alam_hes(1,2)
!write(*,*) alam_hes(2,1),alam_hes(2,2)
!
!write(*,*) 'his_inv'
!write(*,*) alam_hes_inv(1,1),alam_hes_inv(1,2)
!write(*,*) alam_hes_inv(2,1),alam_hes_inv(2,2)



end subroutine



subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer n
double precision a(n,n), c(n,n),aold(N,N)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

aold(:,:)=a(:,:)
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
do i=k+1,n
coeff=a(i,k)/a(k,k)
L(i,k) = coeff
do j=k+1,n
a(i,j) = a(i,j)-coeff*a(k,j)
end do
end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
do i=1,j
U(i,j) = a(i,j)
end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
b(k)=1.0
d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
do i=2,n
d(i)=b(i)
do j=1,i-1
d(i) = d(i) - L(i,j)*d(j)
end do
end do
! Step 3b: Solve Ux=d using the back substitution
x(n)=d(n)/U(n,n)
do i = n-1,1,-1
x(i) = d(i)
do j=n,i+1,-1
x(i)=x(i)-U(i,j)*x(j)
end do
x(i) = x(i)/u(i,i)
end do
! Step 3c: fill the solutions x(n) into column k of C
do i=1,n
c(i,k) = x(i)
end do
b(k)=0.0
end do


a(:,:)=aold(:,:)
end subroutine inverse



!! convert mean parameter (aveold) and uncertainty sdold and mean both from (1) and to (-1) log par and uncertainty,
!and from (2) and to (-2) cosine par with uncertainty

subroutine mystats(iversion,aveold,sdold,avenew,sdnew)


integer,intent(in):: iversion
real,intent(in):: aveold, sdold
real,intent(out)::sdnew, avenew


rln10 = alog(10.)
deg2rad = 3.1415926536/180.

!! if going from log to real
if (iversion == 1) then
a = 10**aveold
avenew = a*(1. + 0.5*sdold*sdold*rln10*rln10)
sdnew  = rln10*a * sdold



!! real to log
else if (iversion == -1) then

a = alog10(aveold)
avenew = a - 0.5*sdold*sdold/(rln10*aveold*aveold)
sdnew  = 1./(rln10*aveold) * sdold


!write(*,*) 'mystats', aveold, sdold, avenew,sdnew, a, 0.5*sdold*sdold, (rln10*aveold*aveold),&
!0.5*sdold*sdold/(rln10*aveold*aveold), a - 0.5*sdold*sdold/(rln10*aveold*aveold)

!read(*,*)




!! cos to real
else if (iversion ==2) then
write(*,*) 'mystats.f90: this doesnt work yet!'
stop






!! real to cos
else if (iversion == -2) then

a = cos(aveold*deg2rad)

avenew = a*(1. - 0.5*sdold*sdold)
sdnew  = sin(aveold*deg2rad)*sdold





!! go from fnu to abmag
else if(iversion .eq. -3) then


sdnew = abs(2.5/rln10/aveold * sdold)

avenew = -2.5*alog10(aveold) + 1.25/(aveold*aveold*rln10) *sdold*sdold


!! go from abmag to fnu
else if(iversion .eq. 3) then

a = 10**(-2./5*(aveold))

sdnew = abs(0.4 * rln10 * a * sdold)

avenew = a* (1. + 2./25 * sdold*sdold * rln10*rln10)
!write(*,*) 'sdsds', sdnew, a, avenew,sdold,rln10,aveold

else


write(*,*) 'mystats.f90: iversion',iversion,' is not an option!!'
stop
endif


end subroutine



!!! fortran routine to pick a colour that represents the wavelength of the light curve
!! INPS wav(NLC)
!! OP   icwav(NLC)

subroutine sciwav(wav,icwav,NLC)

integer NLC, icwav(NLC)
real wav(NLC)

do ilc = 1,NLC
wavnow = wav(ilc)

if (wavnow .lt. 1500) then
icwav(ilc) = 15
cycle

else if (wavnow .lt. 2000) then
icwav(ilc) = 14
cycle

else if (wavnow .lt. 2500) then
icwav(ilc) = 16
cycle

else if (wavnow .lt. 3000) then
icwav(ilc) = 1
cycle

else if (wavnow .lt. 3500) then
icwav(ilc) = 12
cycle

else if (wavnow .lt. 3500) then
icwav(ilc) = 12
cycle

else if (wavnow .lt. 4500) then
icwav(ilc) = 6
cycle

else if (wavnow .lt. 5200) then
icwav(ilc) = 4
cycle

else if (wavnow .lt. 5500) then
icwav(ilc) = 5
cycle

else if (wavnow .lt. 6000) then
icwav(ilc) = 3
cycle

else if (wavnow .lt. 6500) then
icwav(ilc) = 8
cycle

else
icwav(ilc) = 2
cycle
endif

enddo


end subroutine



!! program to evaluate the uncertainties and parameter values of a model by minimising chi squared
!! tested 3rd august with ax**2+bx +c and works (code below)
subroutine hesfit(x,y,sig,p,N,NP,pout,sigpout)

!inputs x(N),y(N),sig(N) x y and uncertainty values
! p(N,NP) the pattern to be fitted i.e fitting y=ax^2 + bsin(x) p(1:N,1) would be x**2 and p(1:N,2) would be an array of sin(x)

! outputs: pout(NP), sigpout(NP parameters and uncertainties


real x(N),y(N),sig(N),p(N,NP),pout(NP),sigpout(NP),c(NP),sig2(N)
double precision hes(NP,NP),cov(NP,NP)

do it=1,N
a=sig(it)
sig2(it)=a*a
enddo



!!! set up the Hessian matrix

do ip2=1,NP

do ip1=1,NP
sum=0.d0
do it=1,N
sum=sum+p(it,ip1)*p(it,ip2)/sig2(it)
enddo

hes(ip1,ip2)=sum
!write(*,*) 'hes,cov',ip1,ip2,hes(ip1,ip2)
enddo!! end ip1
enddo!! end ip2

!!!find the inverse of the hessian matrix (the covariance matrix)
call inverse(hes,cov,np)

!stop
!!! set up the c matrix in:  a * hes = c (where a is a(NP), hes is the hessian matrix and c = sum (y_i p_i / sig^2_i) see ada lecture 9)

do ip=1,NP
sum=0.d0
do it=1,N
sum=sum+y(it)*p(it,ip)/sig2(it)
!write(*,*) it,sum
enddo
c(ip)=sum
enddo
!stop
!! find the solution
do ip0=1,NP
sum=0.d0
do ip1=1,NP
sum=sum+cov(ip1,ip0)*c(ip1)
enddo
pout(ip0)=sum
sigpout(ip0)=sqrt(cov(ip0,ip0)) !! and uncertainties... parameter variance are cov(1,1), cov(2,2) etc cov(n,n) if orthogonal, have no covariances
enddo
!!

end subroutine







































!! program to evaluate the uncertainties and parameter values of a model by minimising chi squared
!! tested 3rd august with ax**2+bx +c and works (code below)
subroutine hesfit2(x,y,sig,p,N,NP,pout,sigpout,val_ignore)

!inputs x(N),y(N),sig(N) x y and uncertainty values
! p(N,NP) the pattern to be fitted i.e fitting y=ax^2 + bsin(x) p(1:N,1) would be x**2 and p(1:N,2) would be an array of sin(x)

! outputs: pout(NP), sigpout(NP parameters and uncertainties


real x(N),y(N),sig(N),p(N,NP),pout(NP),sigpout(NP),c(NP),sig2(N)
double precision hes(NP,NP),cov(NP,NP)

do it=1,N
a=sig(it)
sig2(it)=a*a
enddo



!!! set up the Hessian matrix

do ip2=1,NP

do ip1=1,NP
sum=0.d0
do it=1,N
if (y(it) .eq. val_ignore) cycle !! skip values specified as being ignored in routine call
sum=sum+p(it,ip1)*p(it,ip2)/sig2(it)
enddo

hes(ip1,ip2)=sum
!write(*,*) 'hes,cov',ip1,ip2,hes(ip1,ip2)
enddo!! end ip1
enddo!! end ip2

!!!find the inverse of the hessian matrix (the covariance matrix)
call inverse(hes,cov,np)

!stop
!!! set up the c matrix in:  a * hes = c (where a is a(NP), hes is the hessian matrix and c = sum (y_i p_i / sig^2_i) see ada lecture 9)

do ip=1,NP
sum=0.d0
do it=1,N
if (y(it) .eq. val_ignore) cycle !! skip values specified as being ignored in routine call
sum=sum+y(it)*p(it,ip)/sig2(it)
!write(*,*) it,sum
enddo
c(ip)=sum
enddo
!stop
!! find the solution
do ip0=1,NP
sum=0.d0
do ip1=1,NP
sum=sum+cov(ip1,ip0)*c(ip1)
enddo
pout(ip0)=sum
sigpout(ip0)=sqrt(cov(ip0,ip0)) !! and uncertainties... parameter variance are cov(1,1), cov(2,2) etc cov(n,n) if orthogonal, have no covariances
enddo
!!

end subroutine























!!!! test this using y=ax + b
!real y(10),x(10),sig(10),p(10,3),pout(3),sigpout(3),ymod(10)
!
!n=10
!np=3
!a=1.2
!b=2.8
!c=5.0
!
!do i=1,n
!x(i)=i
!y(i)=a*x(i)+b*x(i)**2+c
!sig(i)=0.2
!enddo
!
!!! set up patterns
!do it=1,n
!p(it,1)=x(it)
!p(it,2)=x(it)**2
!p(it,3)=1
!enddo
!!!
!
!call hesfit(x,y,sig,p,N,NP,pout,sigpout)
!
!
!do it=1,N
!ymod(it)=pout(1)*x(i)+pout(2)*x(i)**2 + pout(3)
!enddo
!
!
!write(*,*) 'model output', p, sigp
!write(*,*) 'compare'
!write(*,*) x
!write(*,*) y
!write(*,*) ymod
!
!end program






!!!!!!!
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
FUNCTION FindDet(matrix, n)
IMPLICIT NONE
REAL, DIMENSION(n,n) :: matrix
INTEGER, INTENT(IN) :: n
REAL :: m, temp,finddet
INTEGER :: i, j, k, l
LOGICAL :: DetExists = .TRUE.
l = 1
!Convert to upper triangular form
DO k = 1, n-1
IF (matrix(k,k) == 0) THEN
DetExists = .FALSE.
DO i = k+1, n
IF (matrix(i,k) /= 0) THEN
DO j = 1, n
temp = matrix(i,j)
matrix(i,j)= matrix(k,j)
matrix(k,j) = temp
END DO
DetExists = .TRUE.
l=-l
EXIT
ENDIF
END DO
IF (DetExists .EQV. .FALSE.) THEN
FindDet = 0
return
END IF
ENDIF
DO j = k+1, n
m = matrix(j,k)/matrix(k,k)
DO i = k+1, n
matrix(j,i) = matrix(j,i) - m*matrix(k,i)
END DO
END DO
END DO

!Calculate determinant by finding product of diagonal elements
FindDet = l
DO i = 1, n
FindDet = FindDet * matrix(i,i)
END DO

END FUNCTION FindDet










function extmag_lmc(wa,ebmv)

!
!    Parameters
!    ----------
!    wa : real
!      One or more wavelengths in Angstroms.
!    ebmv: assumed ebmv (usually use 1)
!    Returns
!    -------
!    Ext : ExtinctionCurve instance
!      A(lambda) at each wavelength and Rv for the
!      LMC Average sample.
!
!    References
!    ----------
!    http://adsabs.harvard.edu/abs/2003ApJ...594..279G
!
!    Notes
!    -----
!    At far IR wavelengths, a MW extinction curve is assumed.


real,intent(in):: wa, ebmv
!real, intent(out):: extmag_lmc

c1 = -0.890    !sig = 0.142
c2 = 0.998     !sig = 0.027
c3 = 2.719     !sig = 0.137
c4 = 0.400     !sig = 0.036
x0 = 4.579     !sig = 0.007
gamma =0.934   !sig = 0.016


wa_micron = wa/1.e4
x = 1/wa_micron



a = c1 + c2*x + c3*smclmc_d(x,gamma,x0) + c4*smclmc_fx(x)
extmag_lmc = ebmv*a


end function






function extmag_smc(wa,ebmv)

!
!    Parameters
!    ----------
!    wa : real
!      One or more wavelengths in Angstroms.
!    ebmv: assumed ebmv (usually use 1)
!    Returns
!    -------
!    Ext : ExtinctionCurve instance
!      A(lambda) at each wavelength and Rv for the
!      LMC Average sample.
!
!    References
!    ----------
!    http://adsabs.harvard.edu/abs/2003ApJ...594..279G
!
!    Notes
!    -----
!    At far IR wavelengths, a MW extinction curve is assumed.

real,intent(in):: wa, ebmv
!real, intent(out):: extmag_smc

c1 = -4.959   !sig = 0.197
c2 = 2.264    !sig = 0.04
c3 = 0.389    !sig = 0.11
c4 = 0.461    !sig = 0.079
x0 = 4.6      !sig = 0
gamma =1.0    !sig = 0


wa_micron = wa/1.e4
x = 1/wa_micron

a = c1 + c2*x + c3*smclmc_d(x,gamma,x0) + c4*smclmc_fx(x)
extmag_smc= ebmv*a


end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! 2 functions that parameterise the extinction laws for smc and lmc

! f(x) =  0.5392*aa + 0.05644*aa*a  (where a = x - 5.9)
function smclmc_fx(x)

real, intent(in):: x
!real, intent(out):: smclmc_fx

a = x - 5.9
aa = a * a

if (x > 5.9) then
smclmc_fx = 0.5392*aa + 0.05644*aa*a
else
smclmc_fx = 0
endif

end function



! D(x,gamma,x0) = x^2 / ((x^2 - x0^2)^2 + x^2gamma^2)
function smclmc_d(x,gamma,x0)

real, intent(in):: x, gamma, x0
!real, intent(out):: smclmc_d

a = x*x - x0*x0
b = x*gamma

smclmc_d = x*x / (a*a + b*b)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!! program to test the extinction curves above
! gflags extmag_smclmc.f90 ../extmag.for ../lookup.for ../extmag_agn.for
!program testsmclmc
!
!
!real,allocatable:: wavs(:), alamlmc(:), alamsmc(:), alammw(:), alamagn(:)
!
!
!ebmv = 1.0
!wavlo = 0.0
!wavhi = 0.001
!dwav  = 0.00001
!
!nwav = (wavhi - wavlo)/dwav + 1
!
!allocate(wavs(Nwav), alamlmc(Nwav), alamagn(Nwav), alamsmc(Nwav),alammw(Nwav))
!
!do i = 1,Nwav
!wavnow = wavlo + (i-1)*dwav
!wavs(i) = wavnow
!alamlmc(i) = extmag_lmc(1./wavnow,ebmv)
!alamsmc(i) = extmag_smc(1./wavnow,ebmv)
!alammw(i)  = extmag(1./wavnow,ebmv)
!alamagn(i)  = extmag_agn(1./wavnow,ebmv)
!!write(*,*) wavs(i), alamlmc(i), alamsmc(i)
!enddo
!
!
!xmin = wavlo
!xmax = wavhi
!ymin = min(minval(alamlmc),minval(alamsmc),minval(alammw),minval(alamagn))
!ymax = max(maxval(alamlmc),maxval(alamsmc),maxval(alammw),maxval(alamagn))
!
!!ier = pgopen('test_extmag_smclmc.ps/CPS')
!!call pgsvp(0.2,0.9,0.2,0.9)
!
!
!!call pgsch(2.0)
!!call pgslw(3)
!
!!call pgswin(xmin,xmax,ymin,ymax*1.1)
!!call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
!
!
!
!
!!call pgline(Nwav,wavs,alammw)
!
!!call pgmtxt('B',2.5, 0.5, 0.5,'1/[Wavelength \A]')
!!call pgmtxt('L',2.5, 0.5, 0.5,'A(\gl)')
!
!!call pgsci(4)
!!call pgline(Nwav,wavs,alamlmc)
!
!!call pgsci(2)
!!call pgline(Nwav,wavs,alamsmc)
!
!!call pgsci(3)
!!call pgline(Nwav,wavs,alamagn)
!
!!call pgend
!
!end program

!
!
!
!
!
!
!
!



!! return number of unique wavelengths
subroutine wav_unique(wav,Ntel,Nwav)
integer Nwav, Ntel
real wav(Ntel)

!!! determine nuber of unique wavlengths
ilcold = 1
idxu = 1
wavold = wav(ilcold)
do ilc = 1,Ntel
wavnow = wav(ilc)
if (wavnow .ne. wavold) then
ilcold = ilc
idxu = idxu + 1
wavold = wavnow
endif
enddo
Nwav = idxu

end subroutine




!!fvg_prep


subroutine fvg_prep_multitel(Ntel, Nwav, wav, idxhi,Ndatmax, t, x, sig, t2, x2, sig2,&
wav2,idxlc2,ilclo2)

integer Ntel, Nwav, idxhi(Ntel),Ndatmax, ilchi2(Nwav),idxlc2(Nwav), ilclo2(Nwav)
real wav(Ntel), t(Ndatmax,Ntel), x(Ndatmax,Ntel), sig(Ndatmax,Ntel), t2(Ndatmax,Nwav), x2(Ndatmax,Nwav),&
sig2(Ndatmax,Nwav), wav2(Nwav) , Ndat(Nwav)
real,allocatable,dimension(:):: xtemp,ttemp,sigtemp,ttemp_s,xtemp_s,sigtemp_s
integer,allocatable,dimension(:):: idxkey




ilclo = 1
ilclo2(1) = 1
idx = 1
wav2(1) = wav(1)
do ilc = 1,Ntel
wavold = wav(ilclo)
if (wav(ilc) == wavold) cycle
ilclo = ilc

idx = idx + 1
ilclo2(idx) = ilclo
wav2(idx) = wav(ilclo)
if (idx .ge. Nwav) exit
enddo


ilclonow = 1


!do ilc = 1,Nwav
!write(*,*) ilc, ilclo2(ilc),idxhi(ilc)
!enddo


do ilc = 1,Nwav
!idx = ilc+1

if (ilc .gt. 1) ilclonow = ilclo2(ilc)
if (ilc .lt. Nwav) then
ilchinow = ilclo2(ilc + 1)-1
else
ilchinow = Ntel
endif


ndatnow = 0
if (ilc .eq. 1) then
ndatnow = idxhi(ilchinow)
else
!write(*,*) 'why is thi shappening', ilclonow,ilchinow,ndatnow,idxhi(ilclonow)
do i = ilclonow,ilchinow
ndatnow = ndatnow + idxhi(i)
enddo
endif

!write(*,*) ndatnow
idxlc2(ilc) = ndatnow
allocate(ttemp(Ndatnow),xtemp(Ndatnow),sigtemp(Ndatnow),idxkey(Ndatnow),&
ttemp_s(Ndatnow),xtemp_s(Ndatnow),sigtemp_s(Ndatnow))
idx = 1
do ilc2 = ilclonow, ilchinow
idxhinow = idxhi(ilc2)
do ic = 1,idxhinow
ttemp(idx)   = t(ic,ilc2)
xtemp(idx)   = x(ic,ilc2)
sigtemp(idx) = sig(ic,ilc2)
!write(*,*) 'test2', ilc,ilc2,idx,idxhinow,ic,t(ic,ilc2),ttemp(idx)
idx = idx + 1

enddo
enddo

!write(*,*) idx,ndatnow,'blalblablbabalb',ilclonow,ilchinow
!read(*,*)
!write(*,*) ndatnow, idx,'test 2',ilc
!! order the light curves using sorrt1d subroutine
call sort1d(ttemp,Ndatnow,ttemp_s,idxkey)
do i = 1,Ndatnow

inow = idxkey(i)
ttemp_s(i) = ttemp(inow)
!write(*,*) 'test2', ttemp(i), ttemp_s(i),ilc
xtemp_s(i) = xtemp(inow)
sigtemp_s(i) = sigtemp(inow)
enddo
!read(*,*)
t2(1:Ndatnow,ilc)    = ttemp_s(1:Ndatnow)
x2(1:Ndatnow,ilc)    = xtemp_s(1:Ndatnow)
sig2(1:Ndatnow,ilc)  = sigtemp_s(1:Ndatnow)
!write(*,*) ilc, xtemp(1), sigtemp(1),'checking prep...',xtemp_s(Ndatnow),sigtemp_s(1)
!read(*,*)
deallocate(ttemp,xtemp,sigtemp,idxkey,ttemp_s,xtemp_s,sigtemp_s)

enddo




end subroutine



!! subroutine to create a top hat transfer function


subroutine tftophat(Ntau, tau, psi, taucent, taufwhm)

integer, intent(in):: Ntau
real, intent(in):: tau(Ntau), taucent, taufwhm
real, intent(out):: psi(Ntau)

tauhwhm = taufwhm/2
taulo = taucent - tauhwhm
tauhi = taucent + tauhwhm
dtaugrid = (tau(Ntau) - tau(1))/(Ntau - 1)



if (dtaugrid .ge. (tauhi - taulo)) then !problem, we must always have at least 1 point being '1'

tausepmin = 1.e6
do itau = 1,Ntau
tsminnow = abs(tau(itau) - taucent)
if (tsminnow .lt. tausepmin) then
idxmin = itau
tausepmin = tsminnow
endif
enddo
do itau =1,Ntau
if (itau .ne. idxmin) then
psi(itau) = 0
else
psi(itau) = 1
endif
enddo

else !otherwise no problem and carry on as normal

do itau = 1,Ntau
taunow = tau(itau)
if (taunow .lt. taulo) then
psinow = 0
else if (taunow .gt. tauhi) then
psinow = 0
else
psinow = 1.
endif
psi(itau) = psinow
enddo

endif

!top = 0
!bot = 0
!do itau = 1,Ntau
!top = top +  tau(itau)*psi(itau)
!bot = bot +  psi(itau)
!enddo
!write(*,*) 'i am making a new tfbox, cent=', taucent, top/bot
!read(*,*)

end subroutine



!******************************************
subroutine avgrms( n, x, bar, rms )
!* unweighted mean and rms of n data points x(i)
!* Input:
!*	n	i4 number of data
!*	x(n)	r4 data
!* output:
!*	bar	r4 mean value
!*	rms	r4 root-mean-square
!* Running update method (B.P.Welford 1962)
!* Seems to be more robust to round-off errors.
!* 2011 Feb Keith Horne @ St.Andrews
real*4 x(*)
real*8 sum, sum2, add
if( n .le. 1 ) then
bar = x(1)
rms = 0.
return
end if
sum1 = x(1)
sum2 = 0.d0
do i=2,n
add = dble( x(i) ) - sum1
sum1 = sum1 + add / i
sum2 = sum2 + add * ( dble( x(i) ) - sum1 )
end do
bar = sum1
var = sum2 / max( 1, n - 1 )
rms = sqrt( var )
return
end

!******************************************
subroutine avgrms1( n, x, bar, rms )
!* unweighted mean and rms of n data points x(i)
!* Input:
!*	n	i4 number of data
!*	x(n)	r4 data
!* output:
!*	bar	r4 mean value
!*	rms	r4 root-mean-square
!* 1999 Sep Keith Horne @ St.Andrews
!* 2011 Feb KDH @ StA - x(i)**2 -> x(i)*x(i)
real*4 x(*)
real*8 sum, sum2
bar = x(1)
rms = 0.
if( n.le.1 ) return
sum = 0.d0
sum2 = 0.d0
do i=1,n
sum = sum + x(i)
sum2 = sum2 + x(i) * x(i)
end do
bar = sum / n
var = ( sum2 - sum * sum / n ) / max( 1, n-1 )
rms = sqrt( max( 0., var ) )
return
end




!! Input: t,1d array of data size N. op output. the sorted, output array
subroutine sort1d(t,N,op,key)
integer N,i,j
real t(N), op(N),temp,op2d(N,2),tempkey
integer key(N)

op2d(1:N,1)=t(1:N)

do i=1,N
op2d(i,2)=i
enddo


do i =1,N-1
do j =i+1,N

if (op2d(i,1) .ge. op2d(j,1)) then
temp = op2d(j,1)
tempkey=op2d(j,2)
op2d(j,1)=op2d(i,1)
op2d(j,2)=op2d(i,2)
op2d(i,1)=temp
op2d(i,2)=tempkey

endif
enddo
enddo

op(1:N)=op2d(1:N,1)
key(1:N)=op2d(1:N,2)
!      write(*,*) op2d(1:N,1)
!      write(*,*) ''
!      write(*,*) key

return
end subroutine



!*************************
subroutine append( str1, str2 )
!* append two strings with a space in between
!* 2003 Jun Keith Horne @ St.Andrews
character*(*) str1, str2
l1 = len1( str1 )
l2 = len1( str2 )
if( l1 .le. 0 ) then
str1 = str2( : l2 )
else
str1 = str1( : l1 ) // ' ' // str2( : l2 )
end if
return
end
!****************************************
subroutine append_i( str, i )
!* append integer to string
!* 2003 Jun Keith Horne @ St.Andrews
character*(*) str
character*40 number
write( number, * ) i
call word1( number, l1, l2 )
call append( str, number( l1 : l2 ) )
return
end
!****************************************
subroutine append_f( str, x, fmt )
!* append formatted real to string
!* 2003 Jun Keith Horne @ St.Andrews
!* 2008 Aug KDH @ St.A - add '*' format
character*(*) str, fmt
character*40 number
if( fmt(1:1) .eq. '*' ) then
write( number, * ) x
else
write( number, fmt ) x
end if
call word1( number, l1, l2 )
call append( str, number( l1 : l2 ) )
return
end
!****************************************
subroutine append_d( str, x, fmt )
!* append formatted real*8 to string
!* 2004 Feb Keith Horne @ St.Andrews
real*8 x
character*(*) str, fmt
character*40 number
write( number, fmt ) x
call word1( number, l1, l2 )
call append( str, number( l1 : l2 ) )
return
end

!****************************************
subroutine append_n( title, val, idig )
!* append val to title with idig significant digits
!* 2008 Jun Keith Horne @ St.Andrews
!* 2008 Jul KDH @ SAAO - idig
character*(*) title
character*20 fmt
i = max( 2, idig )
a= abs( val ) / 10**(i)
if( a .eq. 0.0 ) then
call append_f( title, val, '(f10.1)' )
else if( a .lt. 9.9e-8) then
fmt = '(1pe'
call append_i( fmt, i + 7 )
call append( fmt, '.' )
call append_i( fmt, i - 1 )
call append( fmt, ')' )
call nospace( fmt )
call append_f( title, val, fmt )
else if( a .lt. 9.9e-7) then
call append_f( title, val, '(f15.6)' )
else if( a .lt. 9.9e-6) then
call append_f( title, val, '(f15.5)' )
else if( a .lt. 9.9e-5 ) then
call append_f( title, val, '(f15.4)' )
else if( a .lt. 9.9e-4 ) then
call append_f( title, val, '(f15.3)' )
else if( a .lt. 9.9e-3 ) then
call append_f( title, val, '(f15.2)' )
else if( a .lt. 0.099 ) then
call append_f( title, val, '(f15.1)' )
else
call append_i( title, nint( val ) )
end if
return
end

!********************************************
subroutine append_datsig( title, dat, sig, kfig )
!* Input
!*       title   c* string
!*       dat     r4 data
!*       sig     r4 sigma
!*       kfig    i4 number of significant figures on sigma
!* Output
!*       title   c* string // data(sigma)
!* 2004 Feb Keith Horne @ St-Andrews
!* 2013 May KDH @ StA trap sig not positive
character*(*) title
real*8 hjd
if( sig .le. 0. ) then
call append_n( title, dat, max( 2, kfig ) )
else
hjd = dat
call append_hjdsig( title, hjd, sig, kfig )
end if
return
end

!********************************************
subroutine append_hjdsig( title, hjd, sig, kfig )
!* Input
!*       title   c* string
!*       hjd     r8 data
!*       sig     r4 sigma
!*       kfig    i4 number of significant figures on sigma
!* Output
!*       title   c* string // data(sigma)
!* 2004 Feb Keith Horne @ St-Andrews
!* 2012 Dec KDH @ StA omit E0, skip if sig not positive
!* 2013 May KDH @ StA trap sig not positive
character*(*) title
real*8 hjd
character*50 fmt

if( sig .le. 0. ) then
call append_n( title, real(hjd), max( 2, kfig ) )
return
end if

l1 = len1( title )
mpow = nint( log10( abs( hjd * 2. ) ) )
ipow = nint( alog10( sig * 2. ) )
kpow = ipow - kfig
div = 10. ** kpow

if( kpow .lt. 0 ) then
fmt = '(f20.'
call append_i( fmt, -kpow )
call append( fmt, ')' )
call nospace( fmt )
call append_d( title, hjd, fmt )

if( sig .gt. 0. ) then
sigdiv = sig / div
call append( title, '(' )
if( sig .lt. 1. ) then
call append_i( title, nint( sigdiv ) )
else
call append_f( title, sig, fmt )
end if
call append( title, ')' )
end if

else if( kpow .gt. 3 ) then
lpow = mpow
if( iabs( ipow ) .gt. iabs( mpow ) ) lpow = ipow
div = 10. ** lpow
fmt = '(f20.'
call append_i( fmt, kfig )
call append( fmt, ')' )
call nospace( fmt )
call append_d( title, hjd / div, fmt )

if( sig .gt. 0. ) then
call append( title, '(' )
call append_f( title, sig / div, fmt )
call append( title, ')' )
if( lpow .ne. 0 ) then
call append( title, 'E' )
call append_i( title, lpow )
end if
end if

else
idiv = 10 ** kpow
call append_i( title, nint( hjd / idiv ) * idiv )
if( sig .gt. 0. ) then
call append( title, '(' )
call append_i( title, nint( sig / idiv ) * idiv )
call append( title, ')' )
end if

end if
call nospace( title(l1+2:) )
return
end

!****************************************
subroutine nospace( str )
!* remove spaces from string
!* 2003 Aug Keith Horne @ St.Andrews
character*(*) str
n = len1( str )
if( n .lt. 1 ) return
j = 0
do i=1,n
if( str( i : i ) .ne. ' ' ) then
j = j + 1
str( j : j ) = str( i : i )
end if
end do
str = str( : j )
return
end
!****************************************
subroutine notail( str, chr )
!* remove specified character from end of string
!* 2003 Dec Keith Horne @ St.Andrews
character*(*) str, chr
n = len1( str )
if( n .lt. 1 ) return
do while( n .gt. 0 .and. str(n:n) .eq. chr(1:1) )
n = n - 1
end do
str = str(:n)
return
end
!*------------------------------
subroutine noeqbl( title )
!* change '= ' to '='
!* 2006 Feb Keith Horne @ St.Andrews
character*(*) title
1	l = index( title, '= ' )
if( l .le. 0 ) return
title = title(:l) // title(l+2:)
goto 1
end
!*------------------------------
subroutine nogreek( title )
!* change pgplot greeks to text
!* 2007 Aug Keith Horne @ LaSilla
character*(*) title
parameter ( nfix = 15 )
character*3 old( nfix )
character*5 new( nfix )
data old/ &
'\ga',	'\gb',	'\gg',	'\gg',	'\ge',&
'\gh',	'\gy',	'\gm',	'\gn',	'\gp',&
'\gr',	'\gs',	'\gx',	'\gD',	'\gW'/

data new/ &
'alp',	'bet',	'gam',	'del',	'eps',&
'theta', 'eta',	'mu',	'nu',	'pi',&
'rho',	'sig',	'chi',	'Del',	'Omega'/


do k=1,nfix
i = 1
do while( i .gt. 0 )
i = index( title, old(k) )
if( i .gt. 0 ) then
n = len1( new(k) )
title = title(:i-1) // new(k)(:n) // title(i+3:)
end if
end do
end do

return
end
!*------------------------------
subroutine noupdn( title )
!* change pgplot sub and superscripts to text
!* 2007 Aug Keith Horne @ LaSilla
character*(*) title
m = 1
do while( m .gt. 0 )
i = index( title, '\u' )
j = index( title, '\d' )
if( i .gt. 0 .and. j .gt. i ) then
title = title(:i-1) // '^' // title(i+2:j-1) // title(j+2:)
else if( j .gt. 0 .and. i .gt. j ) then
title = title(:j-1) // '_' // title(j+2:i-1) // title(i+2:)
else
m = 0
end if
end do
return
end





function extmag_mw( wave, ebmv )
extmag_mw = extmag( wave, ebmv )
return
end

!C*EXTMAG ... interstellar extinction function from Seaton(1979) and Nandy(1975)
!C+
FUNCTION EXTMAG( WAVE, EBMV )
!*
!* interstellar extinction law
!*	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
!*	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
!*
!* Input:
!*	WAVE	R4 wavelength (Angstroms)
!*	EBMV	R4 extinction parameter E(B-V) (mags)
!*
!* Output:
!*	EXTMAG	= extinction (mags)
!C--
!* Mar 1986 Keith Horne @ STScI -- adapted from R.WADE routine SEATON
!* May 1989 Keith Horne @ STScI -- improve extrapolation
!* Jun 1991 KDH @ USM - minor changes
!*
!* COMMENTS FROM R.WADE SUBROUTINE:
!c seaton's paper in m.n.r.a.s. vol 187, page 75p (1979).
!c the formulae are based on an adopted value of R = 3.20.
!c
!c note that seaton's representation of of the interstellar reddening law
!c differs substantially from schild's representation (astron. j. 82, 339,
!c table ii, 1977) in the region of overlap.  schild also adopted r = 3.20.
!c
!c for wavelengths > 3704 angstroms, the program interpolates
!c linearly in 1/lambda c in seaton's table 3.
!c for wavelengths < 3704 angstroms, the program uses the formulae
!c from seaton's table 2.
!c the formulae match at the endpoints of their respective intervals.
!c there is a mismatch of 0.009 mag/ebmv at nu=2.7 (lambda=3704 angstroms).
!c seaton's tabulated value of 1.44 mags at 1/lambda = 1.1 may be in error;
!c 1.64 seems more consistent with his other values.
!c
!c wavelength range allowed is 0.1 to 1.0 microns.
!c calls to the subroutine for wavelengths outside this range
!c result in extrapolated values for the extinction being returned.
!*
!* sources:
!*	lambda < 1000		same as lambda = 1000.
!*	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
!*	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
!*	10000 < lambda		extrapolate linearly in 1/lam (can be improved)
!* 2001 Jul KDH @ St.And - PARAMETER()

PARAMETER ( NTABLE = 19 )
REAL*4 XTABLE(NTABLE), ETABLE(NTABLE)

!* tabulated inverse wavelengths (1/micron)
DATA  XTABLE/ 0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,&
1.6, 1.7, 1.8, 1.9, 2.0, 2.1,&
2.2, 2.3, 2.4, 2.5, 2.6, 2.7	/

!* tabulated extinctions at E(B-V)=1.
!c      DATA    ETABLE/ 0., 1.36, 1.44, 1.84, 2.04, 2.24, 2.44,
DATA  ETABLE/ 0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44,&
2.66, 2.88, 3.14, 3.36, 3.56, 3.77,&
3.96, 4.15, 4.26, 4.40, 4.52, 4.64/


EXTMAG = 0.
IF( WAVE.LE.0. ) RETURN
IF( EBMV.EQ.0. ) RETURN
X = 10000. / WAVE

!* infrared - extend optical results linearly to 0 at 1/lam=0
IF( X.LE.1.0 ) THEN
EXTMAG = ETABLE(2) * X * X

!* optical - interpolate linearly in magnitude vs 1/lam from Seaton's Table 3
ELSE IF ( X.LT.2.7 ) THEN
10	CALL LOOKUP( NTABLE, XTABLE, X, ILO, IHI, PART )
EXTMAG = ETABLE(ILO)*(1.-PART) + ETABLE(IHI)*PART

!* ultraviolet - use analytic formulae from Seaton's Table 2
ELSE IF( X.LT.3.65 ) THEN
DIFF = X - 4.6
EXTMAG = 1.56 + 1.048*X + 1.01/( DIFF*DIFF + 0.280)

ELSE IF( X.LT.7.14 ) THEN
DIFF = X - 4.6
EXTMAG = 2.29 + 0.848*X + 1.01/( DIFF*DIFF + 0.280)

ELSE IF( X.LE.10. ) THEN
EXTMAG = 16.17 + X*(-3.20 + 0.2975*X)

!* far-uv - generally should not occur, Lyman edge is at 912 A.
ELSE
X = MIN( X, 50. )
EXTMAG = 16.17 + X*(-3.20 + 0.2975*X)

END IF

EXTMAG = EBMV * EXTMAG
RETURN
END




function extmag_agn( angst, ebmv )
! AGN dust extinction ala Gaskell et al. 2004 ApJ 616,147
! input:
!       angst   r4 wavelength in angstroms
!       ebmv    r4 E(B-V) in magnitudes
! output:
!       extmag_agn      r4 AGN extinction in magnitudes
! 2012 Aug Keith Horne @ St And
x = 1.e4 / angst
! extrapolate
if( x .lt. 1.6 ) then
x0 = 1.6
y0 = x0 * ( x0 * ( x0 * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
dydx = x0 * ( x0 * 3*0.0296 - 2*0.377 ) + 1.5848
y = y0 + dydx * ( x - x0 )
! uv
else if( x .lt. 3.69 ) then
y = x * ( x * ( x * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
! optical
else if( x .lt. 8. ) then
y = 1.3468 + 0.0087 * x
! extrapolate
else
y = 1.3468 + 0.0087 * x
end if
rv = 5.15
av = rv * ebmv
extmag_agn = av*y!2.5 * alog10( av * y )
return
end




!************************************
integer function len1( string )
!* length up to last non-blank character
character*(*) string
do i = len( string ), 1, -1
if( string( i:i ) .ne. ' ' ) then
len1 = i
return
end if
end do
len1 = 0
return
end




!C*WORD1 -- isolate first word of a text string
!C+
SUBROUTINE WORD1( STRING, L1, L2 )
!*
!* Isolates the first word of a character string
!*
!* Input:
!*	STRING	= character string
!* Output:
!*	L1,L2	= Inclusive index limits of first word
!C--
!* Jan 1986 Keith Horne @ STScI
!* Jan 1989 KDH @ STScI - revise to allow BLANK or TAB separators
CHARACTER*(*) STRING
CHARACTER*1 TAB, BLANK, TEST
TAB = CHAR(9)
BLANK = ' '
LMAX = LEN1(STRING)
DO I1=1,LMAX
TEST = STRING(I1:I1)
IF( TEST.NE.BLANK .AND. TEST.NE.TAB ) THEN
L1 = I1
DO I2=L1+1,LMAX
TEST = STRING(I2:I2)
IF( TEST.EQ.BLANK .OR. TEST.EQ.TAB ) THEN
L2 = I2-1
RETURN
END IF
END DO
L2 = LMAX
RETURN
END IF
END DO
L1 = 0
L2 = 0
RETURN
END




!! DFT2 gives frequencies in cycs per day

!! NOTE all x values of input series in time space MUST start at zero if you are
! planning an inverse fourier trasnform. Otherwise there will be a constant time offset
! in the ift data to the input data. (Just subtract the lowest x value before taking dft
! and add it back on at the end. This error cost me most of a day on 10th dec 2014
! Seriously, I almost went berserk!

!.  also does not fit zero frequency component so when ift'ing a function,
! must add on mean level of original time series to get back original time series.


!! fourier transform of an input data set

subroutine  dft(t,x,N,ft,fa)
! IP
!!! t: input times
!!! x: input function of time
!!! N: data points

! OP
!!! ft fourier transform (this is complex)
!!! fa:absolute value of fourier transform


integer N
real x(N), t(N), dt, pi, w,fa(N)
complex ft(N)

pi=4.*atan2(1.,1.)
twopi=2*pi

dt=(t(N)-t(1))/(N-1)
do it = 1,N
ft(it)=(0,0)
do ik = 0,N-1

w=(twopi*ik*it)/N
ft(it)=ft(it)+x(ik+1)*cmplx(cos(w),sin(w))    !!the x(ik+1) was changed from x(ik) on 21nov 2013
enddo

fa(it)=abs(ft(it))
enddo

return
end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! fourier transform of an input data set

subroutine  dft2(t,x,N,ft,fa,freq)
! IP
!!! t: input times
!!! x: input function of time
!!! N: data points

! OP
!!! ft complex fourier transform
!!! fa:absolute value of fourier transform
!!! freq(ceiling(N/2)) array of frequencys corresponding to the first half of the fourier transform array freq(i)=i/(t(N)-t(1)) in cycs/time NOT angular units

integer N
real x(N), t(N), dt, pi, w,fa(N),freq(ceiling(0.5*N))
complex ft(N)

pi=4.*atan2(1.,1.)
twopi=2*pi
NW=ceiling(0.5*N)
ft(1:N) = 0

dt=(t(N)-t(1))/(N-1)
do ik = 1,N
ft(ik)=(0,0)
do it = 0,N-1
w=(twopi*ik*it)/N

ft(ik)=ft(ik)+x(it+1)*cmplx(cos(w),sin(w))    !!the x(ik+1) was changed from x(ik) on 21nov 2013
enddo
fa(ik)=abs(ft(ik))
enddo


tlen=t(N)-t(1)
do iw=1,NW
freq(iw)=1.*iw/tlen
enddo


return
end










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














!! fourier transform of an input data set

subroutine  dft_fast(t,x,N,ft)
! IP
!!! t: input times
!!! x: input function of time
!!! N: data points

! OP
!!! ft complex fourier transform
!!! fa:absolute value of fourier transform
!!! freq(ceiling(N/2)) array of frequencys corresponding to the first half of the fourier transform array freq(i)=i/(t(N)-t(1)) in cycs/time NOT angular units

integer N
real x(N), t(N), dt, pi, w,freq(ceiling(0.5*N))
complex ft(N)


pi=4.*atan2(1.,1.)
twopi=2*pi
NW=ceiling(0.5*N)
ft(1:N) = 0

dt=(t(N)-t(1))/(N-1)
do ik = 1,N
ft(ik)=(0,0)
do it = 0,N-1
w=(twopi*ik*it)/N

ft(ik)=ft(ik)+x(it+1)*cmplx(cos(w),sin(w))    !!the x(ik+1) was changed from x(ik) on 21nov 2013
enddo
enddo


return
end










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!! subroutine to calculate the inverse fourier transform from with an input fourier transform
!!! input t(NT) times desired output times to evaluate x(NT)
!         w(NFT) 2pi*freq REMEMBER THIS - this lost me another day! -bollocksy fucking fortran!
!         dft(NFT) fourier transform complex(NFT)


!!! output x(NT) time series evaluated at times t(1:NT)


subroutine ift2(t,w,dft,NFT,NT,x)

integer NFT,NT
real t(NT),x(NT),w(NFT),dftb(NFT,2)
complex dft(NT)




do i=1,NFT
dftb(i,1)=real(dft(i))
dftb(i,2)=aimag(dft(i))
!       write(*,*) 'real',dftb(i,1),'imag',dftb(i,2),dft(i)
enddo



do it=1,NT

time=t(it)
sumreal=0.d0
sumimag=0.d0
do iw=1,NFT
sumreal=sumreal+dftb(iw,1)*cos(w(iw)*time)
sumreal=sumreal+dftb(iw,2)*sin(w(iw)*time)   !!(a + ib)*(c + id) but d is -ve

sumimag=sumimag+dftb(iw,2)*cos(w(iw)*time)
sumimag=sumimag-dftb(iw,1)*sin(w(iw)*time)
enddo

opreal=sumreal/NFT
opimag=sumimag/NFT
opabs=sqrt(opreal*opreal+opimag*opimag)

x(it)=opreal

enddo

end subroutine













!! test 1tp 1d2 6th jan

!      program test_1tp1d2

!      real y1(10),x1(10),y2(6),x2(6)
!      integer choplo,chophi

!      do i = 1,10
!      x1(i)=i
!      y1(i)=i
!      enddo

!      do i=1,6
!      x2(i)=i+7

!      enddo

!      x2(6)=4.5
!      x2(1)=-2.0



!      call itp1d(y2,choplo,chophi,x1,x2,y1,10,6,0)

!      write(*,*) x1
!      write(*,*) y1
!      write(*,*) ''
!      write(*,*) x2
!      write(*,*) y2

!      end program


!!! subroutine to interpolate an entir array of y values given 1 array of y values,
!!! the corresponding x values and the y values to be interpolated


!!! 6th jan 2014 tested and works with above code!!

!!! inputs 1D array: x1(N1),x2(N2),y1(N1)
!!! integer N1 or N2 elements in each array
!!! output 1D array y(N2)

subroutine itp1d(y,choplo,chophi,x1,x2,y1,N1,N2,set0)

integer i,gtn,lt,gt,set0,choplo,chophi,N1,N2
real y(N2),x1(N1),x2(N2),y1(N1)

i=1
choplo=1
chophi=N2

x11 = x1(1)
x12 = x1(2)
x1N1 = x1(N1)
y1N1 = y1(N1)
y12 = y1(2)
do i=1,N2
!!!!! the first do loop deals with interpolation of x2 points less than x1(1)
if  (x2(i) .le. x11 .and. i .le. N2) then
if (set0 .eq.0) then
y(i) =0.0
else
y(i)=y1(1)+((y12-y1(1))/(x12-x11))*(x2(i)-x11)
endif
!	  i=i+1
choplo=choplo+1



!!! the next loop deals with interpolation of all x2 values between x1 values
elseif (x2(i) < x1N1 .and. i .le. N2) then
lt=gtn(x1,x2(i),1,N1)-1
gt=lt+1
y(i)=y1(lt)+((x2(i)-x1(lt))/(x1(gt)-x1(lt)))*(y1(gt)-y1(lt))


!      i=i+1



!!! the final loop deals with interpolation of x2 points beyond the final x1 point

else !(i .le. N2 .and. x2(i) .ge. x1(N1))


if (set0 .eq. 0) then
y(i)=0
else
y(i)=y1N1+((y1(N1)-y1(N1-1))/(x1N1-x1(N1-1)))*(x2(i)-x1N1)
endif

!      i=i+1
endif
enddo

chophi=chophi-1
return
end







!!!! David Starkey 2013 November
!! function to interpolate linearly between inputted x and y values
!input: x1, x2, y1, y2 The known xand y values either side of the value...
! !      x, to be interpolated
!! output       op the outputted interpolated y value


subroutine itp(x1,x2,x,y1,y2,op)
real:: x1,x2,x,y1,y2,op,p,big

!! which is larger x1 or x2?
if (x2 > x1) then
big=x2
small=x1
else
big=x1
small = x2
endif

if ( x > big .or. x < small) then
write(*,*) ' interpolation error'
endif
p= x - x1
p= p/(x2 - x1)

if (x .lt. small .or. x .gt. big .or. abs(p) > 1) then
write(*,*) ' interpolation error',x1,x2,x,p
stop
endif
op = y1 + p*(y2 - y1)
return
end





!!! return 1st, 2nd... nth, element of an array less than some value
!! yours truly (ps today was the day i fainted in keith's office)
! input x: 1d array, lt the value the array must be less than, n the nth time the array is less than the value

! output ltn: the value of the nth time the array falls below a certain value


function ltn(x,lt,n,num)
integer i,n,n1,idx,ltn,num
real x(num),lt

idx=0
n1=0


do i=1,num

if (x(i) < lt .AND. n1<n) then
idx=i
n1=n1+1
endif
enddo
ltn=idx

return
end






!!! return 1st, 2nd... nth, element of an array greater than some value

!!! gtn is an integer function i.e declare as integer i the calling program
!!! input num is the size of the 1d array x

function gtn(x,gt,n,num)
integer i2,n,n1,idx,gtn,num
real x(num),gt


idx=0
n1=0
!      write(*,*) 'gt',gt
i2=1
do while (i2 .le. num .and. n1<n)
!      write(*,*) x(i2)
if (x(i2) > gt) then
idx=i2
n1=n1+1

endif
i2=i2+1
enddo

gtn=idx
!      write(*,*) 'gt =',gtn,idx
return
end







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ran0(ran,idum)
INTEGER idum
REAL ran,AM
!       PARAMETER (IA=16807,IM=2147483647,AM=1./IM,
!     &  IQ=127773,IR=2836,MASK=123459876)



!!!“Minimal” random number generator of Park and Miller. Returns a uniform random deviate between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK) to initialize the sequence; idum must not be altered between calls for successive deviates in a sequence.
INTEGER k,i,IA,IM,IQ,IR,MASK
IA=16807
IM=2147483647
AM=1./IM
IQ=127773
IR=2836
MASK=123459876
stop
!      write(*,*) idum, MASK
idum=2
stop
i=ieor(idum,MASK)
idum=2
stop
!       write(*,*) idum, MASK,i
stop
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
stop
if (idum.lt.0) idum=idum+IM
ran=AM*idum
idum=ieor(idum,MASK)
return
END





FUNCTION ran2(idum)
INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
REAL ran2,AM,EPS,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,&
IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,&
IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,&
EPS=1.2e-7,RNMX=1.-EPS)
INTEGER idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/
if (idum.le.0) then
idum=max(-idum,1)
idum2=idum
do j=NTAB+8,1,-1
k=idum/IQ1

enddo

iy=iv(1)
endif

k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1)iy=iy+IMM1
ran2=min(AM*iy,RNMX)
return
end



FUNCTION ran3(idum)

INTEGER idum
INTEGER MBIG,MSEED,MZ
!	  REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!	  PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)

INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!	  REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if(idum.lt.0.or.iff.eq.0)then
iff=1
mj=abs(MSEED-abs(idum))
mj=mod(mj,MBIG)
ma(55)=mj
mk=1
do i=1,54
ii=mod(21*i,55)
ma(ii)=mk
mk=mj-mk
if(mk.lt.MZ)mk=mk+MBIG


mj=ma(ii)
enddo
do k=1,4
do i=1,55
ma(i)=ma(i)-ma(1+mod(i+30,55))
if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
enddo
enddo
inext=0
inextp=31
idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
return
END




!!uniform integer random number generator between a and b iseed = iseed
!! input a;r4,b;r4,iseed;i4
function iranu(a,b,iseed)

real:: u,a,b,ranu,ran3
integer iseed,iranu

if (iseed .lt. 0) call system_clock(iseed)
!write(*,*) iseed
ranu=ran3(iseed)*(int(b)-int(a)) + int(a)

!write(*,*) ranu,ran3(iseed),b,a,iranu
iranu=nint(ranu)
!write(*,*) ranu,ran3(iseed),b,a,iranu
return
end





!!uniform real random number generator
function ranu(a,b,iseed)

real:: u,a,b,ranu,ran3
integer iseed

ranu=ran3(iseed)*(b-a) + a

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ************************************************
function rang( avg, rms, iseed )
!! Update 07/08/2015: if iseed -ve use clock to get seed value, else get same random sequence again
!* Gaussian random number generator.
!* Input:
!*	avg	r4	mean value
!*	rms	r4	standard deviation
!*	iseed	i4	seed integer for ran( iseed )
!* Output:
!*	rang	r4	gaussian random number
!*	iseed	i4	seed integer from ran( iseed )
!*
!* Box-muller transform of 2 independent random variables
!* 1978 Peter Young @ Caltech -- original version GAUSS
!* 1995 Keith Horne @ St-Andrews -- adapted from GAUSS.
!* 2001 Sep KDH @ St-And - tidy up
!* changed to incorporate system clock iseed 2014 jan
logical	newpair /.true./
save y, u, newpair
rang = avg

if (iseed .lt. 0) call system_clock(iseed)

if( rms .le. 0. ) return
if ( newpair ) then
r2 = 10.
do while ( r2 .ge. 1. )
x = 2. * ran3( iseed ) - 1.
y = 2. * ran3( iseed ) - 1.
r2 = x * x + y * y
end do
u = sqrt( -2. * alog( r2 ) / r2 )
rang = x * u
else
rang = y * u
end if
newpair = .not. newpair
rang = avg + rms * rang
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! generate gausian random number
!	   subroutine rang(mu,sd,G,iseed)

!	   real:: u1,v,mu,sd,r,rsq,G,ran

!	   rsq = 2
!	   do while (rsq > 1)
!	   rsq = 0
!	   call ranu(0.,1.,ran)
!	   u1 = ran3(iseed)
!	   call ranu(0.,1.,ran)
!	   v=ran
!	   u1 = ran0(234651211)
!	   v = ran0(35165143)
!	   rsq = u1*u1 +v*v
!	   r = sqrt(rsq)
!	   !write(*,*) 'u1',u1,'v',v
!	   end do
!	   !write(*,*) 'rsq',rsq
!	   G=2.*u1/r*sqrt(-1.*log(r))*sd + mu

!	   end subroutine



!!!! 29th nov 2013 function to calculate the median value between adjacent points in a 1, D array

!!!!input t (1d real array), N points

!      program testcode
!      real t(4),op(4),medspac,med
!      integer i,N

!      N=4

!      do i=1,N
!      t(i)=1.*i
!      enddo
!      t(1)=6
!      t(2)=53

!      write(*,*) t

!      write(*,*) med(t,N)

!      end
!     call sort(t,N,op)
!      write(*,*) t
!      write(*,*) op


!      write(*,*) medspac(op,N)


!      end


!! function to calculate the median of 1d array t(N) Tested 19/12/2013

function med(t,N)

integer N
real t(N),op(N),med

call sort(t,N,op)

if (mod(N,2) ==0) then
med = 0.5*(op(N/2)+op(N/2+1))

else

med =op(N/2+1)

endif

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! return the median spacing between elements in a time series
!! input t(N) array of time spacings
function medspac(t,N)

integer i,N
real medspac, t(N), op(N)

call sort(t,N,op)

if (mod(N,2) ==0) then   !!! here there are an odd number of differences between adjacent array elements if there are an even number of elements in the array
medspac=t(N/2+1)-t(N/2)
else
medspac=op(ceiling(0.5*N)+1)-op(ceiling(0.5*N))
medspac=medspac+(op(ceiling(0.5*N))-op(ceiling(0.5*N)-1))
medspac=0.5*medspac
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Input: t,1d array of data size N. op output. the sorted, output array
subroutine sort(t,N,op)
real t(N), op(N),temp
integer N

op(1:N)=t(1:N)

do i =1,N-1
do j =i+1,N

if (op(i) .ge. op(j)) then
temp = op(j)
op(j)=op(i)
op(i)=temp
endif
enddo
enddo

return
end



!! returns the number of lines in a data file

function numline(c)
integer:: n,error,i,nlines,numline,lunit
character (len=*):: c



!	  n=0

!	  open (unit=2,file=c)
!     READ(12,*, iostat = error)
!      IF (error == -1) EXIT
!      n = n + 1
!      END D


!	  write(*,*) c
nlines = 0
OPEN (56, file = c)
DO
READ (56,*, END=10)
nlines = nlines + 1
END DO

10 CLOSE (56)
numline = nlines+1
!      write(*,*)'numline',numline,nlines
return
end



!!!! function to calculate chi sq
!      program testcs
!      real c(6),e(6),cmod(6)
!      integer N
!      real cisq,cs

!      N=6

!!      do i =1,N
!      e(i)=1
!      c(i)=6.
!      cmod(i)=6.

!      enddo
!      cmod(2)=7
!      write(*,*) c,cmod
!      write(*,*) cisq(c,cmod,e,N)
!      end



function cisq(x,xmod,er,N)
real x(N),xmod(N),er(N),chisq,cisq
integer N,i

sum = 0.d0
do i =1, N
a = x(i) - xmod(i)
sig = er(i)

sum = sum + a*a/(sig*sig)
!if (sum .ne. sum .or. 1./sum == 0) then
!write(*,*) 'infinity chisq'
!write(*,*) i,x(i),xmod(i),er(i)
!stop
!endif
enddo

cisq = real(sum)


return
end









!!function to fit "chisquared" but reject outliers by fitting slightly different
!!statistics for outlier > k sd from the model
!!iversion == 1, fit normal chi^2 for all points

!!iversion == 2, Automatically reject outliers more than k-sigma away using:
!!............chi_i -> min( chi_i, k ) [Makes each outlier pull on the model as if it were located at k sigma].

!!iversion == 3, similar to above but use smooth transition at k
!!............. chi_i^2 -> k^2 ( 1 - exp( -(chi_i/k)^2 )

!!DONT DO 4 not worked out yet
!!iversion == 4, make the BOF linear in |chi_i| rather than quadratic for |chi_i| > k
!!............. chi_i^2 => min( chi_i^2/k^2, |chi_i|/k )

!INPUT:: x(NIn)),xmod(Nin),er(Nin), input data model
! iversion 1,2,3
! sigk boundary over which rejection criteria apply (if |chi_i| > sigk apply rejection criteria)
!OUTPUT:: Nout: number of non rejected points
!........ modout(1) if modified 0 if original
subroutine cisqmod(x,xmod,er,Nin,iversion,sigk,Nout,modout,cisqsum)
integer,intent(in):: Nin,iversion
real,intent(in):: x(Nin), xmod(Nin), er(Nin),sigk
real,intent(Out):: cisqsum
integer modout(Nin)
double precision:: top, bot, a

modout(1:Nin) = 0
Nout = Nin
sigk2 = sigk*sigk

if (iversion == 1) then

bot = 0.d0
do it = 1, Nin
xnow = x(it)
xmodnow = xmod(it)
signow = er(it)
top = (xnow - xmodnow)/signow
bot = bot + top*top
enddo
cisqsum = bot

else if (iversion == 2) then

bot = 0.d0
do it = 1, Nin
xnow = x(it)
xmodnow = xmod(it)
signow = er(it)
top = abs((xnow - xmodnow)/signow)
bot = bot + min(top*top,sigk2)
if (top > sigk) then
modout(it) = 1
Nout = Nout - 1
endif
enddo
cisqsum = bot

else if (iversion == 3) then

bot = 0.d0
do it = 1, Nin
xnow = x(it)
xmodnow = xmod(it)
signow = er(it)
top = abs((xnow - xmodnow)/signow)
a = top/sigk
bot = bot + sigk*sigk * ( 1. - exp( -a*a ) )
if (top > sigk) then
modout(it) = 1
Nout = Nout - 1
endif
enddo
cisqsum = bot

else if (iversion == 4) then

bot = 0.d0
do it = 1, Nin
xnow = x(it)
xmodnow = xmod(it)
signow = er(it)
top = abs((xnow - xmodnow)/signow)
t2 = top*top
if (top .le. sigk) then
bot = bot + t2
else
b = sigk*(2*top - sigk)
bot = bot + b + alog(b/(t2))
modout(it) = 1
Nout = Nout - 1
endif
enddo

cisqsum = bot

endif



end subroutine





















!! subroutine to calculate the additional term in the BOF that is the penalty for
! expanding the error bars / applying outlier rejection
! FULL BOF (excluding fourier and additional priors)
! BOF_i = chi_i^2 + ln( 2 pi sigma_i )



!subroutine cisqmod_sigbof(x,xmod,er,Nin,iversion,sigk,Nout,modout,cisqsum)
!integer,intent(in):: Nin,iversion
!real,intent(in):: x(Nin), xmod(Nin), er(Nin),sigk
!real,intent(Out):: cisqsum
!integer modout(Nin)
!double precision:: top, bot, a
!
!
!
!
!
!
!end subroutine









!C*CGS ... physical constants in cgs units
!C+
REAL FUNCTION CGS( REQUEST )
!*
!* supplies value in cgs units of requested constant
!C--
!* 1988 Jul Keith Horne @ STScI
!* 1989 Jul Keith Horne @ STScI - added MH, ME, E
!* 1990 Jan KDH @ STScI - call double-precision version
CHARACTER*(*) REQUEST
REAL*8 DCGS
CGS = DCGS( REQUEST )
RETURN
END

!C*DCGS ... physical constants in double precision cgs units
!C+
REAL*8 FUNCTION DCGS( REQUEST )
!*
!* supplies value in cgs units of requested constant
!C--
!* Jul 1988 Keith Horne @ STScI
!* Jul 1989 Keith Horne @ STScI - added MH, ME, E
!* Jan 1990 KDH @ STScI - additional constants, double precision version
!* 2001 Apr KDH @ Austin - MJY and ANGSTROM
!* 2001 Jul KDH @ StAnd - g77 can't concatenate request
CHARACTER*(*)	REQUEST

!* WARNING: long names must preceed otherwise identical short ones
!* since otherwise wrong identification may be made

IF( REQUEST .EQ. 'PI' ) THEN
DCGS = 4.D0 * DATAN(1.D0)
ELSE IF( REQUEST .EQ. 'C' ) THEN
DCGS = 2.997925D10
ELSE IF( REQUEST .EQ. 'G' ) THEN
DCGS = 6.673D-8
ELSE IF( REQUEST .EQ. 'H' ) THEN
DCGS = 6.6262D-27
ELSE IF( REQUEST .EQ. 'KMS' ) THEN
DCGS = 1.E5
ELSE IF( REQUEST .EQ. 'K' ) THEN
DCGS = 1.3806D-16
ELSE IF( REQUEST .EQ. 'SIGMA' ) THEN
DCGS = 5.66956D-5
ELSE IF( REQUEST .EQ. 'MJY' ) THEN
DCGS = 1.E-26
ELSE IF( REQUEST .EQ. 'ANGSTROM' ) THEN
DCGS = 1.D-8

!* energies
ELSE IF( REQUEST .EQ. 'RYD') THEN
DCGS = 2.17992D-11
ELSE IF( REQUEST .EQ. 'EV') THEN
DCGS = 1.602192D-12
ELSE IF( REQUEST .EQ. 'E') THEN
DCGS = 4.80325D-10

!* earth
ELSE IF( REQUEST .EQ. 'MEARTH' ) THEN
DCGS = 5.976D27
ELSE IF( REQUEST .EQ. 'REARTH' ) THEN
DCGS = 6.378164D8
!* sun
ELSE IF( REQUEST .EQ. 'MSUN' ) THEN
DCGS = 1.989D33
ELSE IF( REQUEST .EQ. 'RSUN' ) THEN
DCGS = 6.9599D10
ELSE IF( REQUEST .EQ. 'LSUN' ) THEN
DCGS = 3.826D33
!* jupiter
ELSE IF( REQUEST .EQ. 'MJUP' ) THEN
DCGS = 3.1783D2 * 5.976D27
ELSE IF( REQUEST .EQ. 'RJUP' ) THEN
DCGS = 7.13D9

!* particles
ELSE IF( REQUEST .EQ. 'ME' ) THEN
DCGS = 9.10956D-28
ELSE IF( REQUEST .EQ. 'MP' ) THEN
DCGS = 1.672661D-24
ELSE IF( REQUEST .EQ. 'MH' ) THEN
DCGS = 1.67352D-24
ELSE IF( REQUEST .EQ. 'AMU' ) THEN
DCGS = 1.660531D-24

!* time
ELSE IF( REQUEST .EQ. 'YR' ) THEN
DCGS = 3.1556925D7

!* distance
ELSE IF( REQUEST .EQ. 'PC' ) THEN
DCGS = 3.085678D18
ELSE IF( REQUEST .EQ. 'AU' ) THEN
DCGS = 1.49597D13

!* help from the user
ELSE
10      WRITE(*,'(A,A,A,$)') &
' Enter CGS value for "', REQUEST, '" : '
READ(*,*,ERR=10 ) DCGS
END IF

!C	WRITE(*,*) 'CGS value for ", REQUEST, " is ', DCGS
RETURN
END




!!! uses avg.f90 from fortcode library

function rms(x,N)


implicit none

integer N,i
real x(N),rms,minv,xtemp,avg
double precision sum


rms=0.0
sum = 0.d0
minv=avg(x,N)
do i=1,N
xtemp=x(i)-minv
sum=xtemp*xtemp+sum
enddo

rms=sqrt(real(sum)/N)

return
end function





!********************************************
function zlum_ml( z, omega_m, omega_l )
!* luminosity distance is d = c * zlum / h0
pow = 3.
q0 = ( 0.5 * pow - 1.) * omega_m - omega_l
zlum_ml = zlum_q( z, q0 )
return
end





!********************************************
function zlum_q( z, q0 )
!* luminosity distance is d = c * zlum / h0
!* 2001 Apr Keith Horne @ Austin
!* from B.M.Peterson's AGN book
zq = z * q0
if( abs(zq) .lt. 0.5 ) then
zl = z * ( 1. + (z-zq) / ( 1. + sqrt(1.+2.*zq) + zq ) )
else if( q0.eq.0. ) then
zl = z * ( 1. + z / 2. )
!* negative q0 case from P.Mannheim paper
!c      else if( q0.eq.-1. ) then
!c        zl = z * ( 1. + z )
!c      else if( q0.lt.0. ) then
!c        opz = 1. + z
!c        zl = opz/(-q0) * ( opz - sqrt( opz**2 * (1. + q0) - q0 ) )
end if
zlum_q = zl
return
end




!C*BNU ... Planck function
!C+
FUNCTION BNU( WAVE, TEMP )
!*
!* input:
!*	WAVE	R4 wavelength in Angstroms
!*	TEMP	R4 Kelvin temperature
!* output:
!*	BNU	R4 blackbody intensity (erg/cm2/s/Hz/ster)
!C--
!* 1989 Jun Keith Horne @ STScI - clean up old subroutine
!* 2002 Aug KDH @ St.And - better long-wavelength limit
DATA C1/1.43883E8/
DATA C2/1.95722E5/
BNU = 0.
X = WAVE * TEMP
IF( X.LE.0. ) RETURN
X = C1 / X
IF( X.LT.1.E-4 ) THEN
FACTOR = 2. / ( X * ( X + 2. ) )
ELSE IF( X.LT.85. ) THEN
FACTOR = 1. / ( EXP( X ) - 1. )
ELSE
bnuln = 3. * alog( ( c1 / wave ) / c2 ) - x
bnu = exp( bnuln )
RETURN
END IF
X = X * TEMP / C2
BNU = FACTOR * X**3
RETURN
END



!*----------------------------------------
subroutine pgbgw
!* set background/foreground colours to white/black
!* 2010 Oct Keith Horne @ St-And
!call pgscr( 0, 1., 1., 1. )
!call pgscr( 1, 0., 0., 0. )
return
end



!C*LOOKUP -- inverse linear interpolation
!C+
SUBROUTINE LOOKUP( NPIX, DATA, TARGET, ILO, IHI, PART )
!*
!* Inverse linear interpolator. Finds ILO, IHI, PART such that
!*   TARGET = DATA(ILO) * (1.-PART) + DATA(IHI) * PART
!*
!* This version uses safe but slow method of one-pixel steps
!*
!* Input:
!*	NPIX	= Number of data points
!*	DATA	= Data values ( must be monotonic )
!*	TARGET	= Target value
!* Output:
!*	ILO	= Low pixel
!*	IHI	= Hi pixel
!*	PART	= fraction (0-1) for interpolation
!C--
!* Nov 1985 KDH @ STScI
!* 2001 Aug KDH @ St.And - single point, IBEG
REAL*4 DATA(*)
DATA NLAST,IHILAST,ILOLAST/1,1,1/

!* Null array
IF( NPIX.LE.0 ) THEN
ILO = 0
IHI = 0
PART = 0.
RETURN
END IF

!* Single point
IF( NPIX .EQ. 1 ) THEN
ILO = 1
IHI = 1
PART = 0.
RETURN
END IF

!* Use previous location
IF( NPIX .EQ. NLAST ) THEN
IF( DATA(IHILAST) .GE. TARGET .AND. &
DATA(ILOLAST) .LE. TARGET ) THEN
!     *	  DATA(ILOLAST) .LE. TARGET ) THEN
IHI = IHILAST
ILO = ILOLAST
GOTO 50
END IF
IBEG = ILOLAST
ELSE
IBEG = ( 1 + NPIX ) / 2
END IF

!* Determine uphill direction

IF( DATA(1) .LE. DATA(NPIX) ) THEN
ILOEND = 1
IHIEND = NPIX
IUP = 1
ELSE
ILOEND = NPIX
IHIEND = 1
IUP = -1
END IF

!* Locate an uphill point

DO I = IBEG, IHIEND, IUP
IF( DATA(I) .GE. TARGET ) GOTO 10
END DO
I = IHIEND
10	IHI = I

!* Locate nearest downhill point

DO I = IHI, ILOEND, -IUP
IF( DATA(I) .LE. TARGET ) GOTO 20
END DO
I = ILOEND
20	ILO = I

!* Locate nearest uphill point

DO I = ILO, IHI, IUP
IF( DATA(I) .GE. TARGET ) GOTO 30
END DO
I = IHI
30	IHI = I

!* Compute fractional part

50	PART = DATA(IHI) - DATA(ILO)
IF( PART.NE.0. ) PART = ( TARGET - DATA(ILO) ) / PART

NLAST = NPIX
ILOLAST = ILO
IHILAST = IHI

RETURN
END




!! program to take affine step in parameter space given covariance matrix
!! note the result is taken assuming mean equal to zero. You will have to add on either the mean value,
! or the current parameter values depending on how you want it.

!stepchange is a factor >0.0 default 1.0 by which the magnitude of the affine step should be altered
!to suit the required acceptance criteria and optimise MCMC
!icovcalc if 1 then covariance matrix to be calculated in this routine, else
!.... cov is an input to the routine and must be calculated prior to calling.
! if icmean is on then step from the mean of all past iterations else step from previous iteration
subroutine affine_step(NP,Niter,p_past,cov,pnew,icovcalc,stepchange,icmean)

integer NP,Niter,icovcalc,icmean
real pnew(NP),p_past(Niter,NP), cov(NP,NP),&
eval(NP), evec(NP,NP), ave(NP),gaus_0_v(NP),u_gaus_0_v(NP,NP),stepchange
iseed=332423532
!icmean = 1

if (icmean .eq. 0) then
!! calculate means
do it = 1,NP
sum = 0.d0
do idx = 1,Niter
sum = sum + p_past(idx,it)
enddo
ave(it) = sum/Niter
enddo
!!
else
do it = 1,NP
ave(it) = p_past(Niter,it)
enddo
endif

!calculate covariance matrix if required
if (icovcalc == 1) then
call covpar_2018(NP,Niter,p_past,cov)
endif


!! calculate Eigen values and vectors of covariance matrix
!! v(np,np) eigen vectors v(1:3,i) are the eigen vectors for eigen value d(i)


call jacobi(cov,np,np,eval,evec,nrot)

!doesnt like negatvie eigen values normally means you have put two identical parameters in
!the parameter saved matrix (e.g stepping tr slope and iradiation slope together)
do idx = 1,np
evnow = eval(idx)
eval(idx)= stepchange*max(0.0,evnow)
enddo



!write(*,*) eval
!write(*,*)
!write(*,*) evec

!evec(1,1) = 0.4985
!evec(2,1) = -0.8669
!evec(1,2) = -0.8669
!evec(2,2) = -0.4985
!stop

!calculate N(0,diag(eval))
do ipc = 1,NP
gaus_0_v(ipc) = rang(0.,sqrt(eval(ipc)),iseed)
enddo


!!note u_gaus_0_v has first entry as columns and second as rows to speed things up
!.. u_gaus_0_v = u*N(0,lam) u is matrix f evecs and lam is diag matrix of evt (evt just vector in this code)
do ipr = 1,NP
evt = gaus_0_v(ipr)
do ipc = 1,NP
u_gaus_0_v(ipc,ipr) = evt*evec(ipc,ipr)
enddo
enddo

!write(*,*) 'evs1', evec(:,1), eval(1)
!write(*,*) 'evs2', evec(:,2), eval(2)


!! add together results
do ipc = 1,NP
sum = 0.d0
do ipr = 1,NP
sum = sum + u_gaus_0_v(ipc,ipr)
enddo
pnew(ipc) = sum + ave(ipc)
if (pnew(ipc) .eq. p_past(niter,ipc)) then
!do it = 1,niter
!write(*,*) p_past(it,ipc)
!enddo
write(*,*) 'affine step parameter didnt move, stuck'
write(*,*) 'Number of iterations used to build covariance matrix', Niter
write(*,*) 'eigen values',eval(ipc)
write(*,*) 'eigen vectors',evec(ipc,1:NP)
write(*,*) 'gaus_0_v',gaus_0_v(1:NP)
!read(*,*)
endif

!write(*,*)
!write(*,*) 'Number of iterations used to build covariance matrix', Niter
!write(*,*) 'eigen values',eval(ipc)
!write(*,*) 'eigen vectors',evec(ipc,1:NP)
!write(*,*) 'average, old -->new params ',ave(ipc),p_past(niter,ipc),pnew(ipc)
!write(*,*) 'gaus_0_v',gaus_0_v(1:NP)
!write(*,*)
enddo


!write(*,*) 'avs,',ave(1:NP)
end subroutine
!if (ipr .eq. ipc) part1(ipr,ipc) =

!a = part1(ipc,ipr)
!a = a +


!!! program to test the affine_step using test rotation matrix
!gflags affine_step.f90 ../../covpar.f90 ../../ran_new2.f90 ../../cgs.f90 ../../itp2.for ../../eigensolve_symetric.f90 ../../numline.for
!
!program affinetest
!
!real,allocatable:: dat(:,:), cov(:,:), covdup(:,:),datnew(:,:), datop(:,:), dat_t(:,:)
!
!npar = 2
!nstep = 100
!nits = numline('affine_step_test.dat') - 1
!
!allocate(dat(nits,npar), cov(npar,npar), covdup(npar, npar), datnew(nstep,npar), dat_t(npar,nits))
!
!
!open(unit = 1,file = 'affine_step_test.dat')
!do i = 1,nits
!read(1,*) dat(i,1:npar)
!dat_t(1:npar,i) = dat(i,1:npar)
!enddo
!close(1)
!
!!! calculate the covariance matrix
!call covpar(npar,10000,dat_t,cov)
!covdup = cov
!
!!covdup(1,1) = sqrt(cov(1,1))
!!covdup(2,2) = sqrt(cov(2,2))
!
!
!!! produce a new sample stepping in the covariance of the parameter (should produce a similar display of results)
!do it = 1,nstep
!!write(*,*) dat(it,1:2)
!call affine_step(npar,nits,dat,cov,datnew(it,1:npar),1,1.0,0)
!cov = covdup
!enddo
!
!write(*,*) 'covmat',cov
!!!
!
!open(unit = 1,file = 'affine_step_test_op.dat')
!do i = 1,nstep
!write(1,*) datnew(i,1:npar)
!enddo
!close(1)
!
!
!end program
!
!
!



!! Input a(np,np) matrix,
!! n set = np
!! Output d(np) eigen values
!! v(np,np) eigen vectors v(1:3,i) are the eigen vectors for eigen value d(i)
!! nrot number of iterations required to solve


SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER n,np,nrot,NMAX
REAL a(np,np),d(np),v(np,np)
PARAMETER (NMAX=500)
!Computes all eigenvalues and eigenvectors of a real symmetric matrix


INTEGER i,ip,iq,j
REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)




do ip=1,n

do iq=1,n
v(ip,iq)=0.
enddo

v(ip,ip)=1.
enddo

do ip=1,n
b(ip)=a(ip,ip)

d(ip)=b(ip)
z(ip)=0.

enddo

nrot=0
do i=1,50
sm=0.
do ip=1,n-1

do iq=ip+1,n
sm=sm+abs(a(ip,iq))
enddo

enddo

if(sm.eq.0.)return

if(i.lt.4)then
tresh=0.2*sm/n**2

else

tresh=0.

endif
do ip=1,n-1
do iq=ip+1,n
g=100.*abs(a(ip,iq))
!After four sweeps, skip the rotation if the o-diagonal element is small.
if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
a(ip,iq)=0.
else if(abs(a(ip,iq)).gt.tresh) then
h=d(iq)-d(ip)

if (abs(h)+g.eq.abs(h)) then
t=a(ip,iq)/h
else
theta=0.5*h/a(ip,iq)
t=1./(abs(theta)+sqrt(1.+theta**2))
if(theta.lt.0.)t=-t
endif

c=1./sqrt(1+t**2)
s=t*c
tau=s/(1.+c)
h=t*a(ip,iq)
z(ip)=z(ip)-h
z(iq)=z(iq)+h
d(ip)=d(ip)-h
d(iq)=d(iq)+h
a(ip,iq)=0.
do j=1,ip-1

g=a(j,ip)
h=a(j,iq)
a(j,ip)=g-s*(h+g*tau)
a(j,iq)=h+s*(g-h*tau)
enddo

do j=ip+1,iq-1

g=a(ip,j)
h=a(j,iq)
a(ip,j)=g-s*(h+g*tau)
a(j,iq)=h+s*(g-h*tau)
enddo

do j=iq+1,n

g=a(ip,j)
h=a(iq,j)
a(ip,j)=g-s*(h+g*tau)
a(iq,j)=h+s*(g-h*tau)
enddo


do j=1,n
g=v(j,ip)
h=v(j,iq)
v(j,ip)=g-s*(h+g*tau)
v(j,iq)=h+s*(g-h*tau)
enddo

nrot=nrot+1
endif
enddo

enddo

do ip=1,n
b(ip)=b(ip)+z(ip)
d(ip)=b(ip)
z(ip)=0.

enddo

enddo

write(*,*) 'too many iterations in jacobi'
read(*,*)
return
END subroutine






!program testabove
!
!real a(3,3), d(3), v(3,3)
!
!
!n = 3
!np = 3
!a(1,1) = 0
!a(2,1) = 1
!a(3,1) = -1
!
!a(1,2) = 1
!a(2,2) = 2
!a(3,2) =1
!
!a(1,3) = -1
!a(2,3) = 1
!a(3,3) = 2
!
!
!
!call jacobi(a,n,np,d,v,nrot)
!
!
!write(*,*) 'eigen values and vectors'
!
!
!do i = 1,3
!write(*,*) i, d(i), v(1:3,i)
!enddo
!
!
!end program



