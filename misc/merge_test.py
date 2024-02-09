import pycecream
import astropy_stark.myfake as mf
import astropy_stark.myedlum as eddington
import matplotlib.pylab as plt
import os
import numpy as np
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import corner
import pandas as pd
import itertools
import multiprocessing as mp

class test_pc:

    def __init__(self, wavelengths = [4720.,7480.0],
                 snr = [10.,10.],
                 cadence = [1.,1,],
                 BHMass = 10000000.0,
                 EddRat=0.1,
                 Mdot = None,
                 BHefficiency = 0.1,
                 time_baseline = 100,
                 line_lag = 20,
                 disk_inc = 0.0,
                 disk_TXslope = 0.75):
        '''input arguments'''
        self.fake_wavelength = wavelengths
        self.fake_snr = snr
        self.fake_cadence = cadence
        self.BHMass = BHMass
        self.EddRat = EddRat
        self.BHefficiency = BHefficiency
        self.time_baseline = time_baseline
        self.line_lag = line_lag
        self.disk_inc = disk_inc
        self.disk_TXslope = disk_TXslope
        if Mdot is None:
            self.Mdot = eddington.ermin_mdotout(BHMass, EddRat, eta=BHefficiency)
        else:
            self.EddRat = eddington.edd(BHMass, Mdot, eta=BHefficiency)



    def gen_fake(self,iseed = -1):
        '''
        make fake lightcurves
        :return:
        mf.myfake arguments are
        wavelengths: enter the wavelengths (-1 indicates an emission line light curve modelled with a top-hat response),
        snr: set the signal-to-noise relative to light curve rms
        cadence:set the mean cadence
        top hat centroid: set the centroid for the top-hat (I think thats what this does but the line lag
        thing is still newish so Im used to just making continuum light curve)
        '''


        synthetic_data = mf.myfake(
            self.fake_wavelength,
            self.fake_snr,
            self.fake_cadence,
            embh = self.BHMass,
            er = self.EddRat,
            eta = self.BHefficiency,
            thcent = self.line_lag,
            degi = self.disk_inc,
            thi=self.time_baseline,
            iseed = iseed,
            sx = self.disk_TXslope
        )
        self.dat = synthetic_data['echo light curves']

    def transform_fake(self, offset, sd,
                       noise_m,
                       noise_var,plot=True):
        '''
        simulate noise
        :return:
        '''
        ncuts = len(offset)
        self.datnorm = {'name':[],
                        'wavelength':[],
                        'light curve':[],
                        'offset true':[],
                        'stretch true':[],
                        'noise m true':[],
                        'noise var true':[],
                        'disk true':[self.Mdot,np.cos(np.pi/180*self.disk_inc), self.disk_TXslope]}

        idxwavelength = 0
        for d in self.dat:
            #randomly select points for each simulated calibration offset
            random_order = np.arange(len(d))
            np.random.shuffle(random_order)
            selected_points = np.array_split(random_order,ncuts)

            #apply artificial offset for each chunk of points
            for i in range(ncuts):
                idx = np.sort(selected_points[i])
                newmean = offset[i]
                newsd = sd[i]
                yold = d[idx,1]
                print(yold.shape, yold.mean(), yold.std(), newsd, newmean)
                ynew = (yold - yold.mean())/yold.std()*newsd + newmean
                sdold = d[idx,2]/yold.std()
                #assume original errors were underestimated by the inputted noise parameters
                sdnew = np.sqrt((sdold**2 - noise_var[i])/noise_m[i]**2)
                d[idx,1] = ynew
                d[idx, 2] = sdnew
                title = 'wavelength'+str(int(self.fake_wavelength[idxwavelength])+1)+'_telescope'+str(i+1)
                self.datnorm['name'].append(title)
                self.datnorm['wavelength'].append(int(self.fake_wavelength[idxwavelength]))
                self.datnorm['light curve'].append(d[idx,:])
                self.datnorm['offset true'].append(newmean)
                self.datnorm['stretch true'].append(newsd)
                self.datnorm['noise m true'].append(noise_m[i])
                self.datnorm['noise var true'].append(noise_var[i])
            idxwavelength += 1
            if plot is True:
                plt.scatter(d[:,0],d[:,1])
        if plot is True:
            plt.show()




    def setup_pycecream(self,test_project_folder = 'test_pycecream_output', dream = False):
        '''
        test pycecream using yasamans script
        :return:
        '''

        #instantiate and remove previous test if present
        os.system('rm -rf '+test_project_folder)
        a = pycecream.pycecream()
        a.project_folder = test_project_folder

        #set the hyperparameters
        a.p_accretion_rate_step = 0.1
        a.bh_mass = self.BHMass
        a.p_inclination = self.disk_inc
        a.bh_efficieny = self.BHefficiency
        a.p_inclination_step = 0.0
        a.hi_frequency = 0.5
        a.N_iterations = 20

        #add the light curves one at a time
        previous_wavelength = np.nan
        wavelengths = self.datnorm['wavelength']
        names = self.datnorm['name']
        LightCurveData = self.datnorm['light curve']


        if dream is False:
            for i in range(len(names)):
                wavelength, lc, name = wavelengths[i],  LightCurveData[i], names[i]
                if wavelength == previous_wavelength:
                    share_previous_lag = True
                else:
                    share_previous_lag = False
                a.add_lc(lc,kind = 'continuum',
                         name=name,wavelength=wavelength,
                         share_previous_lag=share_previous_lag)
                previous_wavelength = wavelength

        elif dream is True:
            wav_tophat = 0.0
            for i in range(len(names)):
                wavelength, lc, name = wavelengths[i],  LightCurveData[i], names[i]
                if wavelength == previous_wavelength:
                    share_previous_lag = True
                else:
                    share_previous_lag = False
                    wav_tophat += -1.0
                a.add_lc(lc,kind='line',name=name,wavelength=np.round(wav_tophat,1), share_previous_lag=share_previous_lag)
                previous_wavelength = wavelength

        return a





    def post_run(self):
        '''
        analyse output
        :return:
        '''
        self.output_chains = self.pc.get_MCMC_chains(location=None)
        self.output_lightcurves = self.pc.get_light_curve_fits(location=None)
        '''
        Check the input settings are ok prior to running
        '''
        print(self.pc.lightcurve_input_params)


def _plotpage(parameter_names, xdatnorm, parameter_nicenames, output_chains):
    '''
    script-specific plotting routines for diagnosis
    :param names:
    :param nicenames:
    :param unique_wavelengths:
    :param output_chains:
    :return:
    '''
    wavelengths = np.array(xdatnorm['wavelength'])
    unique_wavelengths = np.unique(wavelengths)
    nwavs = len(unique_wavelengths)
    fig = plt.figure()
    idx = 1
    npars = len(parameter_names)
    for w in unique_wavelengths:
        lc_idx = np.where(wavelengths == w)[0]
        for i2 in range(npars):
            ax1 = fig.add_subplot(nwavs, npars, idx+i2)
            ax1.set_xlabel('Iteration')
            ax1.set_title(parameter_nicenames[i2]+'\n' + str(w) + 'Å')
            for lci in lc_idx:
                color = next(ax1._get_lines.prop_cycler)['color']
                name = xdatnorm['name'][lci]
                ax1.plot(output_chains[parameter_names[i2] + name].values, label=name, color=color)
                try:
                    true = xdatnorm[parameter_names[i2]+'true'][lci]
                    print(name+' '+parameter_names[i2]+'true value = '+str(true))
                    ax1.plot([0,len(output_chains)],[true]*2,label=None,ls=':',color = color)
                except:
                    print(parameter_names[i2]+'true: key not found in xdatnorm')
                    pass

        idx += npars




def _postsim_plot(pc, x, simulation_name):
    '''
        Post-simulation analysis
        gather merged light curves, compare with inputs, generate diagnostic plots
        '''
    # gather simulation outputs
    output_chains = pc.get_MCMC_chains(location=None)
    output_lightcurves = pc.get_light_curve_fits(location=None)

    # compare the input light curves with merged light curves
    input_lightcurves = {}
    for i in range(len(x.datnorm['wavelength'])):
        input_lightcurves[x.datnorm['name'][i]] = x.datnorm['light curve'][i]
    merged_lightcurves = output_lightcurves['merged data']

    # diagnose number of points in input and output light curves
    # previously there was a bug where the fake data was being
    # regenerated but not the simulation
    # making the input lightcurves unrelated to output lightcurves
    # this checks for same numner of points in each input and output dataset
    names = list(input_lightcurves.keys())
    num_points = {'names': names,
                  'inputs': [],
                  'outputs': []}
    for n in names:
        num_points['inputs'].append(len(input_lightcurves[n]))
        num_points['outputs'].append(len(merged_lightcurves[n]))
    num_points = pd.DataFrame(num_points)

    # plot one wavelength at a time
    wavelengths = np.array(x.datnorm['wavelength'])
    unique_wavelengths = np.unique(wavelengths)
    nwavs = len(unique_wavelengths)

    # save diagnostics as multipage pdf
    with PdfPages(simulation_name + '_merge_test_diagnostic.pdf') as pdf:
        fig = plt.figure()
        idx = 1
        for w in unique_wavelengths:
            ax1 = fig.add_subplot(nwavs, 2, idx)
            ax1.set_ylabel('flux')
            ax1.set_xlabel('time')
            ax1.set_title(str(w) + 'Å\n Input')

            # assemble the input and output (merged) lightcurves for each wavelength
            lc_idx = np.where(wavelengths == w)[0]
            datmerged = np.zeros((0, 3))
            datinput = np.zeros((0, 3))
            for lci in lc_idx:
                name = x.datnorm['name'][lci]
                datmerged = np.vstack([datmerged, merged_lightcurves[name]])
                lcin = x.datnorm['light curve'][lci]
                datinput = np.vstack([datinput, lcin])
                ax1.errorbar(lcin[:, 0], lcin[:, 1], lcin[:, 2], ls='', label=name)

            datmerged = datmerged[np.argsort(datmerged[:, 0]), :]
            datinput = datinput[np.argsort(datinput[:, 0]), :]
            ax1.legend(fontsize='xx-small')

            # plot the merged light curves for each wavelength
            ax2 = fig.add_subplot(nwavs, 2, idx + 1)
            ax2.set_ylabel('flux')
            ax2.set_xlabel('time')
            ax2.set_title(str(w) + 'Å\n Merged')
            ax2.errorbar(datmerged[:, 0], datmerged[:, 1], datmerged[:, 2], ls='', marker=None,
                         color='k', label='merged')
            ax2.errorbar(datmerged[:, 0], datmerged[:, 1], datinput[:, 2], ls='', marker=None,
                         color='grey', label=None)
            ax2.legend(fontsize='xx-small')
            idx += 2
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # parms
        parms = ['offset ', 'stretch ', 'noise m ', 'noise var ', 'disk ']
        parms_nicenames = ['Offset Parameter', 'Vertical Stretch Parameter',
                           'Multiplicative Noise Parameters',
                           'Extra Variance Noise Parameters',
                           'Accretion Disk Parameters']
        # now plot the trace plots for the stretch and offset parameters
        # on a new page
        _plotpage(parms[:2],
                  x.datnorm,
                  parms_nicenames[:2],
                  output_chains)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # now plot the trace plots for the error bar noise parameters
        # on a third page
        _plotpage(parms[2:4],
                  x.datnorm,
                  parms_nicenames[2:4],
                  output_chains)

        # make covariance corner plots
        colnames_output_chains = list(output_chains.columns)
        datnormcolumns = list(x.datnorm.keys())
        for Parm, ParmNicename in zip(parms, parms_nicenames):
            corner_columns = [c for c in colnames_output_chains if Parm in c]
            new_corner_columns = [c.replace(Parm, '') for c in corner_columns]
            df = output_chains[corner_columns].copy()
            df.columns = new_corner_columns

            # Isolate parameters optimised in MCMC chain
            # if non of the current parameter set are varied then skip
            dfstd = df.std()
            names_VariedColumns = list(dfstd[dfstd > 1.e-10].index)
            if len(names_VariedColumns) > 0:

                # add the truths if present
                true_parms_all_columns = np.array([c for c in datnormcolumns if Parm in c])
                if len(true_parms_all_columns) > 0:
                    true_parms_column = true_parms_all_columns[0]
                    truths_Varied = x.datnorm[true_parms_column]
                else:
                    truths_Varied = None

                # make the corner covariance plot and add title
                print('making corner plots for ', names_VariedColumns)
                fig = corner.corner(df[names_VariedColumns], plot_contours=False, truths=truths_Varied)
                fig.suptitle('Covariance Plots: ' + ParmNicename, fontsize=16)
                fig.tight_layout()
                pdf.savefig()
                plt.close()

def run_one_sim(inputs):
    '''
    setup and run new lightcurve-merging test
    '''
    offset, sd, noise_m, noise_var, simulation_name, dream = inputs


    wavelengths = [4720., 7480.0]
    cadence = [1.0,1.0]
    snr = [10., 10.]
    BHMass = 1.e8
    EddRat = 0.1
    Mdot = None
    N_iterations = 20
    BHefficiency = 0.1
    time_baseline = 100
    line_lag = 20
    iseed = 1234
    disk_inc = 0.0
    disk_TXslope = 0.75
    newsim = True

    # only do these steps if running a new simulation (takes times)
    # for diagnostic plots and analysis just reload previous
    if newsim is True:
        # generate fake data. Specify wavelengths, snr, cadence
        x = test_pc(wavelengths=wavelengths,
                    snr=snr,
                    cadence=cadence,
                    BHMass=BHMass,
                    EddRat=EddRat,
                    Mdot=Mdot,
                    BHefficiency=BHefficiency,
                    time_baseline=time_baseline,
                    line_lag=line_lag,
                    disk_inc=disk_inc,
                    disk_TXslope=disk_TXslope)
        x.gen_fake(iseed = iseed)

        # generate three 'telescopes' for each wavelength with vertical offset, vertical scaling,
        # multiplicative error bar and additive error bar (variance) calibration differences
        x.transform_fake(offset=offset, sd=sd,
                         noise_m=noise_m,
                         noise_var=noise_var, plot=False)

        # setup pycecream object
        pc = x.setup_pycecream(test_project_folder=simulation_name,dream = dream)

        # save the fake data for later use
        pc.x = x

        #set number of iterations
        pc.N_iterations = N_iterations

        # run pycecream
        pc.run(ncores=1)

        # save to pickle for reloading later
        file = simulation_name + '.pickle'
        os.system('rm ' + file)
        pickle_out = open(file, "wb")
        pickle.dump(pc, pickle_out)
        pickle_out.close()
    else:
        file = simulation_name + '.pickle'
        pickle_in = open(file, "rb")
        pc = pickle.load(pickle_in)
        x = pc.x


    _postsim_plot(pc, x, simulation_name)


    return pc




if __name__ == '__main__':

    '''
    run several tests varying each of the offset, vertical scaling, multiplicative and variance noise terms
    in turn 
    1) offsets
    2) vertical scaling
    3) multiplicative noise
    4) additive noise
    
    #default values (change globally in run_one_sim function)
    wavelengths = [4720., 7480.0],
    snr = [10., 10.],
    cadence = [1., 1, ],
    offset = [0, 10, 20],
    sd = [1, 1, 1],
    noise_m = [1, 1, 1],
    noise_var = [0, 0, 0],
    BHMass = 1.e8,
    EddRat = 0.1,
    Mdot = None,
    BHefficiency = 0.1,
    time_baseline = 100,
    line_lag = 20,
    iseed = -1,
    disk_inc = 0.0,
    disk_TXslope = 0.75,
    outputdir = 'simulation_results',
    simulation_name = 'sim_1',
    newsim = True
    
    '''


    outputdirectory = 'scratch_dream_tests'
    dream = True

    offset_default = [0,0,0]
    Alloffset = [[0, 10, 20]] +[offset_default]*3

    sd_default = [1,1,1]
    Allsd = [sd_default] +[[1,2,3]] + [sd_default]*2

    noise_m_default = [1,1,1]
    Allnoise_m = [noise_m_default]*2 + [[1,2,3]] + [noise_m_default]

    noise_var_default = [0,0,0]
    Allnoise_var = [noise_var_default]*3 + [[5,10,15]]

    simnames = ['scratch_offset','scratch_stretch','scratch_noisem','scratch_noisevar']


    #run simulations in parallel
    N = 1#len(Alloffset)
    sim_settings = [(Alloffset[i], Allsd[i], Allnoise_m[i], Allnoise_var[i], outputdirectory+'/'+simnames[i], dream) for i in range(N)]

    #'''
    os.system('rm -rf '+outputdirectory)
    os.system('mkdir '+outputdirectory)
    jobs = []
    for i in range(N):
        p = mp.Process(target=run_one_sim, args=(sim_settings[i],))
        jobs.append(p)
        p.start()
    # wait for jobs to finish before continuing
    for j in jobs:
        j.join()
    #'''


    #load the output
    results = []
    for i in range(N):
        file = outputdirectory+'/'+simnames[i] + '.pickle'
        pickle_in = open(file, "rb")
        results.append(pickle.load(pickle_in))

    pc = results[0]
    x = pc.x

    _postsim_plot(pc, x, outputdirectory+'/'+simnames[0])