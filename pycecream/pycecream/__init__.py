import numpy as np
import pandas as pd
import os
import glob


class pycecream:
    '''
    One stop shop for fitting time lags and response functions to AGN
    accretion discs. Fitting continuum light curves, pycecream can infer the
    inclination and temperature profile of the AGN disc by fitting the wavelegnth
    dependent response functions described in
    Starkey et al 2016 ( https://ui.adsabs.harvard.edu/#abs/arXiv:1511.06162 )
    For a full list of creams features, please see the (sadly out of date) mamnual that
    describes these features as applied to the previous fortran version of the code
    CREAM. A more up-to-date manual will follow shortly.

    Global class instance arguments...

    redshift: The target redshift (default 0).
    high_frequency: Large numbers will explore higher frequency variations at the expense
    of computation time.

    '''
    def __init__(self):
        '''
        The default parameters are given below. First up are the global non-fitted input parameters.
        Second are the global fitted parameters whos starting values can be modified below.
        Note that...
        1) Entering zero step sizes indicates the parameter will not be stepped)
        2) Light curve-specific parameters (the error bar scaling, starting lag centroids and widths
        for line light curves etc) are specified in the add_lc function.
        '''


        #config, path, compiler parameters
        self.module_path = os.path.dirname(os.path.realpath(__file__))
        self.fortran_compile_command = 'gfortran cream_f90.f90 -o creamrun.exe'
        print('pycecream path... ' + self.module_path)

        #convention parameters
        self.output_directory = 'pycecream'
        self.append_date_to_output_directory = False

        #global non-fitted parameters
        self.redshift = 0.0
        self.high_frequency = 0.5
        self.bh_mass = 1.e7
        self.bh_efficieny = 0.1
        self.N_iterations = 1000
        self.lag_lims = [-10.0,50.0]

        #fitted parameters
        self.p_inclination = 0.0
        self.p_inclination_step = 0.0
        self.p_inclination_priorcentroid = None
        self.p_inclination_priorwidth = None
        self.p_accretion_rate = 0.1
        self.p_accretion_rate_step = 0.0
        self.p_accretion_rate_priorcentroid = None
        self.p_accretion_rate_priorwidth = None
        self.p_viscous_slope = 0.75
        self.p_viscous_slope_step = 0.0
        self.p_extra_variance_step = 0.1

        #non configureable parameters
        self.count_lightcurves = 0
        self.count_continuum_lightcurves = 0
        self.count_line_lightcurves = 0
        self.p_linelag_centroids_start = 0.0
        self.p_linelag_centroids_step = 5.0
        self.p_linelag_widths_start    = 2.0
        self.p_linelag_widths_step = 0.0
        self.dir_pwd = os.getcwd()

        self.lightcurve_input_params = pd.DataFrame(columns=[
                            'name', 'type', 'wavelength', 'noise model',
                            'share previous lag','temporary file name',
                            'mean', 'standard deviation', 'tophat centroid',
                            'tophat centroid step', 'tophat width',
                            'tophat width step'
        ])
        self.global_input_params = pd.DataFrame(columns = [
            'redshift','BH mass','BH efficiency','upper fourier frequency',''
        ])



    def setup_directory_structure(self):
        '''
        Set up the output directory structure for a cream simulation.
        Shouild be called once as you add the first light curve
        :return:
        '''
        #make a directory into which to store the cream results
        dir_pycecream = self.output_directory
        if self.append_date_to_output_directory is True:
            dir_pycecream = self.output_directory+'_'+str(pd.datetime.today().strftime("%d_%m_%Y"))+'_'
            child_dirs = next(os.walk('.'))[1]
            number_of_pyceream_dirs = len( [c for c in child_dirs if dir_pycecream in c] )
            dir_pycecream = dir_pycecream+str(number_of_pyceream_dirs)

        self.dir_pycecream = dir_pycecream
        self.dir_sim = 'simulation_files'
        os.mkdir(self.dir_pycecream)
        os.mkdir(self.dir_pycecream+'/'+self.dir_sim)

        #copy fortran files to pycecream directory
        os.system('cp '+self.module_path+'/cream_f90.f90 ./'+self.dir_pycecream)
        os.system('cp ' + self.module_path + '/creaminpar.par ./' + self.dir_pycecream)

    def add_lc(self,input,
               kind = 'line',
               wavelength = -1.0,
               expand_errors = ['var','multiplicative'],
               name = None,
               share_previous_lag = False):
        '''
        This is the go to command to add a new light curve into the
        simulation.
        :param input: either a Nx3 numpy array of time,flux,errors or a link to a file in the same format
        :param kind: 'line' or 'continuum'. If continuum, must specify
        wavelegnth
        :param wavelength: The centroid wavelength of the contuum light curve
        :param expand_errors:
        :param share_previous_lag:
        :param name: optional to set name for a light curve to annotate on plots and in data frames
        :return:
        '''
        #set up the directory structure if first call
        if self.count_lightcurves == 0:
            self.setup_directory_structure()

        #count the numnber of line or continuum light curves already specified
        if kind is 'line':
            count = self.count_line_lightcurves
            self.count_line_lightcurves = self.count_line_lightcurves + 1
        elif kind is 'continuum':
            count = self.count_continuum_lightcurves
            self.count_continuum_lightcurves = self.count_continuum_lightcurves + 1
            if wavelength == -1.0:
                raise Exception('Must specify wavelength for a continuum light curve')
        else:
            raise Exception('kind argument must be "line" or "continuum"')

        #configure the naming convention
        if name is None:
            name_ann = kind + ' lightcurve ' + np.str(count)
        else:
            name_ann = name

        #load the data and save in required directory
        if type(input) is str:
            dat = np.loadtxt(input)
        elif type(input) is np.ndarray:
            dat = np.array(input)
        else:
            raise Exception('input to add_lc must be file name or numpy.ndarray')
        fname = kind+'_'+np.str(count)+'.dat'

        #check the data for problems and save
        check_for_bad_values(dat,name_ann)
        np.savetxt(self.dir_pycecream+'/'+self.dir_sim+'/'+fname,dat)


        '''
        configure the line lag settings (if sharing the same response function as the previous line
        should then use the same step size else increment by a small positive number e.g 0.1       
        '''
        tophat_centroid      = self.p_linelag_centroids_start
        tophat_centroid_step = self.p_linelag_centroids_step
        tophat_width         = self.p_linelag_widths_start
        tophat_width_step    = self.p_linelag_widths_step

        if share_previous_lag is False:
            tophat_centroid_step = tophat_centroid_step + 0.1*self.count_lightcurves
        else:
            tophat_centroid_step = self.lightcurve_input_params['tophat centroid step'].values[-1]


        #update the lightcurve_input_params table of records
        df = pd.DataFrame(data = [name_ann,kind,wavelength,expand_errors,share_previous_lag,fname,
                                  np.mean(dat[:,1]), np.std(dat[:,1]),
                                  tophat_centroid,
                                  tophat_centroid_step,
                                  tophat_width,
                                  tophat_width_step
                                  ],
                     index=['name', 'type', 'wavelength', 'noise model',
                            'share previous lag','temporary file name',
                            'mean', 'standard deviation', 'tophat centroid',
                            'tophat centroid step', 'tophat width',
                            'tophat width step']).T

        self.lightcurve_input_params = pd.DataFrame(pd.concat([self.lightcurve_input_params,df]))
        self.lightcurve_input_params['wavelength']= \
            pd.to_numeric(self.lightcurve_input_params['wavelength'],downcast = 'float')
        self.lightcurve_input_params['mean'] = \
            pd.to_numeric(self.lightcurve_input_params['mean'],downcast = 'float')
        self.lightcurve_input_params['standard deviation'] = \
            pd.to_numeric(self.lightcurve_input_params['standard deviation'],downcast = 'float')
        self.lightcurve_input_params['tophat centroid']= \
            pd.to_numeric(self.lightcurve_input_params['tophat centroid'],downcast = 'float')
        self.lightcurve_input_params['tophat centroid step'] = \
            pd.to_numeric(self.lightcurve_input_params['tophat centroid step'],downcast = 'float')
        self.lightcurve_input_params['tophat width'] = \
            pd.to_numeric(self.lightcurve_input_params['tophat width'],downcast = 'float')
        self.lightcurve_input_params['tophat width step'] = \
            pd.to_numeric(self.lightcurve_input_params['tophat width step'],downcast = 'float')
        self.count_lightcurves = self.count_lightcurves + 1




    def set_creaminpar(self):
        '''
        configure the creaminpar.par fortran input file
        :return:
        '''
        # lines for each of the input parameters in creaminpar.par (dont change)
        idcos = 37
        idmdot = 31
        idslope = 70
        idsig = 22
        idnits = 10
        idplot = 4
        idlaglim = 13
        idhifreq = 9

        #modify and write the creaminpar.par file as needed
        with open(self.dir_pycecream+'/creaminpar.par') as f:
            content = f.readlines()
            content = [x.strip() for x in content]
            f.close()
            content[0] = './'+np.str(self.dir_sim)
            f.close()
        content[idcos] = np.str(self.p_inclination)
        content[idcos+1] = np.str(self.p_inclination_step)
        content[idmdot] = np.str(self.p_accretion_rate)
        content[idmdot+3] = np.str(self.p_accretion_rate_step)
        content[idslope] = np.str(self.p_viscous_slope)
        content[idslope+2] = np.str(self.p_viscous_slope_step)
        content[idhifreq] = np.str(self.high_frequency)
        content[idnits] = np.str(self.N_iterations)
        content[idlaglim] = np.str(self.lag_lims[0]) + ' ' + np.str(self.lag_lims[1]) + '       ! lower and upper lag limits'
        content[idplot] = np.str(int(np.ceil(self.N_iterations / 4.0)))

        #deal with multiplicative step sizes
        step_multiplicative = list(self.lightcurve_input_params['noise model'].map(set(['multiplicative']).issubset))
        step = []
        for i in range(self.count_lightcurves):
            if step_multiplicative[i] is True:
                step.append(0.1)
            else:
                step.append(0.0)
        a = ''.join([np.str(step[i]) + ' ' for i in range(self.count_lightcurves)])
        content[idsig] = a
        a = ''.join([np.str(1.0) + ' ' for i in range(self.count_lightcurves)])
        content[idsig - 1] = a

        #write updated creaminpar.par file
        f = open(self.dir_pycecream+'/creaminpar.par', 'w')
        for fn in content:
            f.write(fn + '\n')
        f.close()

    def set_creamnames(self):
        '''
        set the creamnames.dat file summarising the file names and wavelengths
        required by fortran
        :return:
        '''
        f = open(self.dir_pycecream+'/'+self.dir_sim+'/creamnames.dat','w')
        for i in range(self.count_lightcurves):
            f.write("'"+self.lightcurve_input_params['temporary file name'].values[i]+"' "+
                    np.str(self.lightcurve_input_params['wavelength'].values[i])+"\n")

        f.close()


    def set_tophat(self):
        '''
        creates the file instructing cream how to treat the top hat response functions
        (which light curves to use the same response function, starting parameter values etc)
        :return:
        '''
        f = open(self.dir_pycecream+'/'+self.dir_sim + '/cream_th.par', 'w')
        for i in range(self.count_lightcurves):
            f.write(
                np.str(self.lightcurve_input_params['tophat centroid'].values[i]) + ' ' +
            np.str(self.lightcurve_input_params['tophat centroid step'].values[i]) + ' 0.0 -1.0 '+
            np.str(self.lightcurve_input_params['tophat width'].values[i]) + ' ' +
            np.str(self.lightcurve_input_params['tophat width step'].values[i]) + ' 0.0 -1.0\n'
            )
        f.close()




    def set_var(self):
        '''
        creates the file instructing cream how to treat the extra variance parameters
        (starting parameter values and step sizes)
        :return:
        '''
        f = open(self.dir_pycecream+'/'+self.dir_sim + '/cream_var.par', 'w')
        for i in range(self.count_lightcurves):
            f.write('0.0 ')
            f.write(np.str( self.p_extra_variance_step*
                            self.lightcurve_input_params['standard deviation'].values[i]) +'\n')
        f.close()



    def run(self):
        '''
        run the cream code. Make sure input above is correct first
        :return:
        '''

        #set the required fortran inoput files with the information entered
        self.set_creaminpar()
        self.set_creamnames()
        self.set_tophat()
        self.set_var()

        #compile and run
        os.chdir(self.dir_pycecream)
        os.system(self.fortran_compile_command)
        os.system('./creamrun.exe')
        os.chdir(self.dir_pwd)


    def get_simulation_dir(self,location=None):
        '''
        return the specific cream directory containing the
        outputpars.par file (this is one level deeper than
        self.dir_pycecream)
        :return:
        '''
        if location is None:
            simulation_dir = self.dir_pycecream
        else:
            simulation_dir = location
        return(simulation_dir)


    def get_MCMC_chains(self,location = None):
        '''
        Load the results from the previous cream simulation.
        If location is None, load the most recent simulation
        :return:
        '''
        #locate the simulation results
        simulation_dir = self.get_simulation_dir(location=location)
        lcnames = self.lightcurve_input_params['name']
        results_dir = glob.glob(simulation_dir + '/simulation_files/output_2*')[0]
        '''
        Begin loading the simulation results from the various stored locations
        '''

        #load the disk parameters
        dat_output = np.loadtxt(results_dir + '/outputpars.dat')
        p_output_names = ['Mdot','cos i','Tr_alpha']
        p_output = dat_output[:,[2,3,4]]

        #load the stretch, offset and multiplicative error bar rescaling parameters
        dat_output = np.loadtxt(results_dir + '/outputpars2.dat')
        p_output_names = p_output_names + ['stretch '+l for l in lcnames] + \
                         ['offset '+l for l in lcnames] + \
                         ['noise m '+l for l in lcnames]
        p_output = np.hstack((p_output,dat_output[:,:3*self.count_lightcurves]))

        #load the extra variance noise parameters
        dat_output = np.loadtxt(results_dir + '/outputpars_varexpand.dat')
        p_output_names = p_output_names + \
                         ['noise var '+l for l in lcnames]
        p_output = np.hstack((p_output,dat_output))

        #load the top hat centroid parameters
        dat_output = np.loadtxt(results_dir + '/outputpars_th.dat')
        p_output_names = p_output_names + \
                         ['top hat centroid '+l for l in lcnames] + \
                         ['top hat width '+l for l in lcnames]
        p_output = np.hstack((p_output,dat_output[:,:2*self.count_lightcurves]))

        #save all results for each MCMC iteration to a pandas data frame
        self.output_parameters = pd.DataFrame(data = p_output,columns = p_output_names)

        return(self.output_parameters)




    def get_light_curve_fits(self,location=None):
        '''
        Load the fitted light curves
        :return:
        '''
        #locate the simulation results
        simulation_dir = self.get_simulation_dir(location=location)
        results_dir = glob.glob(simulation_dir + '/simulation_files/output_2*')[0]
        '''
        Begin loading the simulation results from the various stored locations
        '''

        '''
        the merged models for each wavelength
        '''
        count = 0
        output_merged_model = {}
        for tf in self.lightcurve_input_params['temporary file name'].values:
            dat = np.loadtxt(results_dir+'/plots/merged_mod_'+tf+'.dat')
            name = self.lightcurve_input_params['name'].values[count]
            if count == 0:
                output_merged_model['time']=dat[:,0]
            output_merged_model[name+' model'] = dat[:,1]
            output_merged_model[name+' uncerts'] = dat[:, 2]
            count = count + 1
        self.output_merged_model = pd.DataFrame(output_merged_model)

        '''
        the merged data points for each file (this just takes the
        original light curves and applies the errorbar rescaling,
        vertical and horizontal shifts to get similar wavelengths
        onto the same scale
        '''
        count = 0
        output_merged_data = {}
        for tf in self.lightcurve_input_params['temporary file name'].values:
            dat = np.loadtxt(results_dir+'/plots/merged_dat_'+tf+'.dat')
            name = self.lightcurve_input_params['name'].values[count]
            if count == 0:
                output_merged_data['time']=dat[:,0]
            output_merged_data[name+' data'] = dat[:,1]
            output_merged_data[name+' uncerts'] = dat[:, 2]
            count = count + 1
        self.output_merged_data = output_merged_data

        '''
        load the unmerged model (easier). This is just the model (with error bar rescalings)
        evaluated at each wavelength without rescaling relative to reference light curve
        '''
        dat = np.loadtxt(results_dir+'/plots/modellc.dat')
        dat_sig = np.loadtxt(results_dir + '/plots/modellc_sig.dat')
        count = 0
        output_model = {'time':dat[:,0]}
        for name in self.lightcurve_input_params['name'].values:
            output_model[name+' model'] = dat[:,count+1]
            output_model[name + ' uncerts'] = dat_sig[:, count + 1]
            count = count + 1
        self.output_model = pd.DataFrame(output_model)

        '''
        condense into a dictionary for output
        '''
        function_output = {'model':self.output_model,
                           'merged model':self.output_merged_model,
                           'merged data': self.output_merged_data}
        return(function_output)









def check_for_bad_values(dat,name_ann):
    '''
    Check an input light curve for bad values and raise exceptions if found
    :param dat:
    :param name_ann:
    :return:
    '''
    time_range = np.max(dat[:, 0]) - np.min(dat[:, 0])
    y = dat[:, 0]
    error_bars = dat[:, 1]
    if time_range == 0:
        raise Exception('light curve '+name_ann+' has no time range')
    if np.std(y) != np.std(y):
        raise Exception('light curve '+name_ann+' has a bad (nan,inf) value')
    if 0 in error_bars:
        raise Exception('light curve '+name_ann+' has a zero error bar')
