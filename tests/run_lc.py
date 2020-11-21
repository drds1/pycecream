import pycecream
import astropy_stark.myfake as mf
import matplotlib.pylab as plt
import os




class test_pc:

    def __init__(self):
        '''input arguments'''
        self.fake_wavelength = [-1.0,4680.0,4686.0,4720.,4720.,7480.0,7760.,7764.0,7764.0]
        self.fake_snr = [10.,50.0,30.0,50.,30.,50.0,50.0,50.0,30.0]
        self.fake_cadence = [1.0]*len(self.fake_snr)

    def gen_fake(self):
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
        thcent = 20.0
        )
        self.dat = synthetic_data['echo light curves']


    def prep_pycecream(self,test_project_folder = 'test_pycecream_output'):
        '''
        test pycecream using yasamans script
        :return:
        '''
        cream_lc0, cream_lc1, cream_lc4, cream_lc2, cream_lc3, cream_lc8, cream_lc5, cream_lc6, cream_lc7 = self.dat

        #instantiate and remove previous test if present
        os.system('rm -rf '+test_project_folder)
        a = pycecream.pycecream()
        a.project_folder = test_project_folder

        #step accretion rate?
        a.p_accretion_rate_step = 0.1
        a.bh_mass = 6.6e8

        # MgII Line lightcurve
        a.add_lc(cream_lc0, name='line 0  (MgII)', kind='line',background_polynomials=[0.1,0.1])
        a.p_linelag_centroids_step = 0.0
        # g-band photometric lightcurves
        a.add_lc(cream_lc1,name='continuum (Bok)', kind='continuum', wavelength = 4680)
        a.add_lc(cream_lc2,name='continuum 4720 (CFHT 1)',kind='continuum', wavelength  = 4720, share_previous_lag=True)
        a.add_lc(cream_lc3,name='continuum 4720 (CFHT 2)',kind='continuum', wavelength = 4720, share_previous_lag=True)
        a.add_lc(cream_lc4,name='continuum 4686  (SynthPhot)',kind='continuum', wavelength = 4686, share_previous_lag=True)
        # i-band photometric lightcurves
        a.add_lc(cream_lc5,name='continuum (Bok)', kind='continuum', wavelength= 7760, share_previous_lag = False)
        a.add_lc(cream_lc6,name='continuum (CFHT 1)',kind='continuum', wavelength = 7764, share_previous_lag=True)
        a.add_lc(cream_lc7,name='continuum (CFHT 2)',kind='continuum', wavelength = 7764, share_previous_lag=True)
        a.add_lc(cream_lc8,name='continuum (SynthPhot)',kind='continuum', wavelength = 7480,share_previous_lag=True)
        a.hi_frequency = 0.5
        a.N_iterations = 20
        self.pc = a


    def run_pycecream(self):
        '''
        run prepared pycecream instance
        :return:
        '''
        self.pc.run(ncores = 4)

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


if __name__ == '__main__':
    x = test_pc()
    x.gen_fake()
    x.prep_pycecream(test_project_folder = 'test_pycecream_output_sub')

    '''
    set step sizes and add priors
    '''
    # accretion rate
    x.pc.p_accretion_rate = 0.1
    x.pc.p_accretion_rate_step = 0.3
    x.pc.p_accretion_rate_priorcentroid = 1.2
    x.pc.p_accretion_rate_priorwidth = 0.004

    # inclination
    x.pc.p_inclination = 0.0
    x.pc.p_inclination_step = 0.1
    x.pc.p_inclination_priorcentroid = 30.0
    x.pc.p_inclination_priorwidth = 0.0001

    # viscous slope
    x.pc.p_viscous_slope = 0.8
    x.pc.p_viscous_slope_step = 0.1
    x.pc.p_viscous_slope_priorcentroid = 0.75
    x.pc.p_viscous_slope_priorcentroid = 0.0001

    # irradiation slope
    x.pc.p_irradiation_slope = 0.81
    x.pc.p_irradiation_slope_step = 0.11
    x.pc.p_irradiation_slope_priorcentroid = 0.75
    x.pc.p_irradiation_slope_priorcentroid = 0.0001

    #dial up number of iterations
    x.pc.N_iterations = 1000

    # turn on priors (once we confirm this works the priors will be turned on by default)
    x.pc.custom_priors = False

    x.run_pycecream()
    x.post_run()
    lcop = x.output_lightcurves
    chains = x.output_chains