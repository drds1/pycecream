import pycecream
import astropy_stark.myfake as mf
import matplotlib.pylab as plt
import unittest
import os


def get_synthetic_data():

    finished = False
    '''input arguments'''
    fake_wavelength = [4680.0,4686.0,4720.,4720.,
                       7480.0,7760.,7764.0,7764.0,-1.0]
    fake_snr = [50.0,30.0,50.,30.,50.0,50.0,50.0,30.0,10.]
    fake_cadence = [1.0]*len(fake_snr)



    synthetic_data = mf.myfake(
    fake_wavelength,
    fake_snr,
    fake_cadence,
    thcent = 20.0
    )
    dat = synthetic_data['echo light curves']

    name = ['continuum 4680',
            'continuum 4686',
            'continuum 4720',
            'continuum 4720  (Telescope 2)',
            'continuum 7480',
            'continuum 7760',
            'continuum 7764',
            'continuum 7764 (Telescope 2)',
            'line 0  (MgII)'
            ]

    share_previous_lag = [False,False,False,True,
                          False,False,False,True,
                          False]

    kind = ['continuum']*(len(dat)-1) + ['line']
    test_data = {'lightcurve':dat,
                 'name':name,
                 'kind':kind,
                 'wavelength':fake_wavelength,
                 'share_previous_lag':share_previous_lag}

    return test_data




class Test_synthetic_data(unittest.TestCase):


    def test_synthetic_data(self):

        finished = False
        '''input arguments'''
        synthetic_data = get_synthetic_data()

        a = pycecream.pycecream()
        a.project_folder = 'test_pycecream'

        #step accretion rate?
        a.p_accretion_rate_step = 0.1
        a.p_linelag_centroids_step = 0.0

        ndata = len(synthetic_data['name'])

        for i in range(ndata):
            a.add_lc(synthetic_data['lightcurve'][i],
                     name=synthetic_data['name'][i],
                     kind = synthetic_data['kind'][i],
                     wavelength=synthetic_data['wavelength'][i],
                     share_previous_lag=synthetic_data['share_previous_lag'][i])
        a.hi_frequency = 0.5
        a.N_iterations = 20
        a.run()


        '''
        post run
        '''
        output_chains = a.get_MCMC_chains(location=None)
        output_lightcurves = a.get_light_curve_fits(location=None)
        '''
        Check the input settings are ok prior to running
        '''
        print(a.lightcurve_input_params)
        finished = True
        self.assertEqual(finished,True)
        os.system('rm -rf test_pycecream')


if __name__ == '__main__':
    unittest.main()