import pycecream
import astropy_stark.myfake as mf
import matplotlib.pylab as plt
import os
import numpy as np
import pickle
from matplotlib.backends.backend_pdf import PdfPages


class test_pc:

    def __init__(self):
        '''input arguments'''
        self.fake_wavelength = [4720.,7480.0]
        self.fake_snr = [10.]*len(self.fake_wavelength)
        self.fake_cadence = [1.0]*len(self.fake_snr)
        self.ncuts = 3

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
        thcent = 20.0,thi=100,
        )
        self.dat = synthetic_data['echo light curves']

    def transform_fake(self, offset = 10, sd = 0,plot=True):
        '''
        simulate noise
        :return:
        '''
        ncuts = self.ncuts
        self.datnorm = {'name':[],
                        'wavelength':[],
                        'light curve':[]}


        idxwavelength = 0
        for d in self.dat:
            ndat = len(d)
            y = d[:, 1]
            y = (y - np.mean(y)) / np.std(y)
            d[:, 1] = y
            d[:, 2] = d[:, 2] / np.std(y)

            #randomly select points for each simulated calibration offset
            random_order = np.arange(ndat)
            np.random.shuffle(random_order)
            selected_points = np.array_split(random_order,ncuts)

            #apply artificial offset for each chunk of points
            for i in range(ncuts):
                idx = selected_points[i]
                newmean = i*offset
                newsd = i*sd + 1
                y = d[idx,1]
                ynew = (y - 0)/1*newsd + newmean
                sdnew = d[idx,2]/1*newsd
                d[idx,1] = ynew
                d[idx, 2] = sdnew
                title = 'wavelength'+str(int(self.fake_wavelength[idxwavelength])+1)+'_telescope'+str(i+1)
                self.datnorm['name'].append(title)
                self.datnorm['wavelength'].append(int(self.fake_wavelength[idxwavelength]))
                self.datnorm['light curve'].append(d[idx,:])
            idxwavelength += 1
            if plot is True:
                plt.scatter(d[:,0],d[:,1])
        if plot is True:
            plt.show()




    def setup_pycecream(self,test_project_folder = 'test_pycecream_output'):
        '''
        test pycecream using yasamans script
        :return:
        '''

        #instantiate and remove previous test if present
        os.system('rm -rf '+test_project_folder)
        a = pycecream.pycecream()
        a.project_folder = test_project_folder

        #step accretion rate?
        a.p_accretion_rate_step = 0.1
        a.bh_mass = 6.6e8

        #add the light curves one at a time
        previous_wavelength = np.nan
        wavelengths = self.datnorm['wavelength']
        names = self.datnorm['name']
        LightCurveData = self.datnorm['light curve']
        for i in range(len(names)):
            wavelength, lc, name = wavelengths[i],  LightCurveData[i], names[i]
            if wavelength == previous_wavelength:
                share_previous_lag = True
            else:
                share_previous_lag = False
            a.add_lc(lc,name=name,wavelength=wavelength, share_previous_lag=share_previous_lag)

        # MgII Line lightcurve
        #a.add_lc(cream_lc0, name='line 0  (MgII)', kind='line',background_polynomials=[0.1,0.1])
        #a.p_linelag_centroids_step = 0.0
        ## g-band photometric lightcurves
        #a.add_lc(cream_lc1,name='continuum (Bok)', kind='continuum', wavelength = 4680)
        #a.add_lc(cream_lc2,name='continuum 4720 (CFHT 1)',kind='continuum', wavelength  = 4720, share_previous_lag=True)
        #a.add_lc(cream_lc3,name='continuum 4720 (CFHT 2)',kind='continuum', wavelength = 4720, share_previous_lag=True)
        #a.add_lc(cream_lc4,name='continuum 4686  (SynthPhot)',kind='continuum', wavelength = 4686, share_previous_lag=True)
        ## i-band photometric lightcurves
        #a.add_lc(cream_lc5,name='continuum (Bok)', kind='continuum', wavelength= 7760, share_previous_lag = False)
        #a.add_lc(cream_lc6,name='continuum (CFHT 1)',kind='continuum', wavelength = 7764, share_previous_lag=True)
        #a.add_lc(cream_lc7,name='continuum (CFHT 2)',kind='continuum', wavelength = 7764, share_previous_lag=True)
        #a.add_lc(cream_lc8,name='continuum (SynthPhot)',kind='continuum', wavelength = 7480,share_previous_lag=True)
        a.hi_frequency = 0.5
        a.N_iterations = 20
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


if __name__ == '__main__':

    newsim = False
    x = test_pc()
    x.gen_fake()
    x.transform_fake(plot=False)


    '''
    setup and run new lightcurve-merging test
    '''
    #only do these steps if running a new simulation (takes times)
    #for diagnostic plots and analysis just reload previous
    if newsim is True:
        #setup pycecream object
        pc = x.setup_pycecream(test_project_folder='test_pycecream_output')

        #run pycecream
        pc.run(ncores = 1)

        #save to pickle for reloading later
        pickle_out = open('merge_test_pcObj.pickle', "wb")
        pickle.dump(pc, pickle_out)
        pickle_out.close()
    else:
        pickle_in = open('merge_test_pcObj.pickle', "rb")
        pc = pickle.load(pickle_in)



    '''
    Post-simulation analysis
    gather merged light curves, compare with inputs, generate diagnostic plots
    '''
    #gather simulation outputs
    output_chains = pc.get_MCMC_chains(location=None)
    output_lightcurves = pc.get_light_curve_fits(location=None)


    #compare the input light curves with merged light curves
    input_lightcurves = {}
    for i in range(len(x.datnorm['wavelength'])):
        input_lightcurves[x.datnorm['name'][i]] = x.datnorm['light curve'][i]
    merged_lightcurves = output_lightcurves['merged data']

    
    #plot one wavelength at a time
    wavelengths = np.array(x.datnorm['wavelength'])
    unique_wavelengths = np.unique(wavelengths)
    nwavs = len(unique_wavelengths)

    #save diagnostics as multipage pdf
    with PdfPages('merge_test_diagnostic.pdf') as pdf:
        fig = plt.figure()
        idx = 1
        for w in unique_wavelengths:
            ax1 = fig.add_subplot(nwavs,1,idx)
            ax1.set_ylabel('flux')
            ax1.set_xlabel('time')
            ax1.set_title(str(w)+'Å')

            #assemble the input and output (merged) lightcurves for each wavelength
            lc_idx = np.where(wavelengths == w)[0]
            datmerged = np.zeros((0,3))
            datinput = np.zeros((0, 3))
            for lci in lc_idx:
                name = x.datnorm['name'][lci]
                datmerged = np.vstack([datmerged,merged_lightcurves[name]])
                lcin = x.datnorm['light curve'][lci]
                datinput = np.vstack([datinput,lcin])
                ax1.scatter(lcin[:, 0], lcin[:, 1], label=name)
            datmerged = datmerged[np.argsort(datmerged[:,0]),:]
            datinput = datinput[np.argsort(datinput[:, 0]), :]

            #plot the merged light curves for each wavelength
            ax1.plot(datmerged[:, 0], datmerged[:, 1], label='merged')
            ax1.legend(fontsize='xx-small')
            idx += 1
        plt.tight_layout()
        pdf.savefig()
        plt.close()



        #now plot the trace plots for the stretch and offset parameters
        #on a new page
        fig = plt.figure()
        idx = 1
        for w in unique_wavelengths:
            ax1 = fig.add_subplot(nwavs,2,idx)
            ax1.set_xlabel('Iteration')
            ax1.set_title('Offset Parameter\n'+str(w)+'Å')
            ax2 = fig.add_subplot(nwavs, 2, idx+1)
            ax2.set_xlabel('Iteration')
            ax2.set_title('Vertical Stretch Parameter\n' + str(w) + 'Å')
            for lci in lc_idx:
                name = x.datnorm['name'][lci]
                ax1.plot(output_chains['offset '+name].values,label=name)
                ax2.plot(output_chains['stretch ' + name].values, label=name)
            idx += 2
        pdf.savefig()
        plt.close()

        #now plot the trace plots for the error bar noise parameters
        #on a third page
        fig = plt.figure()
        idx = 1
        for w in unique_wavelengths:
            ax1 = fig.add_subplot(nwavs, 2, idx)
            ax1.set_xlabel('Iteration')
            ax1.set_title('Multiplicative Noise Parameter\n' + str(w) + 'Å')
            ax2 = fig.add_subplot(nwavs, 2, idx + 1)
            ax2.set_xlabel('Iteration')
            ax2.set_title('Extra Variance Noise Parameter\n' + str(w) + 'Å')
            for lci in lc_idx:
                name = x.datnorm['name'][lci]
                ax1.plot(output_chains['noise m ' + name].values, label=name)
                ax2.plot(output_chains['noise var ' + name].values, label=name)
            idx += 2
        pdf.savefig()
        plt.close()







