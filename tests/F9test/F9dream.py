import pandas as pd
import pycecream as pc
import os
import pickle
import matplotlib.pylab as plt
import numpy as np
import corner
from matplotlib.backends.backend_pdf import PdfPages


if __name__ == '__main__':

    '''
    1: Load F9 light curves and arrange into nice dictionary
    '''
    #load data and assemble into pycecream format
    raw = pd.read_csv('scratch_F9lightcurves.csv')
    raw['MJD'] = raw['MJD'] - 58000
    raw['TelObs'] = raw['Tel'] + raw['Obs']
    groups = {}
    means = {}
    sd_mean = {}
    unique_filters = list(raw['Filter'].unique())

    LightcurveDictionary = {'name':[],
                      'group':[],
                      'lightcurve':[]}
    LightcurveSummary = {'group':[],'names':[],'means':[],'sd_means':[]}
    for g in unique_filters:
        rawgroup = raw[raw['Filter'] == g]
        unique_telobs = list(rawgroup['TelObs'].unique())
        LightcurveSummwary['group'].append(g)
        gnames = []
        gmeans = []
        for u in unique_telobs:
            lightcurve = rawgroup[rawgroup['TelObs'] == u][['MJD', 'Flux', 'Error']].sort_values(by='MJD')
            gnames.append(u)
            gmeans.append(lightcurve.values[:,1].mean())
            LightcurveDictionary['name'].append(u)
            LightcurveDictionary['group'].append(g)
            LightcurveDictionary['lightcurve'].append(lightcurve.values)
        LightcurveSummary['names'].append(gnames)
        LightcurveSummary['means'].append(gmeans)
        LightcurveSummary['sd_means'].append(np.std(gmeans))

    #this tells us which filter groups have the greates misalignment (in terms of means)
    LightcurveSummary = pd.DataFrame(LightcurveSummary).sort_values(by='sd_means',ascending = False)

    print('testing on group...',test_group)
    picklefile = 'test_output.pickle'


    '''
    2: setup and run a pycecream.dream instance for the light curve with
    most misaligned telescope points defined by standard deviation of means
    '''
    test_group = LightcurveSummary['group'].iloc[0]

    # Prepare pycecream
    dream = pc.dream(Niterations = 200)
    for (name, group, lcdat) in zip(LightcurveDictionary['name'],
                                         LightcurveDictionary['group'],
                                         LightcurveDictionary['lightcurve']):
        if group == test_group:
            dream.add_lc(lcdat, name, errorbar_variance=True, errorbar_rescale=True)

    dream.run()

    
    #save simulation as pickle output
    os.system('rm ' + picklefile)
    pickle_out = open(picklefile, "wb")
    pickle.dump(dream, pickle_out)
    pickle_out.close()


    #load previous simulation
    pickle_in = open(picklefile, "rb")
    dream = pickle.load(pickle_in)


    '''
    3: Make diagnostic plots for the merged light curves
    and mcmc parameters
    '''
    #make diagnostic plot
    with PdfPages('dream_diagnostic.pdf') as pdf:
        #combine the merged and input light curve plots into a single figure
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        fig, ax1 = dream.plot_input_individual(fig_in = fig, ax_in = ax1)
        ax1.legend(fontsize='xx-small')
        ax2 = fig.add_subplot(2, 1, 2)
        fig, ax1 = dream.plot_merged_individual(fig_in=fig, ax_in=ax2)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        #plot the mcmc chains and covariances
        output_chains = dream.MCMC_chains
        parms = ['offset','stretch',
                 'noise m','noise var']
        parms_nicenames = ['vertical offset','vertical stretch',
                           'error bar multiplicative correction','error bar extra variance']
        #plot the covariances
        for Parm, ParmNicename in zip(parms,parms_nicenames):
            fig = plt.figure()
            fig = dream.plot_chains_or_covariances(output_chains,Parm, ParmNicename, type = 'covariances', fig_in = fig)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

        #plot the chains
        for Parm, ParmNicename in zip(parms,parms_nicenames):
            fig = plt.figure()
            fig = dream.plot_chains_or_covariances(output_chains,Parm, ParmNicename, type = 'chain', fig_in = fig)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

