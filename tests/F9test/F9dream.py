import pandas as pd
import pycecream as pc
import os
import pickle
import matplotlib.pylab as plt
import numpy as np



if __name__ == '__main__':

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
        LightcurveSummary['group'].append(g)
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
    test_group = LightcurveSummary['group'].iloc[0]
    print('testing on group...',test_group)
    picklefile = 'test_output.pickle'

    '''
    # Prepare pycecream
    pc = pc.dream()
    for (name, group, lcdat) in zip(LightcurveDictionary['name'],
                                         LightcurveDictionary['group'],
                                         LightcurveDictionary['lightcurve']):
        if group == test_group:
            pc.add_lc(lcdat, name, errorbar_variance=True, errorbar_rescale=True)

    pc.run()

    
    #save pickle output
    os.system('rm ' + picklefile)
    pickle_out = open(picklefile, "wb")
    pickle.dump(pc, pickle_out)
    pickle_out.close()
    '''


    pickle_in = open(picklefile, "rb")
    pc = pickle.load(pickle_in)


    #combine the merged and input light curve plots into a single figure
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    fig, ax1 = pc.plot_input(fig_in = fig, ax_in = ax1)

    ax2 = fig.add_subplot(2,1,2)
    fig, ax2 = pc.plot_merged(fig_in = fig, ax_in = ax2)

    plt.show()
