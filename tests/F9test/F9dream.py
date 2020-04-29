import pandas as pd
import pycecream as pc
import os
import pickle
import matplotlib.pylab as plt



if __name__ == '__main__':

    #load data and assemble into pycecream format
    raw = pd.read_csv('scratch_F9lightcurves.csv')
    raw['MJD'] = raw['MJD'] - 58000
    raw['TelObs'] = raw['Tel'] + raw['Obs']
    groups = {}
    unique_filters = list(raw['Filter'].unique())

    LightcurveDictionary = {'name':[],
                      'group':[],
                      'lightcurve':[]}
    for g in unique_filters:
        rawgroup = raw[raw['Filter'] == g]
        unique_telobs = list(rawgroup['TelObs'].unique())
        groups[g] = []
        for u in unique_telobs:
            groups[g].append(u)
            lightcurve = rawgroup[rawgroup['TelObs']==u][['MJD','Flux','Error']].sort_values(by='MJD')
            LightcurveDictionary['name'].append(u)
            LightcurveDictionary['group'].append(g)
            LightcurveDictionary['lightcurve'].append(lightcurve.values)

    picklefile = 'test_output.pickle'


    # Prepare pycecream
    pc = pc.dream()
    for (name, group, lcdat) in zip(LightcurveDictionary['name'],
                                         LightcurveDictionary['group'],
                                         LightcurveDictionary['lightcurve']):
        if group == 'g':
            pc.add_lc(lcdat, name, errorbar_variance=True, errorbar_rescale=True)

    pc.run()
    ax1,fig = pc.plot_merged()
    #pc.add_lc()
    
    #save pickle output
    os.system('rm ' + picklefile)
    pickle_out = open(picklefile, "wb")
    pickle.dump(pc, pickle_out)
    pickle_out.close()


    pickle_in = open(picklefile, "rb")
    pc = pickle.load(pickle_in)
    x = pc.x

    x = pc._dream__combine_individual_output_lightcurves()
    df_merged_out = x['combined_output']
    df_merged_in = x['combined_input']

    plt.close()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(df_merged_in[:,0],df_merged_in[:,1],label='input')
    ax1.scatter(df_merged_out[:, 0], df_merged_out[:, 1],label='merged')
    plt.legend()
    plt.show()



