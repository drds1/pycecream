import pandas as pd
import pycecream as pc





if __name__ == '__main__':

    #load data and assemble into pycecream format
    raw = pd.read_csv('scratch_F9lightcurves.csv')
    raw['MJD'] = raw['MJD'] - 58000
    filter_groups = {
        'HX':1,
        'SX':2,
        'UVW2':3,
        'UVM2':4,
        'UVM1':5,
        'UVW1':6,
        'U':7,
        'B':8,
        'V':9,
        'u':10,
        'g':11,
        'r':12,
        'i':13,
        'z':14
    }

    raw['group'] = raw['Filter'].replace(filter_groups)
    raw['TelObs'] = raw['Tel'] + raw['Obs']

    unique_groups = list(raw['group'].unique())

    LightcurveDictionary = {'name':[],
                      'group':[],
                      'lightcurve':[]}
    for g in unique_groups:
        rawgroup = raw[raw['group'] == g]
        unique_telobs = list(rawgroup['TelObs'].unique())
        for u in unique_telobs:
            lightcurve = rawgroup[rawgroup['TelObs']==u][['MJD','Flux','Error']].sort_values(by='MJD')
            LightcurveDictionary['name'].append(u)
            LightcurveDictionary['group'].append(g)
            LightcurveDictionary['lightcurve'].append(lightcurve.values)




    # Prepare pycecream
    pc = pc.dream()
    for (name, group, lcdat) in zip(LightcurveDictionary['name'],
                                         LightcurveDictionary['group'],
                                         LightcurveDictionary['lightcurve']):
        if group == 1:
            pc.add_lc(lcdat, name, errorbar_variance=True, errorbar_rescale=True)

    pc.run()
    #pc.add_lc()








