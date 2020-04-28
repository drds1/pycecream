import pandas as pd
import pycecream as pc





if __name__ == '__main__':

    #load data and assemble into pycecream format
    raw = pd.read_csv('scratch_F9lightcurves.csv')
    raw['MJD'] = raw['MJD'] - 58000
    filter_wavelengths = {
        'HX':4,
        'SX':5,
        'UVW2':100,
        'UVM2':100,
        'UVM1':100,
        'UVW1':100,
        'U':1,
        'B':2,
        'V':3,
        'u':4,
        'g':5,
        'r':6,
        'i':7,
        'z':2
    }

    raw['wavelengths'] = raw['Filter'].replace(filter_wavelengths)
    raw['TelObs'] = raw['Tel'] + raw['Obs']

    unique_wavs = list(raw['wavelengths'].unique())

    LightcurveDictionary = {'name':[],
                      'wavelength':[],
                      'lightcurve':[]}
    for w in unique_wavs:
        rawwav = raw[raw['wavelengths'] == w]
        unique_telobs = list(rawwav['TelObs'].unique())
        for u in unique_telobs:
            lightcurve = rawwav[rawwav['TelObs']==u][['MJD','Flux','Error']].sort_values(by='MJD')
            LightcurveDictionary['name'].append(u)
            LightcurveDictionary['wavelength'].append(w)
            LightcurveDictionary['lightcurve'].append(lightcurve.values)




    # Prepare pycecream
    pc = pc.pycecream()
    for (name, wavelength, lcdat) in zip(LightcurveDictionary['name'],
                                         LightcurveDictionary['wavelength'],
                                         LightcurveDictionary['lightcurve']):

    pc.add_lc(lcdat,name=name,group=)
    #pc.add_lc()








