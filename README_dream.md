# PyceCREAM


Reverberation Mapping requires telescope images to be converted to light curves. These images are often taken 
from multiple telescopes over several months with unique calibration anomalies between the telescopes. 
These calibration anomalies are often visible in the final combined light curve. 

##Example Figure showing lighcurves with calibration problem

Often the individual light curves are rescaled to a reference mean and standard deviation 
prior to merging.

Here Pycecream is used to combine light curves in a more sophisticated way by rescaling each light curve
to pycecreams random walk model fit. In addition to correcting the calibration fluxes, 
dream also modifies the input error bars with a multiplicative (f) and additive (V) parameter.


\sigma^{2} = \left ( f \sigma_0 \right )^{2} + V


## Installation (requires Python 3)

`pip install pycecream`


## g-band merging example

The example below shows how to use dream to merge example g-band light curves 
from 5 telescopes and access the merged output. 

```python
import pycecream as pc
import pickle

#initialise dream instance
dream = pc.dream(Niterations = 200)

#add each light curve ('dat' should be a N x 3 array of time, flux, errorbar)
#errorbar_variance, errorbar_rescale should be True to optimise the 'f' and 'V' error bar parameters 
dream.add_lc(dat1, 'g-band 1', errorbar_variance=True, errorbar_rescale=True)
dream.add_lc(dat2, 'g-band 2', errorbar_variance=True, errorbar_rescale=True)
dream.add_lc(dat3, 'g-band 3', errorbar_variance=True, errorbar_rescale=True)
dream.add_lc(dat4, 'g-band 4', errorbar_variance=True, errorbar_rescale=True)
dream.add_lc(dat5, 'g-band 5', errorbar_variance=True, errorbar_rescale=True)



#run the simulation
dream.run()


#access the input lightcurves
input = dream.lcinput

#access the combined merged light curve
merged_combined = dream.lc_combined

#access the individual (but rescaled light curves)
merged_individual = dream.lc_merged_individual



#OPTIONAL: Save the output for later
os.system('rm ' + picklefile)
pickle_out = open(picklefile, "wb")
pickle.dump(dream, pickle_out)
pickle_out.close()

```


