# PyceCREAM


Reverberation Mapping requires telescope images to be converted to light curves. These images are often taken 
from multiple telescopes over several months with unique calibration anomalies between the telescopes. 
These calibration anomalies are often visible in the final combined light curve. 

##Example Figure showing lighcurves with calibration problem

Often the individual light curves are rescaled to a reference mean and standard deviation 
prior to merging. WHATS WRONG WITH THIS?

Here Pycecream is used to combine light curves in a more sophisticated way by rescaling each light curve
to pycecreams random walk model fit. In addition to correcting the calibration fluxes, 
dream also modifies the input error bars with a multiplicative (f) and additive (V) parameter.

$\sigma$

\sigma^{2} = \left ( f \sigma_0 \right )^{2} + V


## g-band merging example

The example below shows how to use dream to merge example g-band light curves 
from 5 telescopes and access the merged output. 

```python

import pycecream.pydream as dream

```


