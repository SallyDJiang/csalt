import os
import sys
import multiprocessing #
multiprocessing.set_start_method('fork') #
import numpy as np
import casatools ## need to install first 
from csalt.model import *
from csalt.helpers import *

#export DYLD_LIBRARY_PATH=/Users/sdjiang/.venvs/casa_env/lib/python3.11/site-packages/casatools/__casac__/lib
#pip install h5py tqdm


###################################CHANGEABLE PARAMETERS############################################
####################################################################################################
# Read in the data MS
datapath = "SigOriData/"
obs_measurement_set = f'{datapath}SO1152_12CO.ms.contsub_lsrk_avg'  # path to the measurement set
mdict_name_prefix = f"{datapath}SO1152_12CO_run1"

# Instantiate a csalt model
cm = model('CSALT0_SDJ')

ddict = read_MS(obs_measurement_set) 

# fixed model info (not to be changed, unless needed)
dist    = 391.3     # distance (pc)
rest_freq  = '345.7959899GHz'  # rest frequency of the line (12CO 3-2)
rest_freq_Hz  = 345.7959899e9  # rest frequency of the line (12CO 3-2)
FOV_val = 5.11      # field of view (arcsec)
Npix_val = 512       # number of pixels in the model
Nup_val = 5         
doppcorr = 'approx'  # Doppler correction method ("exact" is faster)
noise = None

fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr, 'online_avg': 4}

####################################################################################################

####################################################################################################
# Read in the data MS
datapath = "SigOriData/SO1153/"
obs_measurement_set = f'{datapath}SO1153_12CO.ms.contsub_lsrk_avg'  # path to the measurement set
mdict_name_prefix = f"{datapath}SO1153_12CO_run1"

# Instantiate a csalt model
cm = model('CSALT0_SDJ_SO1153')

ddict = read_MS(obs_measurement_set) 

# fixed model info (not to be changed, unless needed)
dist    = 390.4     # distance (pc)
rest_freq  = '345.7959899GHz'  # rest frequency of the line (12CO 3-2)
rest_freq_Hz  = 345.7959899e9  # rest frequency of the line (12CO 3-2)
FOV_val = 5.11      # field of view (arcsec)
Npix_val = 512       # number of pixels in the model
Nup_val = 5         
doppcorr = 'approx'  # Doppler correction method ("exact" is faster)
noise = None

fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr, 'online_avg': 4}

_ = cm.sample_posteriors(obs_measurement_set, kwargs=fixed_kw,
                         vra=[7e3, 1.8e4], restfreq=rest_freq_Hz, 
                         Nwalk=75, Nthreads=12, Ninits=20, Nsteps=36000,
                         outpost=mdict_name_prefix+'_posteriors.h5')

####################################################################################################












# Sample the posteriors!
_ = cm.sample_posteriors(obs_measurement_set, kwargs=fixed_kw,
                         vra=[7e3, 1.8e4], restfreq=rest_freq_Hz, 
                         Nwalk=75, Nthreads=12, Ninits=20, Nsteps=6000,
                         outpost=mdict_name_prefix+'_posteriors.h5')

#EB 0 SCOV inverted with direct calculation.
#100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [00:00<00:00, 99.15it/s]
#100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 50/50 [00:00<00:00, 84.77it/s]
#backend run in  0.6036109924316406

#   This run took 0.00 hours

### Nwalkers -> Nwalk

_ = cm.sample_posteriors(obs_measurement_set, kwargs=fixed_kw,
                         vra=[7e3, 1.8e4], restfreq=rest_freq_Hz, 
                         Nwalk=75, Nthreads=10, Ninits=20, Nsteps=6000,
                         outpost=mdict_name_prefix+'_posteriors.h5',
                         append=True)