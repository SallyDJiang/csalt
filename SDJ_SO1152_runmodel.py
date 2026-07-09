import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl


### IMPORTANT THINGS TO CHANGE 

datapath = "SigOriData/"
obs_measurement_set = f'{datapath}SO1152_12CO.ms.contsub_lsrk_.ms'  # path to the measurement set
mdict_name_prefix = f"{datapath}SO1152_12CO"

########----- MODEL PARAMETERS BELOW (CHANGE) -----#########
# Disk geometry
inc     = 65.8      # inclination (degrees)
PA      = 324.9      # position angle (degrees)
mstar   = 0.58       # stellar mass (solar masses)
dist    = 391.3     # distance (pc)

# Disk structure
r_l     = 120       # outer edge (au)
z_10    = 0.3       # emission surface height normalization (" at r = 1")
z_q     = 1.5       # emission surface height gradient

# Temperature structure
Tb_10   = 120       # brightness temperature at 10 au (K)
Tb_q    = -0.5      # brightness temperature gradient
Tmax_b  = 20.       # max brightness temperature (K)

# Velocity structure
dV_10   = 217       # linewidth at 10 au (m/s)

# Optical depth
logtau_10 = 2.2     # log optical depth at 10 au
tau_q     = -1      # optical depth gradient

# Kinematics
vlsr    = 12000      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

# Tclean parameters
imsize     = 500      # image size (pixels)
start_vel  = '7.0km/s'       # starting velocity for imaging
width      = '1.0km/s'    # velocity width for imaging
nchan      = 10      # number of channels to image
rest_freq  = '345.7959899GHz'  # rest frequency of the line (12CO 3-2)
cell_size  = '0.02arcsec'  # cell size for imaging
scales     =  [0, 10, 20, 50]
niter      = 50000
robust_val =  0.5
thres      = '3mJy'  

#  Keplerian mask parameters
r_max  = 1.0       # maximum radius for the mask (arcsec)
nbeam  = 1.0       # number of beams for the mask
zr    = 0.2       # height of the mask (arcsec)

# fixed model info (not to be changed, unless needed)
FOV_val = 5.11      # field of view (arcsec)
Npix_val = 512       # number of pixels in the model
Nup_val = 5         
doppcorr = 'exact'  # Doppler correction method ("exact" is faster)
noise = None


##### COMPARE MODELS ###### 


# Read in the data MS
ddict = read_MS(obs_measurement_set) 

# Instantiate a csalt model
cm = model('CSALT')


tclean_kw = {'imsize': imsize, 'start': start_vel, 'width': width,
             'nchan': nchan, 'restfreq': rest_freq, 'cell': cell_size,
             'scales': scales, 'niter': niter,
             'robust': robust_val, 'threshold': thres}#

# Define some Keplerian mask keywords
kepmask_kw = {'inc': inc, 'PA': PA, 'mstar': mstar, 'dist': dist, 'vlsr': vlsr,
              'r_max': r_max, 'nbeams': nbeam, 'zr': zr}


# Set the CSALT model A parameters # 
inc = 65.8      # inclination (degrees)
pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])
mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(mdict, outfile=f'{mdict_name_prefix}_modelA.ms')

# Set the CSALT model B parameters # 
inc = 160      # inclination (degrees)
pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])


mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(mdict, outfile=f'{mdict_name_prefix}_modelB.ms')



### image both models 

# Image the data, model, and residual cubes
imagecube('testdata/Sz129_CO.ms', 'testdata/Sz129_DATA', 
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

## reuse mask (from modelA)
tclean_kw['mask'] = 'testdata/Sz129_DATA.mask'
imagecube('testdata/Sz129_MODEL.ms', 'testdata/Sz129_MODEL', 
          mk_kepmask=False, tclean_kwargs=tclean_kw)