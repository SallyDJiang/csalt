import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl

##### DOUBLE CHECK THESE PARAMETERS BEFORE RUNNING THE SCRIPT ##### 
tclean_kw = {'imsize': 600, 'start': '7.0km/s', 'width': '1.0km/s',
             'nchan': 10, 'restfreq': '345.7959899GHz', 'cell': '0.02arcsec',
             'scales': [0, 10, 20, 50], 'niter': 50000,
             'robust': 0.5, 'threshold': '3mJy'}#


# Define some Keplerian mask keywords
kepmask_kw = {'inc': 65.8, 'PA': 324.9, 'mstar': 0.58, 'dist': 391.3, 'vlsr': 12000,
              'r_max': 1.0, 'nbeams': 1.0, 'zr': 0.2}


# Image the residual cubes
imagecube('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg', 'SigOriData/SO1152_12CO_PRE',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)


##### MODEL #####

# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 391.3,
            'Nup': 5, 'doppcorr': 'approx', 'noise_inject': 0.005}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = -47.4       # inclination (degrees)
PA      = 162       # position angle (degrees)
mstar   = 0.1       # stellar mass (solar masses)

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
vlsr    = 11000      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])


# Instantiate a csalt model
cm = model('CSALT0')

cdir  = '/Users/sdjiang/.casa/data/alma/simmos/'

cm.template_MS('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg', 
               config = [cdir+'alma.cycle12.6.cfg'],
               t_total=['37.88min'], # time on source (OBS ID and FIELD ID)
               t_integ='6.0s', # 
               observatory='ALMA',
               date=['2025/11/13'], # date of observation
               HA_0=['0.492h'], #  Hour angle at start EB (tstart - ra_angle) 
               V_span=5800, ### added to specify the velocity span (m/s)
               restfreq=345.7959899e9, # 12CO 3-2
               dnu_native=[576725.6325073242], # ChanWid 
               RA='05:39:39.381821', DEC='-02:17:04.50121')



# Get the data dictionary from the empty MS - just created 
ddict = read_MS('SigOriData/SO1152_12CO.ms.contsub_lsrk_.ms') #'SigOriData/SO1152_12CO.ms.contsub_lsrk_.ms') 

pure_mdict, _ = cm.modeldict(ddict, pars, kwargs=fixed_kw)

write_MS(pure_mdict, outfile='SigOriData/SO1152_12CO_PURE.ms')

# Image the residual cubes
imagecube('SigOriData/SO1152_12CO_PURE.ms', 'SigOriData/SO1152_12CO_PURE',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

### Show the results!

cubes = ['SigOriData/SO1152_12CO_PURE']
lbls = ['pure']

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)
      