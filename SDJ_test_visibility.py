import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl





###### 1. Create a template MS for the model (DEFAULT MODEL)

# Instantiate a csalt model
cm = model('CSALT0')

cdir  = '/Users/sdjiang/.casa/data/alma/simmos/'
cm.template_MS('testdata/test.ms', 
               config=[cdir+'alma.cycle8.5.cfg', cdir+'alma.cycle8.8.cfg'],
               t_total=['2min', '7min'], t_integ='30s', observatory='ALMA',
               date=['2023/04/20', '2023/07/20'], HA_0=['0h', '1h'],
               restfreq=230.538e9, dnu_native=[122e3, 122e3],
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS('testdata/test.ms')


fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 161, 
            'restfreq': 230.538e9,
            'Nup': 5, 'doppcorr': 'approx'}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = -33      # inclination (degrees)
PA      = 90       # position angle (degrees)
mstar   = 1.1       # stellar mass (solar masses)

# Disk structure
r_l     = 300       # outer edge (au)
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
vlsr    = 4100      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])

pure_mdict  = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile='testdata/test_template_default.ms')

# Define some tclean keywords to be used in all imaging
tclean_kw = {'imsize': 1500, 'start': '2.35km/s', 'width': '0.35km/s',
             'nchan': 7, 'restfreq': '230.538GHz', 'cell': '0.01arcsec',
             'scales': [0, 10, 30, 60], 'niter': 50000,
             'robust': 1.0, 'threshold': '5mJy', 'uvtaper': '0.04arcsec'}

# Define some Keplerian mask keywords
kepmask_kw = {'inc': 33, 'PA': 150, 'mstar': 0.85, 'dist': 161, 'vlsr': 4100,
              'r_max': 1.1, 'nbeams': 1.5, 'zr': 0.2}

# Image the residual cubes
imagecube('testdata/test_template_default.ms', 'testdata/test_template_default',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)



# 2. Change disk parameters with template model 

# Instantiate a csalt model
cm = model('CSALT0')

cdir  = '/Users/sdjiang/.casa/data/alma/simmos/'
cm.template_MS('testdata/test.ms', 
               config=[cdir+'alma.cycle8.5.cfg', cdir+'alma.cycle8.8.cfg'],
               t_total=['2min', '7min'], t_integ='30s', observatory='ALMA',
               date=['2023/04/20', '2023/07/20'], HA_0=['0h', '1h'],
               restfreq=230.538e9, dnu_native=[122e3, 122e3],
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS('testdata/test.ms')


fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 161, 
            'restfreq': 230.538e9,
            'Nup': 5, 'doppcorr': 'approx'}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = -47.4       # inclination (degrees)
PA      = 162       # position angle (degrees)
mstar   = 0.1       # stellar mass (solar masses)

# Disk structure
r_l     = 300       # outer edge (au)
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
vlsr    = 4100      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])

pure_mdict  = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile='testdata/test_template_change.ms')

# Define some tclean keywords to be used in all imaging
tclean_kw = {'imsize': 1500, 'start': '2.35km/s', 'width': '0.35km/s',
             'nchan': 15, 'restfreq': '230.538GHz', 'cell': '0.01arcsec',
             'scales': [0, 10, 30, 60], 'niter': 50000,
             'robust': 1.0, 'threshold': '5mJy', 'uvtaper': '0.04arcsec'}

# Define some Keplerian mask keywords
kepmask_kw = {'inc': inc, 'PA': PA, 'mstar': mstar, 'dist': 161, 'vlsr': vlsr,
              'r_max': 2.1, 'nbeams': 1.5, 'zr': 0.2}

# Image the residual cubes
imagecube('testdata/test_template_change.ms', 'testdata/test_template_change',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)


# 3. Change disk parameters with SO1152 obs measurement set

# Instantiate a csalt model
cm = model('CSALT0')

# Get the data dictionary from the empty MS
ddict = read_MS('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg')


fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 391.3, 
            'restfreq': 345.7959899e9,
            'Nup': 5, 'doppcorr': 'approx'}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = 65.8      # inclination (degrees)
PA      = 144.9      # position angle (degrees)
mstar   = 0.58       # stellar mass (solar masses)

# Disk structure
r_l     = 300       # outer edge (au)
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

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])

pure_mdict  = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile='SigOriData/SO1152_template.ms')


##### DOUBLE CHECK THESE PARAMETERS BEFORE RUNNING THE SCRIPT ##### 
tclean_kw = {'imsize': 500, 'start': '7.0km/s', 'width': '0.5km/s',
             'nchan': 20, 'restfreq': '345.7959899GHz', 'cell': '0.02arcsec',
             'scales': [0, 10, 20, 50], 'niter': 50000,
             'robust': 0.5, 'threshold': '0.1mJy'}#

# Define some Keplerian mask keywords
kepmask_kw = {'inc': inc, 'PA': PA, 'mstar': mstar, 'dist': 391.3, 'vlsr': vlsr,
              'r_max': 1.2, 'nbeams': 1.0, 'zr': 0.0}

# Image the residual cubes
imagecube('SigOriData/SO1152_template.ms', 'SigOriData/SO1152_template_dirty',dirty_only=True,
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

imagecube('SigOriData/SO1152_template.ms', 'SigOriData/SO1152_template',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

imagecube('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg', 'SigOriData/SO1152_data',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)




### Show the results!

cubes = ['SigOriData/SO1152_template']
lbls = ['model']

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits', 
               velocity=True, overwrite=True) 

for i in range(len(cubes)):
    exportfits(cubes[i]+'.residual', fitsimage=cubes[i]+'_residual.fits', 
               velocity=True, overwrite=True) 