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
obs_measurement_set = f'{datapath}SO1152_12CO.ms.contsub_lsrk_avg'  # path to the measurement set
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
rest_freq_Hz  = 345.7959899e9  # rest frequency of the line (12CO 3-2)

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
doppcorr = 'approx'  # Doppler correction method ("exact" is faster)
noise = None


##### COMPARE MODELS ###### 

tclean_kw = {'imsize': imsize, 'start': start_vel, 'width': width,
             'nchan': nchan, 'restfreq': rest_freq, 'cell': cell_size,
             'scales': scales, 'niter': niter,
             'robust': robust_val, 'threshold': thres}#

# Define some Keplerian mask keywords
kepmask_kw = {'inc': inc, 'PA': PA, 'mstar': mstar, 'dist': dist, 'vlsr': vlsr,
              'r_max': r_max, 'nbeams': nbeam, 'zr': zr}


# Define some fixed attributes for the model
if noise is None:
    fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr}
else: 
    fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr, 'noise_inject': noise}


# Read in the data MS
ddict = read_MS(obs_measurement_set) 

# Instantiate a csalt model
cm = model('CSALT0')


### change r of outer edge 


for r_var in np.array([50, 120, 150, 200]):  # outer edge (au)
    pars = np.array([inc, PA, mstar, r_var, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                    dV_10, logtau_10, tau_q, vlsr, dx, dy])
    mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
    write_MS(mdict, outfile=f'{mdict_name_prefix}_model_R{r_var}.ms')

for r_var in np.array([50, 120, 150, 200]):  # outer edge (au)
    imagecube(f'{mdict_name_prefix}_model_R{r_var}.ms', f'{mdict_name_prefix}_MODEL_R{r_var}',
            kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

cubes = [f'{mdict_name_prefix}_MODEL_R{r_var}' for r_var in np.array([50, 120, 150, 200])]
lbls = np.array([50, 120, 150, 200]).astype(str)

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)
    


### change mass of the star and see how it affects the model

# Set the CSALT model multi parameters # 
for mstar_var in np.array([0.1, 0.2, 0.5, 1.0]):  
    pars = np.array([inc, PA, mstar_var, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                    dV_10, logtau_10, tau_q, vlsr, dx, dy])
    mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
    write_MS(mdict, outfile=f'{mdict_name_prefix}_model_M{mstar_var}.ms')

for mstar_var in np.array([0.1, 0.2, 0.5, 1.0]):  
    imagecube(f'{mdict_name_prefix}_model_M{mstar_var}.ms', f'{mdict_name_prefix}_MODEL_M{mstar_var}',
            kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

cubes = [f'{mdict_name_prefix}_MODEL_M{mstar_var}' for mstar_var in np.array([0.1, 0.2, 0.5, 1.0])]
lbls = np.array([0.1, 0.2, 0.5, 1.0]).astype(str)

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)
    


### change parameter sweep function to make it easier to run multiple models with varying parameters


import numpy as np
from casatasks import exportfits

# Order must match what cm.modeldict expects
PARAM_ORDER = ['inc', 'PA', 'mstar', 'r_l', 'z_10', 'z_q', 'Tb_10', 'Tb_q',
               'Tmax_b', 'dV_10', 'logtau_10', 'tau_q', 'vlsr', 'dx', 'dy']

# Baseline/default values for every parameter (edit to match your fiducial model)
base_params = {
    'inc': inc, 'PA': PA, 'mstar': mstar, 'r_l': r_l, 'z_10': z_10, 'z_q': z_q,
    'Tb_10': Tb_10, 'Tb_q': Tb_q, 'Tmax_b': Tmax_b, 'dV_10': dV_10,
    'logtau_10': logtau_10, 'tau_q': tau_q, 'vlsr': vlsr, 'dx': dx, 'dy': dy,
}

def run_param_sweep(param_name, values, tag, base_params=base_params,
                     ddict=ddict, fixed_kw=fixed_kw, kepmask_kw=kepmask_kw,
                     tclean_kw=tclean_kw, mdict_name_prefix=mdict_name_prefix):

    values = np.asarray(values)
    labels = values.astype(str)
    cubes = []

    # --- 1. Write model MS files ---
    for val in values:
        pars_dict = dict(base_params)
        pars_dict[param_name] = val
        pars = np.array([pars_dict[k] for k in PARAM_ORDER])

        mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
        outfile = f'{mdict_name_prefix}_model_{tag}{val}.ms'
        write_MS(mdict, outfile=outfile)

    # --- 2. Image each MS ---
    for val in values:
        msfile = f'{mdict_name_prefix}_model_{tag}{val}.ms'
        cube_name = f'{mdict_name_prefix}_MODEL_{tag}{val}'
        imagecube(msfile, cube_name, kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)
        cubes.append(cube_name)

    # --- 3. Export to FITS ---
    for cube in cubes:
        exportfits(cube + '.image', fitsimage=cube + '.fits',
                   velocity=True, overwrite=True)

    return cubes, labels


r_cubes, r_lbls = run_param_sweep('r_l', [50, 120, 150, 200], tag='R')
m_cubes, m_lbls = run_param_sweep('mstar', [0.1, 0.2, 0.5, 1.0], tag='M')
zq_cubes, zq_lbls = run_param_sweep('z_q', [0.5, 1.0, 1.5, 2.0], tag='Zq')



#### TERMINAL 

# Set these for the parameter sweep
param_name = 'z_q'
values = np.asarray([0.5, 1.0, 1.5, 2.0])
tag = 'Zq'

labels = values.astype(str)
cubes = []

# --- 1. Write model MS files ---
for val in values:
    pars_dict = dict(base_params)
    pars_dict[param_name] = val
    pars = np.array([pars_dict[k] for k in PARAM_ORDER])
    mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
    outfile = f'{mdict_name_prefix}_model_{tag}{val}.ms'
    write_MS(mdict, outfile=outfile)

# --- 2. Image each MS --- (this takes a while)...
for val in values:
    msfile = f'{mdict_name_prefix}_model_{tag}{val}.ms'
    cube_name = f'{mdict_name_prefix}_MODEL_{tag}{val}'
    imagecube(msfile, cube_name,
              kepmask_kwargs=kepmask_kw,
              tclean_kwargs=tclean_kw)
    cubes.append(cube_name)

# --- 3. Export to FITS ---
for cube in cubes:
    exportfits(cube + '.image',
               fitsimage=cube + '.fits',
               velocity=True,
               overwrite=True)

# Save the results
result_cubes = cubes
result_labels = labels