import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl
import json


#datapath          = 'SigOriData/SO844/'
#source_name      = "SO844"
#priors = importlib.import_module('priors_CSALT0_SDJ_SO844_B')
cm = model('CSALT0_SDJ_SO844_B')

test_Hz = 345.7959899e9

cdir  = '/Users/sdjiang/.casa/data/alma/simmos/'

cm.template_MS('testdata/test.ms', 
               config=[cdir+'alma.cycle8.5.cfg', cdir+'alma.cycle8.8.cfg'],
               t_total=['2min', '7min'], t_integ='30s', observatory='ALMA',
               date=['2023/04/20', '2023/07/20'], HA_0=['0h', '1h'],
               restfreq=test_Hz, dnu_native=[122e3, 122e3],
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS('testdata/test.ms') 

# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 161,
            'Nup': 5, 'doppcorr': 'approx'}#, 'noise_inject': 0.005}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = 10.0   # inclination (degrees)
PA      = 162       # position angle (degrees)
mstar   = 0.1       # stellar mass (solar masses)

# Disk structure
Rout     = 120       # outer edge (au)
z_10    = 0.3       # emission surface height normalization (" at r = 1")
z_q     = 1.5       # emission surface height gradient

# Temperature structure
Tb_10   = 120       # brightness temperature at 10 au (K)
Tb_q    = -0.5      # brightness temperature gradient
Tmax_b  = 20.       # max brightness temperature (K)

# Velocity structure
dV_10   = 217       # linewidth at 10 au (m/s)

# Optical depth
#logtau_10 = 2.2     # log optical depth at 10 au
tau_0 = 15
tau_q     = -1      # optical depth gradient

# Kinematics
vlsr    = 12000      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, Rout, Tb_10, Tb_q, 
                tau_0, vlsr, dx, dy])

pure_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)

write_MS(pure_mdict, outfile='test_.ms')

# Define velocity axis — same way modelset does it internally
# (you need to pick a velocity range covering your emission)
vlsr   = 12000                         # systemic velocity [m/s]
dv     = 350.                              # channel width [m/s]
nchan  = 7                               # number of channels
velax  = vlsr + dv * (np.arange(nchan) - nchan // 2)  # [m/s]

# Call cube() directly — returns SkyImage before vis_sample touches it
icube = cm.cube(velax, pars=pars,
                restfreq=test_Hz,
                FOV=5.11,
                Npix=512,
                dist=161)

# icube is a vis_sample SkyImage with these attributes:
print(type(icube))
print('data shape:    ', icube.data.shape)   # (npix, npix, nchan)
print('RA axis:       ', icube.ra.shape)     # arcsec offsets
print('Dec axis:      ', icube.dec.shape)
print('freq axis:     ', icube.freqs.shape)  # Hz


# Write it to FITS yourself
cube_to_fits(icube, 'demodisk_SkyImage.fits',
             restfreq=test_Hz)





pure_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)

write_MS(pure_mdict, outfile='test_.ms')

# Define velocity axis — same way modelset does it internally
# (you need to pick a velocity range covering your emission)
vlsr   = 12000                         # systemic velocity [m/s]
dv     = 350.                              # channel width [m/s]
nchan  = 7                               # number of channels
velax  = vlsr + dv * (np.arange(nchan) - nchan // 2)  # [m/s]

# Call cube() directly — returns SkyImage before vis_sample touches it
icube = cm.cube(velax, pars=pars,
                restfreq=test_Hz,
                FOV=5.11,
                Npix=512,
                dist=161)

# icube is a vis_sample SkyImage with these attributes:
print(type(icube))
print('data shape:    ', icube.data.shape)   # (npix, npix, nchan)
print('RA axis:       ', icube.ra.shape)     # arcsec offsets
print('Dec axis:      ', icube.dec.shape)
print('freq axis:     ', icube.freqs.shape)  # Hz


# Write it to FITS yourself
cube_to_fits(icube, 'demodisk_SkyImage.fits',
             restfreq=test_Hz)