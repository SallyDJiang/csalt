import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl


# Instantiate a csalt model
cm = model('CSALT0')
# Define some fixed attributes for the model
# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 391.3, 'restfreq': 345.7959899e9,
            'Nup': 5, 'doppcorr': 'approx'}

# Get the data dictionary from the empty MS - just created 
ddict = read_MS('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg') #'SigOriData/SO1152_12CO.ms.contsub_lsrk_.ms') 
#ddict = read_MS('test.ms') #'SigOriData/SO1152_12CO.ms.contsub_lsrk_.ms') 

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = 65.8       # inclination (degrees)
PA      = 342       # position angle (degrees)
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
vlsr    = 12000      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])

pure_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)

write_MS(pure_mdict, outfile='SigOriData/SO1152_12CO_model.ms')

##### DOUBLE CHECK THESE PARAMETERS BEFORE RUNNING THE SCRIPT ##### 
tclean_kw = {'imsize': 600, 'start': '7.0km/s', 'width': '1.0km/s',
             'nchan': 10, 'restfreq': '345.7959899GHz', 'cell': '0.02arcsec',
             'scales': [0, 10, 20, 50], 'niter': 50000,
             'robust': 0.5, 'threshold': '3mJy'}#


# Define some Keplerian mask keywords
kepmask_kw = {'inc': 65.8, 'PA': 324.9, 'mstar': 0.58, 
              'dist': 391.3, 'vlsr': 12000,
              'r_max': 1.0, 'nbeams': 1.0, 'zr': 0.2}


# Image the residual cubes
imagecube('SigOriData/SO1152_12CO_model.ms', 'SigOriData/SO1152_12CO_model_image',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)


### Show the results!

cubes = ['SigOriData/SO1152_12CO_model_image']
lbls = ['pure']

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)
    






### tests 

# 1. CSALT VISIBILITY + SO1152 DISK PARAMETER
# Instantiate a csalt model
cm = model('CSALT0')

# Create an empty MS from scratch
#cdir = '/pool/asha0/casa-release-5.7.2-4.el7/data/alma/simmos/'
cdir  = '/Users/sdjiang/.casa/data/alma/simmos/' 
cm.template_MS('testdata/test.ms', 
               config=[cdir+'alma.cycle8.5.cfg', cdir+'alma.cycle8.8.cfg'],
               t_total=['2min', '7min'], t_integ='30s', observatory='ALMA',
               date=['2023/04/20', '2023/07/20'], HA_0=['0h', '1h'],
               restfreq=230.538e9, dnu_native=[122e3, 122e3],
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS('testdata/test.ms')

# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 391.3,
            'Nup': 5, 'doppcorr': 'approx'}

#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry
inc     = 65.8       # inclination (degrees)
PA      = 342       # position angle (degrees)
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
vlsr    = 4100      # systemic velocity (m/s)
dx      = 0         # RA offset (arcsec)
dy      = 0         # Dec offset (arcsec)

pars = np.array([inc, PA, mstar, r_l, z_10, z_q, Tb_10, Tb_q, Tmax_b, 
                 dV_10, logtau_10, tau_q, vlsr, dx, dy])

pure_mdict  = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile='testdata/test_PURE.ms')

# Define some tclean keywords to be used in all imaging
tclean_kw = {'imsize': 1500, 'start': '2.35km/s', 'width': '0.35km/s',
             'nchan': 7, 'restfreq': '230.538GHz', 'cell': '0.01arcsec',
             'scales': [0, 10, 30, 60], 'niter': 50000,
             'robust': 1.0, 'threshold': '2mJy', 'uvtaper': '0.04arcsec'}


# Define some Keplerian mask keywords
kepmask_kw = {'inc': 65.8, 'PA': 324.9, 'mstar': 0.58, 
              'dist': 391.3, 'vlsr': 4100,
              'r_max': 1.0, 'nbeams': 1.0, 'zr': 0.2}


# Image the residual cubes
imagecube('testdata/test_PURE.ms', 'testdata/test_PURE',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)






### test 2: my visibilities, old  model parameters

ddict = read_MS('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg') 


# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': 161,
            'Nup': 5, 'doppcorr': 'exact'}

# Set the CSALT model parameters
pars = np.array([-33, 90, 1.1, 120, 0.3, 1.5, 120, -0.5, 20., 217,
                 2.2, -1, 4100, 0, 0])

# Calculate a model dictionary; insert it to model MS files
pure_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile='testdata/test_PURE.ms')

# Define some tclean keywords to be used in all imaging
tclean_kw = {'imsize': 1500, 'start': '7.0km/s', 'width': '0.5km/s',
             'nchan': 23,  'restfreq': '345.7959899GHz', 'cell': '0.01arcsec',
             'scales': [0, 10, 30, 60], 'niter': 50000,
             'robust': 1.0, 'threshold': '5mJy', 'uvtaper': '0.04arcsec'}

# Define some Keplerian mask keywords
kepmask_kw = {'inc': 33, 'PA': 150, 'mstar': 0.85, 'dist': 161, 'vlsr': 4100,
              'r_max': 1.1, 'nbeams': 1.5, 'zr': 0.2}


imagecube('testdata/test_PURE.ms', 'testdata/test_PURE',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

cubes = ['testdata/test_PURE', 'testdata/test_NOISY']
lbls = ['pure', 'noisy']

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)












##### DOUBLE CHECK THESE PARAMETERS BEFORE RUNNING THE SCRIPT ##### 
tclean_kw = {'imsize': 600, 'start': '7.0km/s', 'width': '1.0km/s',
             'nchan': 10, 'restfreq': '345.7959899GHz', 'cell': '0.02arcsec',
             'scales': [0, 10, 20, 50], 'niter': 50000,
             'robust': 0.5, 'threshold': '3mJy'}#


# Define some Keplerian mask keywords
kepmask_kw = {'inc': 65.8, 'PA': 324.9, 'mstar': 0.58, 
              'dist': 391.3, 'vlsr': 12000,
              'r_max': 1.0, 'nbeams': 1.0, 'zr': 0.2}



## TEST: clean SO1152_image: 
# Image the residual cubes
imagecube('SigOriData/SO1152_12CO.ms.contsub_lsrk_avg', 'SigOriData/SO1152_12CO_contsub_lsrk_avg',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)



### interative cleaning of the avg 12CO data

################################### make clean images 
combined_obs_name_12CO = 'SO1152_12CO_model.ms'
image_name = 'SO1152_12CO_model_image_cleaned'
os.system(f'rm -rf {image_name}*')
tclean(vis=combined_obs_name_12CO, #replace with whatever your measurement set is called
           imagename=image_name, #replace with whatever you're calling your image
           specmode='cube',
           restfreq='345.7959899GHz', ##### CHANGE
           deconvolver = 'multiscale',
           scales = [0, 10, 20, 50],
           weighting='briggs',
           niter = 50000,
           robust=0.5,
           cell='0.01arcsec', #replace with whatever cellsize you're using
           imsize=500, #replace with whatever imsize you're using
           restoringbeam='common',  
		   threshold='3mJy', #StdDev * 2-3 ### CHANGE ###
           interactive=True,
           start='7km/s',  
           width='0.5km/s', 
           nchan=23)
           mask=prefix+"_12CO.ms.contsub_image.mask.image")


################################### make clean images 
combined_obs_name_12CO = 'SO1152_12CO.ms.contsub_lsrk_avg'
image_name = 'SO1152_12CO_contsub_lsrk_avg'
os.system(f'rm -rf {image_name}*')
tclean(vis=combined_obs_name_12CO, #replace with whatever your measurement set is called
           imagename=image_name, #replace with whatever you're calling your image
           specmode='cube',
           restfreq='345.7959899GHz', ##### CHANGE
           deconvolver = 'multiscale',
           scales = [0, 10, 20, 50],
           weighting='briggs',
           niter = 50000,
           robust=0.5,
           cell='0.01arcsec', #replace with whatever cellsize you're using
           imsize=500, #replace with whatever imsize you're using
           restoringbeam='common',  
		   threshold='3mJy', #StdDev * 2-3 ### CHANGE ###
           interactive=True,
           start='7km/s',  
           width='0.5km/s', 
           nchan=23)