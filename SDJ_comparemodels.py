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


## Best Model Parameters from MCMC run (from SDJ_infer.py)

########----- MODEL PARAMETERS BELOW (CHANGE) -----#########
# Disk geometry
inc     = 66.405    # inclination (degrees)          +0.185 -0.280
PA      = 324.434    # position angle (degrees)       +0.179 -0.172
mstar   = 0.534      # stellar mass (solar masses)    +0.005 -0.004
dist    = 391.3      # distance (pc)

# Disk structure
r_l     = 120        # outer edge (au)
z_10    = 0.157      # emission surface height normalization (" at r = 1")   +0.003 -0.004
z_q     = 0.691       # emission surface height gradient                     +0.012 -0.015
z_tap   = 111.732     # emission surface height taper radius (au)            +9.701 -10.115
z_exp   = 4.945       # emission surface height taper exponent               +0.042 -0.155

# Temperature structure
Tb_10   = 56.888     # brightness temperature at 10 au (K)   +1.906 -1.533
Tb_q    = -0.188      # brightness temperature gradient       +0.011 -0.012
Tmax_b  = 20.         # max brightness temperature (K)

# Velocity structure
dV_10   = 217         # linewidth at 10 au (m/s)

# Optical depth
logtau_10 = 4.870     # log optical depth at 10 au        +0.081 -0.272
tau_q     = 0.937      # optical depth gradient             +0.047 -0.151
tau_tap   = 148.639    # optical depth taper radius (au)    +3.958 -4.201
tau_exp   = 9.533      # optical depth taper exponent       +1.243 -1.221

# Kinematics
vlsr    = 12000   # systemic velocity (m/s)   +5.304 -5.031
dx      = -0.017      # RA offset (arcsec)         +0.001 -0.001
dy      = 0.007       # Dec offset (arcsec)        +0.001 -0.001

theta = np.array([66.458, 324.447, 0.533, 0.158, 
                  0.693, 56.513, -0.186, 110.083, 4.965, 4.909, 0.957, 
                  149.775, 9.984, 12162.681,  -0.016, 0.007])

# Tclean parameters
imsize     = 500      # image size (pixels)
start_vel  = '7.0km/s'       # starting velocity for imaging
width      = '0.5km/s'    # velocity width for imaging
nchan      = 23      # number of channels to image
rest_freq  = '345.7959899GHz'  # rest frequency of the line (12CO 3-2)
rest_freq_Hz  = 345.7959899e9  # rest frequency of the line (12CO 3-2)

cell_size  = '0.02arcsec'  # cell size for imaging
scales     =  [0, 10, 20, 50]
niter      = 50000
robust_val =  0.5
thres      = '1mJy'  

#  Keplerian mask parameters
r_max  = 2.0       # maximum radius for the mask (arcsec)
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
            'Nup': Nup_val, 'doppcorr': doppcorr, 'online_avg': 4}
else: 
    fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr, 'noise_inject': noise, 'online_avg': 4}


### Mak



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
            'Nup': Nup_val, 'doppcorr': doppcorr, 'online_avg': 4}
else: 
    fixed_kw = {'FOV': FOV_val, 'Npix': Npix_val, 'dist': dist, 'restfreq': rest_freq_Hz,
            'Nup': Nup_val, 'doppcorr': doppcorr, 'noise_inject': noise, 'online_avg': 4}


# Read in the data MS
imagecube(f'{obs_measurement_set}', f'{obs_measurement_set}',
        kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)





# Read in the data MS
ddict = read_MS(obs_measurement_set) 

# Instantiate a csalt model
cm = model('CSALT0_SDJ')

pars = theta
mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(mdict, outfile=f'{mdict_name_prefix}_Bestmodel.ms')
write_MS(mdict, outfile=f'{mdict_name_prefix}_Bestmodel_RESD.ms', resid=True)

imagecube(f'{mdict_name_prefix}_Bestmodel.ms', f'{mdict_name_prefix}_Bestmodel',
        kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

tclean_kw['mask'] = f'{mdict_name_prefix}_Bestmodel.mask'
imagecube(f'{mdict_name_prefix}_Bestmodel.ms', f'{mdict_name_prefix}_Bestmodel', 
          mk_kepmask=False, tclean_kwargs=tclean_kw)

imagecube(f'{mdict_name_prefix}_Bestmodel_RESD.ms', f'{mdict_name_prefix}_Bestmodel_RESD', 
          mk_kepmask=False, tclean_kwargs=tclean_kw)


cubes = [f'{obs_measurement_set}', f'{mdict_name_prefix}_Bestmodel', f'{mdict_name_prefix}_Bestmodel_RESD']
lbls = ['data', 'model', 'residual']


# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)


# Plot a subset of the cube as a direct comparison
import importlib
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
from matplotlib.patches import Ellipse
from astropy.io import fits
import scipy.constants as sc
import cmasher as cmr

nchans = 13

fig, axs = plt.subplots(nrows=3, ncols=nchans, figsize=(7.5, 3.2))
fl, fr, fb, ft, hs, ws = 0.06, 0.88, 0.13, 0.99, 0.12, 0.03
xlims = [1.1, -1.1]
ylims = [-1.1, 1.1]

for ia in range(len(cubes)):

    # Load the cube and header data
    hdu = fits.open(cubes[ia]+'.fits')
    Ico, h = np.squeeze(hdu[0].data), hdu[0].header
    hdu.close()

    # coordinate grids, beam parameters
    dx = 3600 * h['CDELT1'] * (np.arange(h['NAXIS1']) - (h['CRPIX1'] - 1))
    dy = 3600 * h['CDELT2'] * (np.arange(h['NAXIS2']) - (h['CRPIX2'] - 1))
    ext = (dx.max(), dx.min(), dy.min(), dy.max())
    bmaj, bmin, bpa = h['BMAJ'], h['BMIN'], h['BPA']
    bm = (np.pi * bmaj * bmin / (4 * np.log(2))) * (np.pi / 180)**2

    # spectral information
    vv = h['CRVAL3'] + h['CDELT3'] * (np.arange(h['NAXIS3']) - (h['CRPIX3']-1))
    ff = 230.538e9 * (1 - vv / sc.c)
    j0 = np.argmin(np.abs(vv - 9.0e3))   # vv is in m/s

    for i in range(nchans):
        
        # in-cube index
        j = i + j0 ## begin at 9 km/s

        # convert intensities to brightness temperatures
        Tb = (1e-26 * np.squeeze(Ico[j,:,:]) / bm) * sc.c**2 / \
             (2 * sc.k * ff[j]**2)

        # allocate the panel
        ax = axs[ia, i]
        pax = ax.get_position()

        # plot the channel map
        im = ax.imshow(Tb, origin='lower', cmap='cmr.cosmic', extent=ext,
                       aspect='auto', vmin=0, vmax=25)

        # set map boundaries
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        # overlay beam dimensions
        beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                        -xlims[0] - 0.1*np.diff(xlims)), 
                       3600 * bmaj, 3600 * bmin, angle=90-bpa)
        beam.set_facecolor('w')
        ax.add_artist(beam)

        # labeling
        if i == 0:
            ax.text(0.02, 0.90, lbls[ia], transform=ax.transAxes, ha='left',
                    va='center', style='italic', fontsize=8, color='w')
        if ia == 2:
            if np.abs(vv[j]) < 0.001: vv[j] = 0.0
            if np.logical_or(np.sign(vv[j]) == 1, np.sign(vv[j]) == 0):
                pref = '+'
            else:
                pref = ''
            vstr = pref+'%.2f' % (1e-3 * vv[j])
            ax.text(0.97, 0.08, vstr, transform=ax.transAxes, ha='right',
                    va='center', fontsize=7, color='w')
        if np.logical_and(ia == 2, i == 0):
            ax.set_xlabel('RA offset  ($^{\prime\prime}$)')
            ax.set_ylabel('DEC offset  ($^{\prime\prime}$)')
            ax.spines[['right', 'top']].set_visible(False)
        else:
            ax.axis('off')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])

# colorbar
cbax = fig.add_axes([fr+0.01, fb+0.01, 0.02, ft-fb-0.02])
cb = Colorbar(ax=cbax, mappable=im, orientation='vertical',
              ticklocation='right')
cb.set_label('$T_{\\rm b}$  (K)', rotation=270, labelpad=13)

fig.subplots_adjust(left=fl, right=fr, bottom=fb, top=ft, hspace=hs, wspace=ws)
fig.savefig('testdata/demo_compare.pdf')
fig.clf()


nchans = 13

fig, axs = plt.subplots(nrows=3, ncols=nchans, figsize=(7.5*nchans/10, 3.2))
fl, fr, fb, ft, hs, ws = 0.06, 0.88, 0.13, 0.99, 0.12, 0.03
xlims = [1.1, -1.1]
ylims = [-1.1, 1.1]

# --- first pass: load all cubes and compute robust color scales ---
# --- first pass: load all cubes and compute robust color scales ---
Tb_all, vv_list, j0_list = [], [], []
for ia in range(len(cubes)):
    hdu = fits.open(cubes[ia]+'.fits')
    Ico, h = np.squeeze(hdu[0].data), hdu[0].header
    hdu.close()

    bmaj, bmin = h['BMAJ'], h['BMIN']
    bm = (np.pi * bmaj * bmin / (4 * np.log(2))) * (np.pi / 180)**2

    vv = h['CRVAL3'] + h['CDELT3'] * (np.arange(h['NAXIS3']) - (h['CRPIX3']-1))
    ff = 230.538e9 * (1 - vv / sc.c)

    Tb_cube = (1e-26 * Ico / bm) * sc.c**2 / (2 * sc.k * ff[:, None, None]**2)
    Tb_all.append(Tb_cube)
    vv_list.append(vv)
    j0_list.append(np.argmin(np.abs(vv - 9.0e3)))

# only look at the channels actually being plotted (9-15 km/s), not the full cube
data_stack = Tb_all[0][j0_list[0]:j0_list[0]+nchans]
model_stack = Tb_all[1][j0_list[1]:j0_list[1]+nchans]

vmax_data = np.nanmax([np.nanpercentile(data_stack, 99.9),
                        np.nanpercentile(model_stack, 99.9)])
vmax_data *= 2.15   # headroom so the true peak isn't clipped at the very top

print(f"vmax_data = {vmax_data:.2f} K")   # sanity check against your old vmax=25

resid = Tb_all[2]
rms = np.std(resid)#1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
vmax_resid = 3 * rms

for ia in range(len(cubes)):

    hdu = fits.open(cubes[ia]+'.fits')
    h = hdu[0].header
    hdu.close()

    dx = 3600 * h['CDELT1'] * (np.arange(h['NAXIS1']) - (h['CRPIX1'] - 1))
    dy = 3600 * h['CDELT2'] * (np.arange(h['NAXIS2']) - (h['CRPIX2'] - 1))
    ext = (dx.max(), dx.min(), dy.min(), dy.max())
    bmaj, bmin, bpa = h['BMAJ'], h['BMIN'], h['BPA']

    vv = vv_list[ia]
    j0 = np.argmin(np.abs(vv - 9.0e3))

    for i in range(nchans):
        j = i + j0
        Tb = Tb_all[ia][j]

        ax = axs[ia, i]

        if ia == 2:
            # diverging scale + rms contours for residuals
            im_r = ax.imshow(Tb, origin='lower', cmap='RdBu_r', extent=ext,
                              aspect='auto')#, vmin=-vmax_resid, vmax=vmax_resid)
            ax.contour(Tb, levels=[-5*rms, -3*rms, 3*rms, 5*rms],
                       colors='k', linewidths=0.4,
                       linestyles=['dashed', 'dashed', 'solid', 'solid'],
                       extent=ext, origin='lower')
        else:
            im = ax.imshow(Tb, origin='lower', cmap='cmr.cosmic', extent=ext,
                           aspect='auto', vmin=0, vmax=vmax_data)

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                        -xlims[0] - 0.1*np.diff(xlims)),
                       3600 * bmaj, 3600 * bmin, angle=90-bpa)
        beam.set_facecolor('w' if ia != 2 else 'k')
        ax.add_artist(beam)

        if i == 0:
            ax.text(0.02, 0.90, lbls[ia], transform=ax.transAxes, ha='left',
                    va='center', style='italic', fontsize=8,
                    color='w' if ia != 2 else 'k')
        if ia == 2:
            if np.abs(vv[j]) < 0.001: vv[j] = 0.0
            pref = '+' if np.sign(vv[j]) >= 0 else ''
            vstr = pref + '%.2f' % (1e-3 * vv[j])
            ax.text(0.97, 0.08, vstr, transform=ax.transAxes, ha='right',
                    va='center', fontsize=7, color='k')
        if np.logical_and(ia == 2, i == 0):
            ax.set_xlabel('RA offset  ($^{\prime\prime}$)')
            ax.set_ylabel('DEC offset  ($^{\prime\prime}$)')
            ax.spines[['right', 'top']].set_visible(False)
        else:
            ax.axis('off')
            ax.set_xticks([])
            ax.set_yticks([])

# two colorbars: one for data/model, one for residual
cbax1 = fig.add_axes([fr+0.01, fb+0.35*(ft-fb)+0.02, 0.02, 0.65*(ft-fb)-0.02])
cb1 = Colorbar(ax=cbax1, mappable=im, orientation='vertical', ticklocation='right')
cb1.set_label('$T_{\\rm b}$  (K)', rotation=270, labelpad=13)

cbax2 = fig.add_axes([fr+0.01, fb, 0.02, 0.35*(ft-fb)-0.01])
cb2 = Colorbar(ax=cbax2, mappable=im_r, orientation='vertical', ticklocation='right')
cb2.set_label('resid.  (K)', rotation=270, labelpad=13)

fig.subplots_adjust(left=fl, right=fr, bottom=fb, top=ft, hspace=hs, wspace=ws)
fig.savefig('testdata/demo_compare.png')
fig.clf()