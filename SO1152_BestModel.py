import os
import sys
### IMPORT CASATOOLS
import casatools
import numpy as np
from csalt.model import *
from csalt.helpers import *
import matplotlib as mpl
import json


datapath          = 'SigOriData/'
source_name      = "SO1152"
#priors = importlib.import_module('CSALT0_SDJ_SO1153')
cm = model('CSALT0_SDJ_SO1152')
#-----------------------------------

obs_measurement_set = f'{datapath}{source_name}_12CO.ms.contsub_lsrk_avg'

with open("SigOriData.json", "r") as file:
    loaded_data = json.load(file)
source_data = loaded_data[source_name]
obs_measurement_set = f'{datapath}{source_name}_12CO.ms.contsub_lsrk_avg'
mdict_name_prefix = f"{datapath}{source_name}_12CO"
outname           = f'{datapath}{source_name}_init_walker'
rest_freq_Hz      = 345.7959899e9
rest_freq_GH_str  = '345.7959899GHz'
dist              = loaded_data[source_name]["dist"]       # pc
vra               = [7e3, 1.8e4] # velocity range for fitdata [m/s]
Nwalk             = 1 # number of walkers for MCMC (1 for test)


#########-----CHANGE MODEL PARAMETERS BELOW-----#########
# Disk geometry      # stellar mass (solar masses)

#inc     = loaded_data[source_name]["incl"]   # inclination (degrees)
#PA      = loaded_data[source_name]["PA"]+180       # position angle (degrees)
#mstar   = loaded_data[source_name]["mass"]  
#vlsr     = loaded_data[source_name]["vlsr"] 


inc = 66.406
PA = 324.438
mstar = 0.534
z_1 = 0.157
z_q = 0.691
Tb_10 = 56.880
Tb_q = -0.188
z_tap = 111.756
z_exp = 4.945
logtau_10 = 4.869
tau_q = 0.938
tau_tap = 148.631
tau_exp = 9.535
vlsr = 12162.790
dx = -0.017
dy = 0.007


pars = np.array([inc, PA, mstar, z_1, z_q, Tb_10, Tb_q, z_tap, z_exp, logtau_10, tau_q,tau_tap, tau_exp, vlsr, dx, dy])


# Tclean parameters
imsize     = 500      # image size (pixels)
start_vel  = '7.0km/s'       # starting velocity for imaging
width      = '0.5km/s'    # velocity width for imaging
nchan      = 30      # number of channels to image
rest_freq  = '345.7959899GHz'  # rest frequency of the line (12CO 3-2)
rest_freq_Hz  = 345.7959899e9  # rest frequency of the line (12CO 3-2)
cell_size  = '0.02arcsec'  # cell size for imaging
scales     =  [0, 10, 20, 50]
niter      = 50000
robust_val =  0.5
thres      = '2mJy'  

#  Keplerian mask parameters
r_max  = 1.0       # maximum radius for the mask (arcsec)
nbeam  = 1.0       # number of beams for the mask
zr    = 0.0      # height of the mask (arcsec)


tclean_kw = {'imsize': imsize, 'start': start_vel, 'width': width,
            'nchan': nchan, 'restfreq': rest_freq, 'cell': cell_size,
            'scales': scales, 'niter': niter,
            'robust': robust_val, 'threshold': thres}#

# Define some Keplerian mask keywords
kepmask_kw = {'inc': inc, 'PA': PA,
                'mstar': mstar, 'dist': dist, 'vlsr': vlsr,
            'r_max': r_max, 'nbeams': nbeam, 'zr': zr}

# Define some fixed attributes for the model
fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': dist,'restfreq': rest_freq_Hz,
            'Nup': 5, 'doppcorr': 'approx', 'online_avg': 4}

### clean the observations 
imagecube(f'{obs_measurement_set}', f'{mdict_name_prefix}_Observation',
        kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

cubes = [f'{mdict_name_prefix}_Observation']
lbls = ['Observation']

# Export the cubes to FITS format
from casatasks import exportfits
for i in range(len(cubes)):
    exportfits(cubes[i]+'.image', fitsimage=cubes[i]+'.fits',
               velocity=True, overwrite=True)
    

      
# Define velocity axis — same way modelset does it internally
# (you need to pick a velocity range covering your emission)
#vlsr   = vlsr                        # systemic velocity [m/s]
dv     = 500.                              # channel width [m/s]
nchan  = 23                               # number of channels
velax  = np.arange(7000, 21500, 500) # [m/s]

# Call cube() directly — returns SkyImage before vis_sample touches it
icube = cm.cube(velax, pars,
                restfreq=rest_freq_Hz,
                FOV=5.11,
                Npix=512,
                dist=dist)

# icube is a vis_sample SkyImage with these attributes:
print(type(icube))
print('data shape:    ', icube.data.shape)   # (npix, npix, nchan)
print('RA axis:       ', icube.ra.shape)     # arcsec offsets
print('Dec axis:      ', icube.dec.shape)
print('freq axis:     ', icube.freqs.shape)  # Hz


# Write it to FITS yourself
cube_to_fits(icube, 'SO1152_InitModel_skyImage.fits',
            restfreq=rest_freq_Hz)


# Get dataset measurement set
ddict = read_MS(obs_measurement_set) 
pure_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(pure_mdict, outfile=f'{mdict_name_prefix}_InitMODEL.ms')

imagecube(f'{mdict_name_prefix}_InitMODEL.ms', f'{mdict_name_prefix}_InitMODEL',
        kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)



write_MS(pure_mdict, outfile=f'{mdict_name_prefix}_InitMODEL_RESID.ms', resid=True)

imagecube(f'{mdict_name_prefix}_InitMODEL_RESID.ms', 
          outfile=f'{mdict_name_prefix}_InitMODEL_RESID',
          mk_kepmask=False, tclean_kwargs=tclean_kw)

cubes = [f'{mdict_name_prefix}_InitMODEL', f'{mdict_name_prefix}_InitMODEL_RESID']
lbls = ['InitMODEL', 'InitMODEL_RESID']

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
import numpy as np 
#_ = importlib.import_module('plot_setups')
#plt.style.use(['default', '/home/sandrews/mpl_styles/nice_img.mplstyle'])
import cmasher as cmr


cubes = [f'{mdict_name_prefix}_Observation', 
         f'{mdict_name_prefix}_InitMODEL', 
         f'{mdict_name_prefix}_InitMODEL_RESID']

lbls = ["Observations", 'InitMODEL', "Residuals"]

# per-row color limits (None -> set from the data)
vmins = [0, 0, None]
vmaxs = [25, 25, None]

# residual contour settings
nsig = 3                         # contour level in units of sigma (drawn at +/- nsig)
sigma_manual = None              # set a value in K to override the measured sigma
r_noise = 3.0                    # arcsec; pixels beyond this radius are used for the noise

nrows = len(cubes)               # rows per block (Obs / Model / Resid)
nchans = 20                      # total channels to show
ncols = 10                       # columns per block (top block = first ncols channels)
nblocks = int(np.ceil(nchans / ncols))   # 2 blocks: top and bottom

# figure layout (fractions of the figure)
figw, figh = 12, 7.0
fl, fr, fb, ft = 0.06, 0.88, 0.09, 0.99
ws = 0.03                        # horizontal space between channels
gap = 0.05                       # vertical gap between the top and bottom blocks
blk_h = (ft - fb - gap) / nblocks
row_h = blk_h / nrows

xlims = [1.1, -1.1]
ylims = [-1.1, 1.1]

fig = plt.figure(figsize=(figw, figh))

# one GridSpec per block; hspace=0 -> no space between rows within a block
gss = []
for b in range(nblocks):
    top = ft - b * (blk_h + gap)
    bot = top - blk_h
    gss.append(fig.add_gridspec(nrows=nrows, ncols=ncols, left=fl, right=fr,
                                top=top, bottom=bot, hspace=0, wspace=ws))

ims = [None] * nrows             # one mappable per row for the colorbars
sigma = None                     # residual noise (K), set in the last row

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

    # convert all displayed channels to brightness temperatures
    Tb_all = (1e-26 * Ico[:nchans] / bm) * sc.c**2 / \
             (2 * sc.k * ff[:nchans, None, None]**2)

    # row color limits: fixed if given, otherwise from the data in this row
    vmin = vmins[ia] if vmins[ia] is not None else np.nanmin(Tb_all)
    vmax = vmaxs[ia] if vmaxs[ia] is not None else np.nanmax(Tb_all)

    # noise estimate for the residual row from the emission-free outskirts
    if ia == nrows - 1:
        if sigma_manual is not None:
            sigma = sigma_manual
        else:
            rr = np.hypot(dx[None, :], dy[:, None])
            noise_mask = rr > r_noise
            if noise_mask.sum() == 0:    # image too small: fall back to a robust global estimate
                med = np.nanmedian(Tb_all)
                sigma = 1.4826 * np.nanmedian(np.abs(Tb_all - med))
            else:
                sigma = np.nanstd(Tb_all[:, noise_mask])
        levels = np.array([-nsig, nsig]) * sigma
        print(f'Residual noise: {sigma:.3f} K  ->  contours at +/-{nsig} sigma = '
              f'+/-{nsig*sigma:.3f} K')

    for b in range(nblocks):
        for i in range(ncols):

            # in-cube index
            j = b * ncols + i

            # allocate the panel
            ax = fig.add_subplot(gss[b][ia, i])

            # blank out unused panels in the last block
            if j >= nchans:
                ax.axis('off')
                continue

            Tb = Tb_all[j]

            # plot the channel map
            im = ax.imshow(Tb, origin='lower', cmap='cmr.cosmic', extent=ext,
                           aspect='auto', vmin=vmin, vmax=vmax)
            ims[ia] = im

            # +/- nsig*sigma contours on the residual row (solid = +, dashed = -)
            if ia == nrows - 1:
                lev = levels[(levels > np.nanmin(Tb)) & (levels < np.nanmax(Tb))]
                if lev.size > 0:
                    ax.contour(dx, dy, Tb, levels=lev, colors='w',
                               linewidths=0.5, negative_linestyles='dashed')

            # set map boundaries
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

            # overlay beam dimensions
            beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                            -xlims[0] - 0.1*np.diff(xlims)),
                           3600 * bmaj, 3600 * bmin, angle=90-bpa)
            beam.set_facecolor('w')
            ax.add_artist(beam)

            # labeling: row labels at the start of each block
            if i == 0:
                ax.text(0.02, 0.90, lbls[ia], transform=ax.transAxes,
                        ha='left', va='center', style='italic', fontsize=8,
                        color='w')

            # velocity labels on the bottom row of each block
            if ia == nrows - 1:
                if np.abs(vv[j]) < 0.001: vv[j] = 0.0
                if np.logical_or(np.sign(vv[j]) == 1, np.sign(vv[j]) == 0):
                    pref = '+'
                else:
                    pref = ''
                vstr = pref+'%.2f' % (1e-3 * vv[j])
                ax.text(0.97, 0.08, vstr, transform=ax.transAxes, ha='right',
                        va='center', fontsize=7, color='w')

            # axes only on the bottom-left panel (last row, last block)
            if (ia == nrows - 1) and (b == nblocks - 1) and (i == 0):
                ax.set_xlabel(r'RA offset  ($^{\prime\prime}$)')
                ax.set_ylabel(r'DEC offset  ($^{\prime\prime}$)')
                ax.spines[['right', 'top']].set_visible(False)
            else:
                ax.axis('off')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])

# colorbars: for each block, one for Obs+Model rows and one for the Residuals row
cbw = 0.015
for b in range(nblocks):
    btop = ft - b * (blk_h + gap)

    # Tb colorbar spanning all rows except the last
    h_data = (nrows - 1) * row_h
    cbax1 = fig.add_axes([fr+0.01, btop - h_data + 0.005, cbw, h_data - 0.01])
    cb1 = Colorbar(ax=cbax1, mappable=ims[0], orientation='vertical',
                   ticklocation='right')
    cb1.set_label(r'$T_{\rm b}$  (K)', rotation=270, labelpad=13)

    # residual colorbar aligned with the last row
    cbax2 = fig.add_axes([fr+0.01, btop - blk_h + 0.005, cbw, row_h - 0.01])
    cb2 = Colorbar(ax=cbax2, mappable=ims[-1], orientation='vertical',
                   ticklocation='right')
    cb2.set_label(r'$T_{\rm b}$  (K)', rotation=270, labelpad=13)

# contour note
fig.text(fr+0.01, 0.01, rf'Residual contours: $\pm{nsig}\sigma$, '
         rf'$\sigma = {sigma:.2f}$ K', ha='left', va='bottom', fontsize=7)

fig.savefig(f'SigOriData/{source_name}_compare.png')
fig.savefig(f'SigOriData/{source_name}_compare.pdf')
fig.clf()