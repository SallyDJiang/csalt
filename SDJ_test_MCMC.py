import multiprocessing
multiprocessing.set_start_method('fork')
import numpy as np
import copy
from csalt.model import *
from csalt.helpers import *
import importlib
import os
import casatools
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
from matplotlib.patches import Ellipse
from astropy.io import fits
import scipy.constants as sc
import cmasher as cmr
from casatasks import exportfits

np.random.seed(100)

# SETUP — change these to match your source
datapath          = 'SigOriData/'
obs_measurement_set            = f'{datapath}SO1152_12CO.ms.contsub_lsrk_avg'
mdict_name_prefix = f"{datapath}SO1152_12CO"
outname           = f'{datapath}SO1152_init_walker'
rest_freq_Hz      = 345.7959899e9
rest_freq_GH_str  = '345.7959899GHz'
dist              = 391.3        # pc
vra               = [7e3, 1.8e4] # velocity range for fitdata [m/s]
Nwalk             = 1 # number of walkers for MCMC (1 for test)

fixed_kw = {'FOV': 5.11, 'Npix': 512, 'dist': dist,
            'restfreq': rest_freq_Hz, 'Nup': 5, 'doppcorr': 'approx'}

# tclean parameters
tclean_kw = {'imsize': 500, 'start': '7.0km/s', 'width': '1.0km/s',
             'nchan': 10, 'restfreq': rest_freq_GH_str, 'cell': '0.02arcsec',
             'scales': [0, 10, 20, 50], 'niter': 50000,
             'robust': 0.5, 'threshold': '3mJy'}

# set up prior 



priors = importlib.import_module('priors_CSALT0')
Ndim = len(priors.pri_pars)
p0 = np.empty((Nwalk, Ndim))

name_to_idx = {name: i for i, name in enumerate(priors.priors_dict)} 
        
for ix in range(Ndim):
    if priors.pri_types[ix] == "linewidth":
                temp_idx = name_to_idx["Tb_10"] # what id is Tb_10, in case not ix=9

                Tmin, Tmax = priors.pri_pars[temp_idx]
            
                thermal_width_min = np.sqrt(
                    2 * sc.k * Tmin /
                    (28 * (sc.m_p + sc.m_e))
                    )
                thermal_width_max = np.sqrt(
                    2 * sc.k * Tmax /
                    (28 * (sc.m_p + sc.m_e))
                    )

                #sample randomly from deltav_min to deltav_max 
                p0[:, ix] = np.random.uniform(
                    low=thermal_width_min,
                    high=thermal_width_max,
                    size=Nwalk)
                
#            if ix == 9:
#                p0[:,ix] = np.sqrt(2 * sc.k * p0[:,6] / (28 * (sc.m_p+sc.m_e)))
    else:
                _ = [str(priors.pri_pars[ix][ip])+', '
                     for ip in range(len(priors.pri_pars[ix]))]
                cmd = 'np.random.'+priors.pri_types[ix]+ \
                      '('+"".join(_)+str(Nwalk)+')'
                try: 
                    p0[:,ix] = eval(cmd)
                except AttributeError: 
                    pass 

theta = p0[0]


# Print the init parameters clearly
param_names = list(priors.priors_dict.keys())
print("--- First walker initial parameters ---")
for name, val in zip(param_names, theta):
    print(f"  {name:>12s} = {val:.4f}")


# Prepare inference data so log_likelihood works
cm = model('CSALT0')
infdata = cm.fitdata(obs_measurement_set, vra=[7e3, 1.8e4], vcensor=None, 
                     restfreq=rest_freq_Hz, online_avg=1,
                     chbin=2, well_cond=300)


# Set globals that log_likelihood needs
import csalt.model as cm_module
cm_module.fdata = copy.deepcopy(infdata)
cm_module.kw    = copy.deepcopy(fixed_kw)
cm_module.kw.setdefault('chpad', 2)
cm_module.kw.setdefault('noise_inject', None)
cm_module.kw.setdefault('SRF', 'ALMA')
cm_module.kw.setdefault('online_avg', 1)

lnT = np.sum(priors.logprior(theta)) * infdata['Nobs']
lnL = cm.log_likelihood(theta, fdata=infdata, kwargs=cm_module.kw)
print(f"ln prior:      {lnT:.4f}")
print(f"ln likelihood: {lnL:.4f}")
print(f"ln posterior:  {lnT + lnL:.4f}\n")

#if not np.isfinite(lnT):
#    print("WARNING: ln prior is -inf — this walker starts outside the prior.")
#    print("         Check your priors_CSALT0.py bounds.\n")
#if not np.isfinite(lnL):
#    print("WARNING: ln likelihood is -inf — model failed at these parameters.\n")

# Generate the forward model at theta
ddict = read_MS(obs_measurement_set)
mdict = cm.modeldict(ddict, theta, kwargs=fixed_kw)
write_MS(mdict, outfile=outname + '_MODEL.ms')
print(f"Model MS written to {outname}_MODEL.ms")

# Keplerian mask keywords from theta0
kepmask_kw = {'inc':   theta[0],
              'PA':    theta[1],
              'mstar': theta[2],
              'dist':  dist,
              'vlsr':  theta[12],
              'r_max': 1.0, 'nbeams': 1.0, 'zr': 0.2}

# Image the model
imagecube(outname + '_MODEL.ms', outname + '_MODEL',
          kepmask_kwargs=kepmask_kw, tclean_kwargs=tclean_kw)

# Export to FITS
exportfits(outname + '_MODEL.image', fitsimage=outname + '_MODEL.fits',
           velocity=True, overwrite=True)
print(f"FITS written to {outname}_MODEL.fits")


# Plot the channel maps
nchans = tclean_kw['nchan']
fig, axs = plt.subplots(nrows=1, ncols=nchans, figsize=(1.1 * nchans, 1.5))
fl, fr, fb, ft = 0.06, 0.88, 0.18, 0.99
xlims = [1.5, -1.5]
ylims = [-1.5, 1.5]

hdu = fits.open(outname + '_MODEL.fits')
Ico, h = np.squeeze(hdu[0].data), hdu[0].header
hdu.close()

dx_ax = 3600 * h['CDELT1'] * (np.arange(h['NAXIS1']) - (h['CRPIX1'] - 1))
dy_ax = 3600 * h['CDELT2'] * (np.arange(h['NAXIS2']) - (h['CRPIX2'] - 1))
ext   = (dx_ax.max(), dx_ax.min(), dy_ax.min(), dy_ax.max())
bmaj, bmin, bpa = h['BMAJ'], h['BMIN'], h['BPA']
bm = (np.pi * bmaj * bmin / (4 * np.log(2))) * (np.pi / 180)**2

vv = h['CRVAL3'] + h['CDELT3'] * (np.arange(h['NAXIS3']) - (h['CRPIX3'] - 1))
ff = rest_freq_Hz * (1 - vv / sc.c)

for i in range(nchans):
    Tb = (1e-26 * np.squeeze(Ico[i, :, :]) / bm) * sc.c**2 / \
         (2 * sc.k * ff[i]**2)

    ax = axs[i]
    im = ax.imshow(Tb, origin='lower', cmap='cmr.cosmic', extent=ext,
                   aspect='auto', vmin=0, vmax=25)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    beam = Ellipse((xlims[0] + 0.1 * np.diff(xlims),
                    -xlims[0] - 0.1 * np.diff(xlims)),
                   3600 * bmaj, 3600 * bmin, angle=90 - bpa)
    beam.set_facecolor('w')
    ax.add_artist(beam)

    vstr = ('+' if vv[i] >= 0 else '') + '%.2f' % (1e-3 * vv[i])
    ax.text(0.97, 0.08, vstr, transform=ax.transAxes,
            ha='right', va='center', fontsize=7, color='w')

    if i == 0:
        ax.text(0.02, 0.90, 'init walker', transform=ax.transAxes,
                ha='left', va='center', style='italic', fontsize=7, color='w')
        ax.set_xlabel('RA offset ($^{\\prime\\prime}$)')
        ax.set_ylabel('DEC offset ($^{\\prime\\prime}$)')
        ax.spines[['right', 'top']].set_visible(False)
    else:
        ax.axis('off')

cbax = fig.add_axes([fr + 0.01, fb + 0.01, 0.02, ft - fb - 0.02])
cb = Colorbar(ax=cbax, mappable=im, orientation='vertical', ticklocation='right')
cb.set_label('$T_{\\rm b}$ (K)', rotation=270, labelpad=13)

fig.subplots_adjust(left=fl, right=fr, bottom=fb, top=ft, hspace=0.12, wspace=0.03)
fig.savefig(outname + '_MODEL.pdf', bbox_inches='tight')
fig.savefig(outname + '_MODEL.png', dpi=150, bbox_inches='tight')
print(f"Plot saved to {outname}_MODEL.png")




#### TEST 2: PRIORS ARE WELL_DSTRIBUTED, AND THE LIKELIHOOD IS FINITE
datapath          = 'SigOriData/'
obs_measurement_set            = f'{datapath}SO1152_12CO.ms.contsub_lsrk_avg'
mdict_name_prefix = f"{datapath}SO1152_12CO"
outname           = f'{datapath}SO1152_init_walker'
rest_freq_Hz      = 345.7959899e9
rest_freq_GH_str  = '345.7959899GHz'
dist              = 391.3        # pc
vra               = [7e3, 1.8e4] # velocity range for fitdata [m/s]
Nwalk             = 10000 # number of walkers for MCMC (1 for test)

# set up prior 
def sample_prior(priors, Nwalk):
    Ndim = len(priors.pri_pars) 
    p0 = np.empty((Nwalk, Ndim))
    name_to_idx = {name: i for i, name in enumerate(priors.priors_dict)} 

    for ix in range(Ndim):
        if priors.pri_types[ix] == "linewidth":
            temp_idx = name_to_idx["Tb_10"] # what id is Tb_10, in case not ix=9

            Tmin, Tmax = priors.pri_pars[temp_idx]
        
            thermal_width_min = np.sqrt(
                2 * sc.k * Tmin /
                (28 * (sc.m_p + sc.m_e))
                )
            thermal_width_max = np.sqrt(
                2 * sc.k * Tmax /
                (28 * (sc.m_p + sc.m_e))
                )

            #sample randomly from deltav_min to deltav_max 
            p0[:, ix] = np.random.uniform(
                low=thermal_width_min,
                high=thermal_width_max,
                size=Nwalk)
                
#            if ix == 9:
#                p0[:,ix] = np.sqrt(2 * sc.k * p0[:,6] / (28 * (sc.m_p+sc.m_e)))
    else:
        _ = [str(priors.pri_pars[ix][ip])+', '
                for ip in range(len(priors.pri_pars[ix]))]
        cmd = 'np.random.'+priors.pri_types[ix]+ \
                '('+"".join(_)+str(Nwalk)+')'
        try: 
            p0[:,ix] = eval(cmd)
        except AttributeError: 
            pass 
    return p0


priors = importlib.import_module('priors_CSALT0')
Nwalks=10000
sample_prior(priors, Nwalk)


theta = p0[0]

import corner
param_names = list(priors.priors_dict.keys())

fig = corner.corner(
    p0,
    labels=param_names,   # optional
    show_titles=True
)

fig.savefig("SigOriData/initial_walker_distribution_test.png", dpi=300, bbox_inches="tight")
