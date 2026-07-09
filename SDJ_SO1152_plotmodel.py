
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


cubes = ['SigOriData/SO1152_template']
lbls = ['model']

nchans = 20

fig, axs_2 = plt.subplots(nrows=2, ncols=int(nchans/2), figsize=(7.5, 2.1), dpi=300)
axs = axs_2.flatten()
fl, fr, fb, ft, hs, ws = 0.06, 0.88, 0.18, 0.99, 0.12, 0.03
xlims = [1.1, -1.1]
ylims = [-1.1, 1.1]

# Load the cube and header data
hdu = fits.open(cubes[0]+'.fits')
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

for i in range(nchans):

    # in-cube index
    j = i + 0

    # convert intensities to brightness temperatures
    Tb = (1e-26 * np.squeeze(Ico[j,:,:]) / bm) * sc.c**2 / \
         (2 * sc.k * ff[j]**2)

    # allocate the panel
    ax = axs[i]
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
    ax.text(0.02, 0.90, lbls[0], transform=ax.transAxes, ha='left',
            va='center', style='italic', fontsize=8, color='w')

    if np.abs(vv[j]) < 0.001: vv[j] = 0.0
    if np.logical_or(np.sign(vv[j]) == 1, np.sign(vv[j]) == 0):
        pref = '+'
    else:
        pref = ''
    vstr = pref+'%.2f' % (1e-3 * vv[j])
    ax.text(0.97, 0.08, vstr, transform=ax.transAxes, ha='right',
            va='center', fontsize=7, color='w')

    if i == 0:
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
fig.savefig('SigOriData/demo_res.png')
fig.savefig('SigOriData/demo_res.pdf')
fig.clf()




