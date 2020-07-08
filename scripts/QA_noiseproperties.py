from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import glob
import sys
sys.path.append('/mnt/space/erosolow/phangs_imaging_scripts')
from phangsPipeline import scMaskingRoutines as scm
import numpy as np
import scipy.stats as ss
import os
import astropy.units as u
import astropy.wcs as wcs
import copy
degas_data = '/mnt/bigdata/erosolow/surveys/DEGAS/'


fl = glob.glob(degas_data + '/IR5/*fits')

thco = [f for f in fl if '13CO' in f]
c18o = [f for f in fl if 'C18O' in f]
hcn = [f for f in fl if 'HCN' in f]
hcop = [f for f in fl if 'HCOp' in f]

for thisfile in fl:
    fig, [[ax1, ax2],[ax3, ax4]] = plt.subplots(2, 2, constrained_layout=True)
    fig.set_size_inches(6.5, 6.5)
    cube = SpectralCube.read(thisfile)
    gal = os.path.split(thisfile)[-1].split('_')[0]
    shortfile = os.path.split(thisfile)[-1]

    try:
        mask = SpectralCube.read(degas_data + '/masks/{0}_mask.fits'.format(gal))
        mask = mask.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        hdr = mask.header
        hdr['CTYPE3'] = 'VRAD'
        mask = SpectralCube(mask.filled_data[:], wcs.WCS(hdr), header=hdr)
        mask = mask.with_mask(mask.filled_data[:] > 0)
        mask = mask.spectral_interpolate(cube.spectral_axis)
        mask = mask.reproject(cube.header)

    except FileNotFoundError:
        print('No mask for '+thisfile)
        mask = SpectralCube(np.ones(cube.shape), wcs=cube.wcs, header=cube.header)
        
    rmsamplitude = scm.noise_cube(cube.filled_data[:].value,
                                  box=5,
                                  spec_box=5, iterations=5,
                                  bandpass_smooth_window=21)
    s2n = cube.filled_data[:].value / rmsamplitude

    s2n_in_mask = s2n[(mask.filled_data[:].value == 1)]
    s2n = s2n.ravel()
    bins = np.linspace(-7,15,221)

    s2nhist, _ = np.histogram(s2n, bins)
    npts = np.nansum(s2nhist)
    s2nhist = s2nhist / npts / (bins[1] - bins[0])

    binmid = 0.5*(bins[0:-1] + bins[1:])
    
    ax1.fill_between(binmid, s2nhist, label='All data')
    ax1.plot(binmid, s2nhist, drawstyle='steps', color='black')
    ax1.set_yscale('log')

    s2ninmask_hist, _  = np.histogram(s2n_in_mask, bins) 
    s2ninmask_hist  = s2ninmask_hist / npts / (bins[1] - bins[0])

    ax1.fill_between(binmid, s2ninmask_hist, label='In Mask')
    ax1.plot(binmid, s2ninmask_hist, drawstyle='steps', color='black')

    ax1.plot(binmid, ss.norm.pdf(binmid), color='red', linestyle='--', label='Normal Dist.')
    ax1.set_ylim(1e-5,1)
    ax1.set_xlim(-7, 15)
    ax1.set_xlabel(r'$T_\mathrm{A}^*/ \sigma_T$')
    ax1.set_ylabel('PDF')
    ax1.legend()
    parts = shortfile.split('_')
    ax1.set_title(parts[0] + ' ' + parts[1])
    img = ax2.imshow(rmsamplitude[int(cube.shape[0]/2), :, :]*1e3)
    cb = plt.colorbar(img,ax=ax2)
    cb.set_label(r'$\sigma(x,y)\ (\mathrm{mK})$')
    ax2.set_xlabel(r'$X\ (\mathrm{pix})$')
    ax2.set_ylabel(r'$Y\ (\mathrm{pix})$')

    ax3.plot(cube.spectral_axis.to(u.km/u.s),
             rmsamplitude[:, int(cube.shape[1]/2), int(cube.shape[2]/2)] * 1e3,
             drawstyle='steps', label='Noise amplitude')
    ax3.set_xlabel(r'$V_\mathrm{LSR}\ (\mathrm{km/s})$')
    ax3.set_ylabel(r'$\sigma(v)\ (\mathrm{ mK})$')

    avspec = np.nanmean(cube.filled_data[:].value, axis=(1,2))
    ax4.plot(cube.spectral_axis.to(u.km/u.s),
             avspec * 1e3, 
             drawstyle='steps', label='Mean Spectrum')

    blanked = np.zeros(cube.shape)
    idx = (mask.filled_data[:].value == 1)
    blanked[idx] = cube.filled_data[idx].value
    bad = np.isnan(cube.filled_data[:])
    blanked[bad] = np.nan
    avspec_in_mask = np.nanmean(blanked, axis=(1,2))
    ax4.plot(cube.spectral_axis.to(u.km/u.s),
             avspec_in_mask * 1e3, 
             drawstyle='steps', label='Mean Masked Spectrum')
    ax4.set_xlabel(r'$V_\mathrm{LSR}\ (\mathrm{km/s})$')
    ax3.set_ylabel(r'$T_\mathrm{A}^*\ (\mathrm{ mK})$')
    ax4.legend()

    fig.savefig('_'.join((parts[0], parts[1], 'QAnoiseplot')) + '.png', dpi=300)
    plt.close()
    plt.clf()

