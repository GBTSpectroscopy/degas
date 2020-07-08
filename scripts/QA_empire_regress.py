from spectral_cube import SpectralCube
from astropy.table import Table, hstack, vstack
from astropy.io import fits
from astropy.wcs import wcs
import numpy as np
import astropy.units as u
import corner
from astropy.stats import median_absolute_deviation as mad
import matplotlib.pyplot as plt
import os
import glob
from astropy.convolution import Box1DKernel

sf = 1e3
eta_mb_dict = {'HCN':0.3,
               'HCOp':0.3,
               '13CO':0.2,
               'C18O':0.2}

def sanitize(hdu):
    bunit_check = {'K*KM/S':'K km / s',
                   'K.KM/S':'K km / s',
                   'K km / s':'K km / s',
                   'NONE':'NONE',
                   'K':'K',
                   'KM/S':'km / s',
                   'km /s':'km / s',
                   'km / s':'km / s', 'M/S':'m / s'}
    
    kwlist = ['CTYPE3', 'CRVAL3', 
              'CDELT3', 'CRPIX3', 'CUNIT3', 
              'PC3_1', 'PC3_2', 'PC1_3', 'PC2_3', 'PC3_3']
    data = hdu.data
    hdr = hdu.header
    data = data.squeeze()
    if len(data.shape) == 2:
        for kw in kwlist:
            if kw in hdr:
                hdr.remove(kw, ignore_missing=True)
    else:
        hdr['CUNIT3'] = bunit_check[hdr['CUNIT3']]
        hdr['BUNIT'] = bunit_check[hdr['BUNIT']]
        
    output_hdu = fits.PrimaryHDU(data, header=hdr)
    return(output_hdu)


t = Table.read('degas_base.fits')

empire_gal = ['ngc3627','ngc6946','ngc5194','ngc5055','ngc0628',
              'ngc2903','ngc4321','ngc4254','ngc3184']
empire_gal = ['ngc2903','ngc4321','ngc5055','ngc6946'] 
empire_data = '/mnt/bigdata/erosolow/surveys/EMPIRE/EMPIRE_cubes/'
degas_data = '/mnt/bigdata/erosolow/surveys/DEGAS/'

subtable = vstack([row for row in t if row['NAME'].lower() in empire_gal])


molecule_list = ['13CO','C18O','HCN','HCOp']

for molecule in molecule_list:
    fig, axlist = plt.subplots(2, 2, constrained_layout=True)
    fig.set_size_inches(6.5, 6.5)
    for gal, ax in zip(subtable['NAME'], axlist.flatten()):
        eta_mb = eta_mb_dict[molecule]
        degas_cube_name = 'IR5/{0}_{1}_rebase3_smooth1.3_hanning1.fits'.format(gal, molecule)
        if not os.path.exists(degas_data + degas_cube_name):
            ax.set_title(gal)
            ax.text(0.5, 0.5, 'No data',ha='center')
            continue
        cube_degas = SpectralCube.read(degas_data + degas_cube_name)
        cube_degas = cube_degas / eta_mb
        mask = SpectralCube.read(degas_data + '/masks/{0}_mask.fits'.format(gal))
        mask = mask.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        hdr = mask.header
        hdr['CTYPE3'] = 'VRAD'
        mask = SpectralCube(mask.filled_data[:], wcs.WCS(hdr), header=hdr)
        mask = mask.with_mask(mask.filled_data[:] > 0)
        mask = mask.spectral_interpolate(cube_degas.spectral_axis)
        mask = mask.reproject(cube_degas.header)

        fl = glob.glob(empire_data + 'EMPIRE_{0}_{1}_*.fits'.format(gal.lower(),molecule.lower()))
        hdulist = fits.open(fl[0])
        hdu = sanitize(hdulist[0])
        cube_empire = SpectralCube(data=hdu.data, header=hdu.header, wcs=wcs.WCS(hdu.header))
        empire_mask = np.isfinite(hdu.data)
        cube_empire = cube_empire.with_mask(empire_mask)

        channel_ratio = ((cube_degas.spectral_axis[1]
                          - cube_degas.spectral_axis[0]) / 
                         (cube_empire.spectral_axis[1]
                          - cube_empire.spectral_axis[0])).to(u.dimensionless_unscaled).value
        
        kernel = Box1DKernel(channel_ratio)
        cube_empire = cube_empire.spectral_smooth(kernel)
        cube_empire = cube_empire.spectral_interpolate(cube_degas.spectral_axis)
        cube_empire = cube_empire.reproject(cube_degas.header)
        cube_degas = cube_degas.convolve_to(cube_empire.beam)

        noise_empire = mad(cube_empire.filled_data[:].value, ignore_nan=True)
        noise_degas = mad(cube_degas.filled_data[:].value, ignore_nan=True)

        bools = mask.filled_data[:] > 1
        degas_vals = cube_degas.filled_data[bools].value
        empire_vals = cube_empire.filled_data[bools].value
        idx = np.isfinite(degas_vals) * np.isfinite(empire_vals)

        ax.plot(empire_vals * sf, degas_vals * sf,'ro', alpha=0.5, markeredgecolor='k',
                label=molecule)
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()

        ax.plot([np.min([xlims[0],ylims[0]]), np.max([xlims[1], ylims[1]])],
                [np.min([xlims[0],ylims[0]]), np.max([xlims[1], ylims[1]])])
        ax.set_xlabel('EMPIRE (mK)')
        ax.set_ylabel('DEGAS (mK)')
        ax.set_title(gal)
        ax.grid(linestyle='dotted')
        ax.set_aspect('equal')
        ax.legend()

        ax.fill_between(ax.get_xlim(), -noise_degas * sf, noise_degas * sf, alpha=0.4, color='gray')
        ax.fill_betweenx(ax.get_ylim(), -noise_empire * sf, noise_empire * sf, alpha=0.4, color='blue')
    fig.savefig('QA_empirecompre_{0}.pdf'.format(molecule))
    plt.close()
    plt.clf()
    
