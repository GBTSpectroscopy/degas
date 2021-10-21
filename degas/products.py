import glob
import os
import numpy as np

from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy import units as u


def makeMap(cubeFile, outDir, baseName=None,maskFile='',maptype='peakIntensity',order=0 ):
    """ Create peak intensity map from a cube

    Parameters

    cubeFile: str
       FITS filename of the cube

    outDir: str
        Directory in which to place the output FITS file

    maskFile: str 
        FITS filename of the mask. Should have same geometry as
        cube. If not use routines in analysis_setup.py to regrid.

    maptype: str
        type of map to make. Options are peakIntensity, peakVelocity, linewidth, 
        moment, and mask2D.

    order: int
        order for moment map

    """

    if not baseName:
        baseName=os.path.basename(cubeFile).replace('.fits','')

    cube = SpectralCube.read(cubeFile).with_spectral_unit(u.km/u.s, velocity_convention='radio')

    if maskFile:
        maskCube = SpectralCube.read(maskFile).with_spectral_unit(u.km/u.s, velocity_convention='radio')
        mask = (maskCube.unitless_filled_data[:,:,:] > 0.5)
        maskedCube = cube.with_mask(mask)
    else:
        maskedCube = cube

    if maptype is 'peakIntensity':
        mymap = maskedCube.max(axis=0)
        outname = baseName+'_peakInt.fits'

    elif maptype is 'mask2D':
        mymap = maskedCube.max(axis=0)
        outname = baseName+'_mask2D.fits'

    elif maptype is 'linewidth':
        mymap = maskedCube.linewidth_fwhm()
        outname = baseName+'_linewidth.fits'

    elif maptype is 'moment':
        mymap = maskedCube.moment(order=order)
        outname = baseName+'_mom'+str(order)+'.fits'

    elif maptype is 'peakVelocity':
        velocity = maskedCube.spectral_axis
        mask = np.sum(maskedCube.get_mask_array(),axis=0) > 1.0
        newwcs = maskedCube.wcs.dropaxis(2)

        peakind = maskedCube.argmax(axis=0)
        peakvel = velocity[peakind]         
        peakvel[peakvel==velocity[0]] = np.nan
        
        mymap = Projection(peakvel, wcs=newwcs, mask=mask)
        outname = baseName + '_peakVelocity.fits'
    else:
        print('Map type unknown. Options are peakIntensity, linewidth, moment')
        return

    mymap.meta['MASK'] = maskFile
    mymap.write(os.path.join(outDir,outname),overwrite=True)

