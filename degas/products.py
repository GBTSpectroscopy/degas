import glob
import os
import numpy as np

from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import units as u


def makeMap(cubeFile, outDir, maskFile='',maptype='peakIntensity',order=0):
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
        type of map to make. Options are peakIntensity, linewidth, and moment.

    order: int
        order for moment map

    """

    cube = SpectralCube.read(cubeFile)
    
    if maskFile:
        maskCube = SpectralCube.read(maskFile)
        mask = (maskCube.unitless_filled_data[:,:,:] > 0.5)
        maskedCube = cube.with_mask(mask)
    else:
        maskedCube = cube

    if maptype is 'peakIntensity':
        mymap = maskedCube.max(axis=0)
        outname = os.path.basename(cubeFile).replace('.fits','_peakInt.fits')

    elif maptype is 'linewidth':
        mymap = maskedCube.linewidth_fwhm()
        outname = os.path.basename(cubeFile).replace('.fits','_linewidth.fits')
    elif maptype is 'moment':
        mymap = maskedCube.moment(order=order)
        outname = os.path.basename(cubeFile).replace('.fits','_mom'+str(order)+'.fits')
    else:
        print('Map type unknown. Options are peakIntensity, linewidth, moment')
        return

    mymap.meta['MASK'] = maskFile
    mymap.write(os.pathjoin(outDir,outname),overwrite=True)

