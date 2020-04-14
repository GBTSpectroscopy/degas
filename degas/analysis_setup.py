import glob
import os
import numpy as np

from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import units as u

from reproject import reproject_interp
from astropy import wcs

def fixBIMA(bimafile):
    '''

    Purpose: fix up the bima cubes

    The smoothed BIMA data doesn't have the outer regions
    appropriately masked. This routine copies the mask from the
    original data to the smoothed data. 

    They also don't have all the header information correct.

    '''

    #get smoothed bima data, bima cube has 4 axes
    cube = fits.getdata(bimafile)[0,:,:,:]
    cubehead = fits.getheader(bimafile)

    #get the non-smoothed data, bima cubes have 4 axes
    oldfile = bimafile.replace('_gauss15.fits','_k.fits')
    oldcube = fits.getdata(oldfile)[0,:,:,:]
    oldmask = oldcube==0.0 #extract the original '0'valued pb mask
    cube[oldmask] = np.nan #where the '0' mask was, mask as nan

    #now data is 3d delete info on 4th axies ('stokes')
    cubehead.remove('CDELT4')
    cubehead.remove('CRPIX4')
    cubehead.remove('CRVAL4')
    cubehead.remove('CTYPE4')

     # convert BIMA HEADERS from m/s to km/s
    if 'CUNIT3' not in cubehead:
        cubehead['CDELT3'] = (cubehead['CDELT3']/1000.0)
        cubehead['CRVAL3'] = (cubehead['CRVAL3']/1000.0)
        cubehead['CUNIT3'] = 'km/s'

    # add CUNIT to RA and DEC axis for BIMA headers.
    if 'CUNIT2' not in cubehead:
        cubehead['CUNIT2'] = 'deg'

    if 'CUNIT1' not in cubehead:
        cubehead['CUNIT1'] = 'deg'

    # Set BUNIT to appropriate value. Original in Jy/beam, but I think
    # it's mis-labeled.
    cubehead['BUNIT'] = 'K'

    fixedcube = fits.PrimaryHDU(cube, cubehead)#fixed cube will have 3axes only
    fixedcube.writeto(bimafile.replace('.fits','_fixed.fits'),overwrite=True)

def fixHeracles(fitsimage):
    '''
 
    [XXX IS THIS ROUTINE STILL USEFUL?]

    fixes up image headers for missing values etc so that we can
    process things appropriately en mass.

 
    '''

    # Using the astropy.io.fits to modify headers 
    # since SpectralCube disapproves of modifying headers as far 
    # as I can tell.
    f = fits.open(fitsimage)

    # Convert VELOCITY to VRAD for Heracles headers.
    if f[0].header['CTYPE3'] == 'VELOCITY':
        f[0].header['CTYPE3'] = 'VRAD' 

    f.writeto(fitsimage.replace('.fits','_fixed.fits'),overwrite=True)

def smoothCube(fitsimage,beam=15.0):
    '''
    Smoothes input cube to given resolution
    
    beam is the beam to smooth to in arcsec.

    '''

    from radio_beam import Beam

    newFits = fitsimage.replace(".fits","_smooth.fits")

    newBeam = Beam(beam*u.arcsec)
    
    cube = SpectralCube.read(fitsimage)
    smoothCube = cube.convolve_to(newBeam)
    
    smoothCube.write(newFits,overwrite=True)

def regridData(baseCubeFits, otherDataFits):
    '''
    regrids one data set to match the wcs of the base data set, which
    is assumed to be a cube. The regridded data set can be either 2d
    or 3d.
    '''

    # open the base cube
    try:
        baseCube = SpectralCube.read(baseCubeFits)
    except:
        print("Can't read in " + baseCubeFits + ".")

    # determine how many dimensions the other data sets have
    f = fits.open(otherDataFits)
    ndim = f[0].header['NAXIS']
    f.close()

    # output image name
    newFits = otherDataFits.replace('.fits','_regrid.fits')

    # now regrid images appropriately
    if ndim == 3:

        # read in cube
        otherCube = SpectralCube.read(otherDataFits)

        # interpolate velocity axis. This needs to be done first.
        regridCube = otherCube.spectral_interpolate(baseCube.spectral_axis)
        
        # interpolate the spatial axes
        newCube = regridCube.reproject(baseCube.header)
        newCube.write(newFits,overwrite=True)

    elif ndim == 2:

        # regrid image
        newcube = reproject_interp(otherDataFits,
                                   output_projection=baseCube.wcs.dropaxis(2),
                                   shape_out=baseCube.wcs.dropaxis(2).pixel_shape, ## Is this right might need to flip 
                                   order='nearest-neighbor',
                     
                                   return_footprint=False)

        # write out regridded data
        fits.writeto(newFits,newcube,baseCube.wcs.dropaxis(2).to_header(),overwrite=True)
        

    else:
        print("Number of dimensions of other data set is not 2 or 3.")
