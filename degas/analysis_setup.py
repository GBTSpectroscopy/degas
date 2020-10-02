import glob
import os
import numpy as np

from astropy.io import fits
from astropy import units as u

from spectral_cube import SpectralCube
from radio_beam import Beam
from astropy.convolution import Box1DKernel

from reproject import reproject_interp
from astropy import wcs

import re

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

    #now data is 3d delete info on 4th axis ('stokes')
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
    fixedcube.writeto(bimafile.replace('.fits','_fixhdr.fits'),overwrite=True)

    # Now trying to change VOPT to VRADIO via SpectralCube
    cube2 = SpectralCube.read(bimafile.replace('.fits','_fixhdr.fits'))
   
    ## change to frequency and interpolate to equal frequency intervals
    cube2freq = cube2.with_spectral_unit(u.GHz,rest_value=115.27120180*u.GHz)
    freq_axis = cube2freq.spectral_axis
    chan_width = np.mean(freq_axis[0:-2] - freq_axis[1:-1])
    new_freq_axis = freq_axis[0] - chan_width * np.arange(0,len(freq_axis))
    cubenewfreq = cube2freq.spectral_interpolate(new_freq_axis)
    
    # convert to radio
    cube2radio = cubenewfreq.with_spectral_unit(u.km/u.s,rest_value=115.27120180*u.GHz,velocity_convention='radio')
    cube2radio.write(bimafile.replace('.fits','_fixed.fits'),overwrite=True)

def fixHeracles(fitsimage):
    '''

    fixes up image headers for missing values etc so that we can
    regrid appropriately.
 
    '''

    # Using the astropy.io.fits to modify headers 
    # since SpectralCube disapproves of modifying headers as far 
    # as I can tell.
    f = fits.open(fitsimage)

    # Convert VELOCITY to VRAD for Heracles headers.
    if f[0].header['CTYPE3'] == 'VELOCITY':
        f[0].header['CTYPE3'] = 'VRAD' 

    f.writeto(fitsimage.replace('.fits','_fixed.fits'),overwrite=True)

def fixOVRO(fitsimage,beam=15.0):
    '''
    fix up OVRO data.
    '''
    
    # open image
    cube = SpectralCube.read(fitsimage)    
    
    # convert to K
    kcube = cube.to(u.K)
    
    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = kcube.convolve_to(newBeam)
    
    # write out
    smoothCube.write(fitsimage.replace('.fits','_gauss15_fixed.fits'),
                     overwrite=True)
    

def fixNRO(fitsimage):
    '''
    fix up NRO data
    '''
    
    cube = SpectralCube.read(fitsimage)
    # add beam
    beam_cube = cube.with_beam(Beam(15.0*u.arcsec))
    beam_cube.write(fitsimage.replace('.FITS','_fixed.fits'),overwrite=True)
    

def fixExtraHERA(fitsimage,beam=15.0):
    '''
    Fix the extra heracles data.

    TBD: need to add an efficiency correction to these data.
    '''
    
    f = fits.open(fitsimage)
    f[0].header['CUNIT3'] = 'm/s'
    newimage = fitsimage.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)
    f.close()
    
    # open image
    cube = SpectralCube.read(newimage)    
    
    if re.search('ngc3631',newimage):
        # extra channels with data
        subCube = cube[72:379,:,:]
    elif re.search('ngc4030',newimage):
        subCube = cube[78:384,:,:]
    else:
        subCube = cube[:,:,:]

    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = subCube.convolve_to(newBeam)
    
    # smooth velocity
    smoothFactor = 4.0
    spSmoothCube = smoothCube.spectral_smooth(Box1DKernel(smoothFactor))

    # Interpolate onto a new axis
    spec_axis = spSmoothCube.spectral_axis
    chan_width = spec_axis[1]-spec_axis[0] # channels are equally spaced in velocity
    new_axis = np.arange(spec_axis[0].value,spec_axis[-1].value,smoothFactor*chan_width.value) * u.m/u.s

    interpCube = spSmoothCube.spectral_interpolate(new_axis,
                                                   suppress_smooth_warning=True)

    # write out
    interpCube.write(newimage.replace('.fits','_10kms_gauss15.fits'),
                     overwrite=True)

def smoothCube(fitsimage,outDir, beam=15.0):
    '''
    Smoothes input cube to given resolution
    
    beam is the beam to smooth to in arcsec.

    '''

    from radio_beam import Beam

    newFits = os.path.join(outDir,os.path.basename(fitsimage).replace(".fits","_smooth.fits"))

    newBeam = Beam(beam*u.arcsec)
    
    cube = SpectralCube.read(fitsimage)
    smoothCube = cube.convolve_to(newBeam)
    
    smoothCube.write(newFits,overwrite=True)
    return newFits

def regridData(baseCubeFits, otherDataFits, outDir, mask=False):
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
    newFits = os.path.join(outDir,os.path.basename(otherDataFits).replace('.fits','_regrid.fits'))

    # now regrid images appropriately
    if ndim == 3:

        # read in cube
        otherCube = SpectralCube.read(otherDataFits)

        # interpolate velocity axis. This needs to be done first.
        regridCube = otherCube.spectral_interpolate(baseCube.spectral_axis)
                
        newCube = regridCube.reproject(baseCube.header)

        if mask:
            # if greater than 0.0 set value to 1. otherwise 0.
            newdata = np.where(newCube.unitless_filled_data[:,:,:] > 0.0, 
                               1, 0)  
            finalCube = SpectralCube(newdata,newCube.wcs,
                                     mask=baseCube.mask)
        
        else:
            newdata = newCube.filled_data[:,:,:]
            finalCube = SpectralCube(newdata,newCube.wcs,mask=baseCube.mask)
            
        finalCube.write(newFits,overwrite=True)

    elif ndim == 2:

        # regrid image
        newcube = reproject_interp(otherDataFits,
                                   output_projection=baseCube.wcs.dropaxis(2),
                                   shape_out=baseCube.wcs.dropaxis(2).array_shape, 
                                   order='nearest-neighbor',
                     
                                   return_footprint=False)        

        if mask:
            # set to 1 where greater than 0.0.
            newdata = np.where(newcube > 0.0,1,0)
        else:
            newdata = newcube

        # get mask on original data
        baseMask = baseCube.get_mask_array()
        totalBaseMask = baseMask.max(axis=0)
        
        newdata[np.invert(totalBaseMask)] = np.nan
        
        # write out regridded data
        fits.writeto(newFits,newdata,baseCube.wcs.dropaxis(2).to_header(),overwrite=True)
        

    else:
        print("Number of dimensions of other data set is not 2 or 3.")
