# Purpose: fix up the ancillary data so that they have appropriate
# header info and beams. This script shouldn't need to be re-run
# often, just once at the beginning of the analysis

# Date          Programmmer             Description of Changes
#----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
import glob
import re

from degas.analysis_setup import fixBIMA
from degas.analysis_setup import fixHeracles
from degas.analysis_setup import smoothCube

from radio_beam import Beam
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.convolution import Box1DKernel
import numpy as np

# set up the relevant directories
otherDataDir = os.environ['OTHERDATA']

# common beam
beam = 15.0 #arcsec

## fix BIMA data
##----------------

bima_list = glob.glob(os.path.join(otherDataDir,'bima_song','*_gauss15.fits'))

for image in bima_list:
    fixBIMA(image)

## fix heracles images 
## ---------------------

## -- what needs to be done here?? They look okay for masks, but
## probably also need want to check the regridding to make sure
## everything works out properly.


## fix up OVRO data
## -----------------

ovro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*.co.cmmsk.fits'))
for image in ovro_list:

    # open image
    cube = SpectralCube.read(image)    
    # convert to K
    kcube = cube.to(u.K)

    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = kcube.convolve_to(newBeam)

    # write out
    smoothCube.write(image.replace('.fits','_gauss15_fixed.fits'),
                     overwrite=True)

## fix up IC0342 12CO data from Jialu
# -----------------------------------

image = os.path.join(otherDataDir,'jialu','ic0342_12co_cube_Tmb_8arcsec.fits')

# open image
cube = SpectralCube.read(image)    
beamcube = cube.with_beam(Beam(8.0*u.arcsec))

# smooth
newBeam = Beam(beam*u.arcsec)
smoothCube = beamcube.convolve_to(newBeam)

# write out
smoothCube.write(image.replace('8arcsec.fits','gauss15.fits'),
                 overwrite=True)

## fix up NRO data 
## ---------------

## From https://www.nro.nao.ac.jp/~nro45mrt/html/COatlas/?  Website says
## see Kuno et al. 2007, PASJ 59, 117 for details.  Beam size is given
## as 15". BUNIT is "K", so maybe only adding Beamsize.  RA and Deg
## coordinates are weird (RA---GLS and DEC--GLS)
nro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*RD.FITS'))
for image in nro_list:
    cube = SpectralCube.read(image)
    # add beam
    beam_cube = cube.with_beam(Beam(15.0*u.arcsec))
    beam_cube.write(image.replace('.FITS','_fixed.fits'),overwrite=True)



## fix up NGC3631, NGC4030, and NGC4501 HERA data from Andreas
## ---------------------------------------------------------

# NGC3631 and NGC4030 are the same cubes as from adam below.

### NEED TO ADD EFFICIENCY CORRECTIONS. ####

extra_hera_list = [os.path.join(otherDataDir,'co_from_andreas',image) for image in ['ngc3631_hera_co21.cube.fits','ngc4030_hera_co21.cube.fits','ngc4501_hera_co21.cube.fits']]

for image in extra_hera_list:

    #image = os.path.join(otherDataDir,'co_from_andreas','ngc3631_hera_co21.cube.fits')

    f = fits.open(image)
    f[0].header['CUNIT3'] = 'm/s'
    newimage = image.replace('.fits','_fixed.fits')
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
    

## fix up NGC3631 and NGC4030 HERA data from Adam
## ----------------------------------------------

# PI of proposal is Andreas Schruba. Instrument HERA, so should have properties similar ot Heracles data.

## ADAM ADDED EFFICIENCY CORRECTIONS TO THESE IMAGES.

extra_hera_list = [os.path.join(otherDataDir,'extra_co_from_adam',image) for image in ['ngc3631_hera_co21_native.fits','ngc4030_hera_co21_native.fits']]

for image in extra_hera_list:

    f = fits.open(image)

    f[0].header['CUNIT3'] = 'm/s'
    newimage = image.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)
    f.close()

    # open image
    cube = SpectralCube.read(newimage)    
    
    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = cube.convolve_to(newBeam)
    
    # write out
    smoothCube.write(newimage.replace('.fits','_gauss15.fits'),
                     overwrite=True)

## NGC4038
## -------

image = os.path.join(otherDataDir,'temp_co','ngc4038_bigiel_carma_co.fits')
cube = SpectralCube.read(image)

newBeam = Beam(beam*u.arcsec)
smoothCube = cube.convolve_to(newBeam)
smoothCube.write(image.replace('.fits','_gauss15.fits'),overwrite=True)

