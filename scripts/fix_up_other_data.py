# Purpose: fix up the ancillary data so that they have appropriate
# header info and beams. This script shouldn't need to be re-run
# often, just once at the beginning of the analysis

# Date          Programmmer             Description of Changes
#----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
import glob
import re

from degas.analysis_setup import fixBIMA, fixHeracles, fixOVRO, fixNRO, fixExtraHERA, smoothCube

from radio_beam import Beam
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.convolution import Box1DKernel
import numpy as np

# set up the relevant directories
otherDataDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data')

# common beam
beam = 15.0 #arcsec

#12CO rest frequency
rest_freq_12co = 115.27120180*u.GHz


## fix BIMA data
##----------------

print("processing BIMA data")

bima_list = glob.glob(os.path.join(otherDataDir,'bima_song','*_gauss15.fits'))

for image in bima_list:
    fixBIMA(image)

## fix heracles images 
## ---------------------

print("processing heracles data")

heracles_list = glob.glob(os.path.join(otherDataDir,'heracles','*_gauss15.fits'))
for image in heracles_list:
    fixHeracles(image)

## fix up OVRO data
## -----------------

print("processing OVRO data")

ovro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*.co.cmmsk.fits'))
for image in ovro_list:
    fixOVRO(image,beam=beam)

## fix up IC0342 12CO data from Jialu
# -----------------------------------

print("processing IC342 12CO from Jialu")

image = os.path.join(otherDataDir,'jialu','ic0342_regrid_12co_cube_Tmb_8arcsec.fits')

f = fits.open(image)

f[0].header['BUNIT'] = 'K'
f[0].header['RESTFRQ'] = float(rest_freq_12co.to(u.Hz).value)    

f.writeto(image.replace('.fits','_fixed.fits'), overwrite=True)

image =  image.replace('.fits','_fixed.fits')

# open image
cube = SpectralCube.read(image)   

# focusing this work. meta['RESTFREQ'] only puts it in as a string
#cube._header['RESTFRQ'] = float(rest_freq_12co.to(u.Hz).value)    

#cube.meta['BUNIT'] = 'K'

cube_ms = cube.with_spectral_unit(u.m / u.s) 

#, rest_value=rest_freq_12co

# chop off bad edge
subcube = cube_ms.subcube(xlo=0,xhi=134,ylo=0,yhi=150)

# add beam
beamcube = subcube.with_beam(Beam(8.0*u.arcsec))

# smooth
newBeam = Beam(beam*u.arcsec)
smoothCube = beamcube.convolve_to(newBeam)

# smooth spectrally.
smoothFactor = 3.0
spSmoothCube = smoothCube.spectral_smooth(Box1DKernel(smoothFactor))
spec_axis = spSmoothCube.spectral_axis
chan_width = spec_axis[1]-spec_axis[0] # channels are equally spaced in velocity
new_axis = np.arange(spec_axis[0].value,spec_axis[-1].value,smoothFactor*chan_width.value) * u.m/u.s

interpCube = spSmoothCube.spectral_interpolate(new_axis,
                                               suppress_smooth_warning=False)


# write out
interpCube.write(image.replace('8arcsec_fixed.fits','10kms_gauss15.fits'),
                 overwrite=True)

## fix up NRO data 
## ---------------

print("processing NRO data")

## From https://www.nro.nao.ac.jp/~nro45mrt/html/COatlas/?  Website says
## see Kuno et al. 2007, PASJ 59, 117 for details.  Beam size is given
## as 15". BUNIT is "K", so maybe only adding Beamsize.  RA and Deg
## coordinates are weird (RA---GLS and DEC--GLS)
nro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*RD.FITS'))
for image in nro_list:
    fixNRO(image)

## fix up every HERA data from Andreas
## ---------------------------------------------------------

# NGC3631 and NGC4030 are the same cubes as from adam below.

print("processing extra HERA data from Andreas")

### NEED TO ADD EFFICIENCY CORRECTIONS. ####

extra_hera_list = glob.glob(os.path.join(otherDataDir,'co_from_andreas','*hera_co21.cube.fits'))

for image in extra_hera_list:
    fixExtraHERA(image,beam=beam)


## NGC4038
## -------

print("Processing NGC4038 data from  Chris Wilson.")

# ALMA 7m from Chris Wilson

image = os.path.join(otherDataDir,'ngc4038_from_chris','ngc_4038_4039_7m_co10.fits')
cube = SpectralCube.read(image)

cube.beam
#Beam: BMAJ=15.044642430220799 arcsec BMIN=7.606865194646399 arcsec BPA=-86.90099904919 deg
 
newBeam = Beam(15.1*u.arcsec) ## beam is slightly too big in one
                              ## direction to smooth to 15.0, so
                              ## making the smoothing slightly larger.
smoothCube = cube.convolve_to(newBeam)
smoothCube.write(image.replace('.fits','_gauss15.fits'),overwrite=True)

# ## fix up NGC3631 and NGC4030 HERA data from Adam
# ## ----------------------------------------------

# # PI of proposal is Andreas Schruba. Instrument HERA, so should have properties similar ot Heracles data.

# ## ADAM ADDED EFFICIENCY CORRECTIONS TO THESE IMAGES.

# extra_hera_list = [os.path.join(otherDataDir,'extra_co_from_adam',image) for image in ['ngc3631_hera_co21_native.fits','ngc4030_hera_co21_native.fits']]

# for image in extra_hera_list:

#     f = fits.open(image)

#     f[0].header['CUNIT3'] = 'm/s'
#     newimage = image.replace('.fits','_fixed.fits')
#     f.writeto(newimage,overwrite=True)
#     f.close()

#     # open image
#     cube = SpectralCube.read(newimage)    
    
#     # smooth
#     newBeam = Beam(beam*u.arcsec)
#     smoothCube = cube.convolve_to(newBeam)
    
#     # write out
#     smoothCube.write(newimage.replace('.fits','_gauss15.fits'),
#                      overwrite=True)

