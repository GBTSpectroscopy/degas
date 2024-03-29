# Purpose: fix up the ancillary data so that they have appropriate
# header info and beams. This script shouldn't need to be re-run
# often, just once at the beginning of the analysis

# Date          Programmmer             Description of Changes
#----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
import glob
import re

import degas.analysis_setup as analysis_setup

from radio_beam import Beam
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.convolution import Box1DKernel
import numpy as np

# set up the relevant directories
otherDataDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data')

# common beam
common_beam = 15.0 #arcsec


## fix BIMA data
##----------------

print("processing BIMA data")

bima_list = glob.glob(os.path.join(otherDataDir,'bima_song','*_gauss15.fits'))

for image in bima_list:
    analysis_setup.fixBIMA(image)

## fix heracles images 
## ---------------------

print("processing heracles data")

heracles_list = glob.glob(os.path.join(otherDataDir,'heracles','*_gauss15.fits'))
for image in heracles_list:
    analysis_setup.fixHeracles(image)

## fix up OVRO data
## -----------------

# print("processing OVRO data")

# ovro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*.co.cmmsk.fits'))
# for image in ovro_list:
#     analysis_setup.fixOVRO(image,beam=common_beam)

## fix up IC0342 12CO data from Jialu
# -----------------------------------

print("processing IC342 12CO from Jialu")

image = os.path.join(otherDataDir,'jialu','ic0342_regrid_12co_cube_Tmb_8arcsec.fits')

analysis_setup.fixJialu(image,beam=common_beam)

## fix up NRO data 
## ---------------

# print("processing NRO data")

# ## From https://www.nro.nao.ac.jp/~nro45mrt/html/COatlas/?  Website says
# ## see Kuno et al. 2007, PASJ 59, 117 for details.  Beam size is given
# ## as 15". BUNIT is "K", so maybe only adding Beamsize.  RA and Deg
# ## coordinates are weird (RA---GLS and DEC--GLS)
# nro_list = glob.glob(os.path.join(otherDataDir,'temp_co','*RD.FITS'))
# for image in nro_list:
#     analysis_setup.fixNRO(image)

## fix up every HERA data from Andreas
## ---------------------------------------------------------

# NGC3631 and NGC4030 are the same cubes as from adam below.

## this data doesn't include an eta_mb efficiency correction. The correction is applied below vis fixExtraHERA

print("processing extra HERA data from Andreas")

extra_hera_list = glob.glob(os.path.join(otherDataDir,'co_from_andreas','*hera_co21.cube.fits'))

for image in extra_hera_list:
    analysis_setup.fixExtraHERA(image,beam=common_beam)

## fix up extra HERA data from Adam 
## --------------------------------------------------

### Includes eta_mb efficiency correction

extra_hera_adam_list = glob.glob(os.path.join(otherDataDir,'everyHeracles_fromadam_20210318','*hera_co21_native.fits'))

for image in extra_hera_adam_list:
    analysis_setup.fixExtraHERAfromAdam(image,beam=common_beam)


## NGC4038
## -------

print("Processing NGC4038 data from  Chris Wilson.")

# ALMA 7m from Chris Wilson

image = os.path.join(otherDataDir,'ngc4038_from_chris','ngc_4038_4039_7m_co10.fits')

## TODO -- move the code below to a function in analysis_setup.py

## getting rid of some suspicious headers
hdu = fits.open(image)
hdr = hdu[0].header
data = hdu[0].data

hdr.remove('ALTRPIX')
hdr.remove('ALTRVAL')
hdr.remove('OBSGEO-X')
hdr.remove('OBSGEO-Y')
hdr.remove('OBSGEO-Z')
#hdr.remove('MJD-OBS')
hdr.remove('DATE-OBS')
hdr.remove('DISTANCE')

newimage  = os.path.join(otherDataDir,'ngc4038_from_chris','ngc_4038_4039_7m_co10_fixed.fits')
fits.writeto(newimage,data,hdr,overwrite=True)

# regridding
cube = SpectralCube.read(newimage)

cube_kms =  cube.with_spectral_unit(u.km / u.s)

cube_kms.beam
#Beam: BMAJ=15.044642430220799 arcsec BMIN=7.606865194646399 arcsec BPA=-86.90099904919 deg
 
newBeam = Beam(15.1*u.arcsec) ## beam is slightly too big in one
                              ## direction to smooth to 15.0, so
                              ## making the smoothing slightly larger.
smoothCube = cube_kms.convolve_to(newBeam)

smoothCubeinK = smoothCube.to(u.K)

smoothCubeinK.write(newimage.replace('.fits','_gauss15.fits'),overwrite=True)

## PHANGS data
## -----------

print("Processing PHANGS data")

phangsList = glob.glob(os.path.join(otherDataDir,'phangs',"*_co21_7p5as.fits"))
for image in phangsList:
    analysis_setup.fixPhangs(image,beam=common_beam)
    

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

