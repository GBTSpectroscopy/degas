# Purpose: fix up the ancillary data so that they have appropriate
# header info and beams. This script shouldn't need to be re-run
# often, just once at the beginning of the analysis

# Date          Programmmer             Description of Changes
#----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
import glob
from degas.analysis_setup import fixBIMA
from degas.analysis_setup import fixHeracles
from degas.analysis_setup import smoothCube
import os
import glob
from radio_beam import Beam
import astropy.units as u
from spectral_cube import SpectralCube

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


## NGC4038
## -------
image = os.path.join(otherDataDir,'temp_co','ngc4038_bigiel_carma_co.fits')
cube = SpectralCube.read(image)

newBeam = Beam(beam*u.arcsec)
smoothCube = cube.convolve_to(newBeam)
smoothCube.write(image.replace('.fits','_gauss15.fits'),overwrite=True)

## NGC4030
## -------

# Beam: BMAJ=17.399999999998798 arcsec BMIN=17.399999999998798 

# So maybe want to smooth our DEGAS data to this?? this is our largest
# resolution.
