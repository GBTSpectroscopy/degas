# Purpose: Rescale the 12CO21 data to the most likely 12CO10
# values. This work draws heavily on Leroy+ 2021

# Date          Programmer      Description of Changes
#----------------------------------------------------------------------
# 10/20/21      A.A. Kepley     Original Code


import os

from astropy.table import Table, Column
import glob
import numpy as np
from spectral_cube import SpectralCube
from astropy.io import fits
from degas.analysis_setup import simpleR21scale

release = 'IR6p1'

# set up relevant directories
analysisDir = os.environ['ANALYSISDIR']
scriptDir = os.environ['SCRIPTDIR']
dataDir = os.path.join(analysisDir,release+'_regrid')

# read in degas table
degas_table = Table.read(os.path.join(scriptDir,"degas_base.fits"))
idx_dr1 = degas_table['DR1'] == 1

# set fiducial value
r21 = 0.65 # mean R21 value from Leroy+ 2021
r21_ref = 'Leroy+2021' 

for galaxy in degas_table[idx_dr1]:

    if ( (galaxy['MASK'] == 'HERACLES') |
         (galaxy['MASK'] == 'everyHERACLES') |
         (galaxy['MASK'] == 'PHANGS') | 
         (galaxy['MASK'] == 'everyHERACLES_Andreas') ):
        
        # Do the simple scaling by a single value
        cube = os.path.join(dataDir,galaxy['NAME']+"_12CO21_regrid.fits")
        mom0 = os.path.join(dataDir,galaxy['NAME']+"_12CO21_mom0_regrid.fits")
        peakInt = os.path.join(dataDir,galaxy['NAME']+"_12CO21_peakInt_regrid.fits")
        sigmaSFR = os.path.join(dataDir,galaxy['NAME']+"_sfr_fuvw4_gauss15_regrid.fits")

        simpleR21scale(cube,r21,r21_ref=r21_ref)
        simpleR21scale(mom0,r21,r21_ref=r21_ref)
        simpleR21scale(peakInt,r21,r21_ref=r21_ref)

        # now do the more complex scaling by sigma_SFR
        ## TODO 
        ##      -- get re from Adam
        ##      -- confirm what sfr to use
        ##      -- decide whether to use galaxy specific r21
        #sfrR21scale(cube, sigmaSFR, galaxy=galaxy['NAME'],
        #            re=re, sfr=sfr, r21=r21, r21_ref=r21_ref)
