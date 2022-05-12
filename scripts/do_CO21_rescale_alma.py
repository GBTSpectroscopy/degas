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
from degas.analysis_setup import simpleR21scale, sfrR21scale
import astropy.units as u
import ipdb

release = 'alma_hcn'

# set up relevant directories
analysisDir = os.environ['ANALYSISDIR']
scriptDir = os.environ['SCRIPTDIR']
dataDir = os.path.join(analysisDir,release+'_regrid')

# read in degas table
degas = Table.read(os.path.join(scriptDir,"degas_base.fits"))
idx_dr1 = (degas['DR1'] == 1) & ((degas['NAME'] == 'NGC2903') | (degas['NAME'] == 'NGC4321'))

# set fiducial value
r21 = 0.65 # mean R21 value from Leroy+ 2021
r21_ref = 'Leroy+2021' 

# read in R21 table
r21_table = Table.read(os.path.join(analysisDir,'database','int_corat_tab.ecsv'))

for galaxy in degas[idx_dr1]:

    if ( (galaxy['MASK'] == 'HERACLES') |
         (galaxy['MASK'] == 'everyHERACLES') |
         (galaxy['MASK'] == 'PHANGS') | 
         (galaxy['MASK'] == 'everyHERACLES_Andreas') ):
        
        # Do the simple scaling by a single value
        cube = os.path.join(dataDir,galaxy['NAME']+"_12CO21_smooth_regrid.fits")
        mom0 = os.path.join(dataDir,galaxy['NAME']+"_12CO21_smooth_regrid_mom0.fits")
        peakInt = os.path.join(dataDir,galaxy['NAME']+"_12CO21_smooth_regrid_peakInt.fits")
        sigmaSFR = os.path.join(dataDir,galaxy['NAME']+"_sfr_fuvw4_gauss15_smooth_regrid.fits")

        simpleR21scale(cube,r21,r21_ref=r21_ref)
        simpleR21scale(mom0,r21,r21_ref=r21_ref)
        simpleR21scale(peakInt,r21,r21_ref=r21_ref)

        # pick R21
        idx = (r21_table['GALAXY'] == galaxy['NAME'].lower()) & (r21_table['RATNAME'] == 'R21') 
        galaxy_ratios = r21_table[idx]

        if 'PHANGSCOMING' in galaxy_ratios['SURVEY_PAIR']:
            idx = galaxy_ratios['SURVEY_PAIR'] == 'PHANGSCOMING'
            r21_gal = 10**galaxy_ratios[idx]['LOGRAT'][0]
        elif 'PHANGSNROATLAS' in galaxy_ratios['SURVEY_PAIR']:
            idx = galaxy_ratios['SURVEY_PAIR'] == 'PHANGSNROATLAS'
            r21_gal = 10**galaxy_ratios[idx]['LOGRAT'][0]
        elif 'HERACOMING' in galaxy_ratios['SURVEY_PAIR']:
            idx = galaxy_ratios['SURVEY_PAIR'] == 'HERACOMING'
            r21_gal = 10**galaxy_ratios[idx]['LOGRAT'][0]
        elif 'HERANROATLAS' in galaxy_ratios['SURVEY_PAIR']:
            idx = galaxy_ratios['SURVEY_PAIR'] == 'HERANROATLAS'
            r21_gal = 10**galaxy_ratios[idx]['LOGRAT'][0]
        else:
            print('No galaxy-specific R21 value found for '+galaxy['NAME']+ ". Using fiducial value.")
            r21_gal = r21

        # convert Re to to kpc
        # The input map is in Mstar/yr/kpc^2, so I think I want the output unit to be Mstar/yr/kpc^2.
        re_kpc = (( galaxy['RE_ARCSEC']*u.arcsec) * (galaxy['DIST_MPC']*u.Mpc)).to(u.kpc,equivalencies=u.dimensionless_angles())


        # now do more complex spatial scaling
        sfrR21scale(cube,sigmaSFR, 
                    galaxy=galaxy['NAME'],
                    re = re_kpc.value, 
                    sfr = 10**galaxy['LOGSFR'],
                    r21 = r21_gal, 
                    r21_ref = r21_ref)

        sfrR21scale(mom0,sigmaSFR, 
                    galaxy=galaxy['NAME'],
                    re = re_kpc.value, 
                    sfr = 10**galaxy['LOGSFR'],
                    r21 = r21_gal, 
                    r21_ref = r21_ref)


        sfrR21scale(peakInt,sigmaSFR, 
                    galaxy=galaxy['NAME'],
                    re = re_kpc.value, 
                    sfr = 10**galaxy['LOGSFR'],
                    r21 = r21_gal, 
                    r21_ref = r21_ref)

