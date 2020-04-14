# Purpose: make masks based on 12CO ancillary data for data reduction
# and analysis.

# Date          Programmer              Description of Changes
# ----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
from degas.masking import cubemask
#import degas
from astropy.table import Table
import glob

# set desired mask parameters
peak_cut = 5.0
low_cut = 3.0

# set up the relevant directories
maskDir = os.environ['MASKDIR'] 
otherDataDir = os.environ['OTHERDATA']
scriptDir = os.environ['SCRIPTDIR']

# get list of galaxies in degas DR1
degas_table = Table.read(os.path.join(scriptDir,"degas_base.fits"))

idx_dr1 = degas_table['DR1'] == 1

# Extract list of galaxies via fancy list comprehension

# heracles
heracles_list =  [os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'heracles','*gauss15.fits'))]

# bima song
bima_list = [ os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'bima_song','*gauss15_fixed.fits'))]

# OVRO
ovro_list = [ os.path.basename(image).split('.')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'temp_co','*co.cmmsk_gauss15_fixed.fits'))]

# NRO 
nro_list = [ os.path.basename(image)[0:5].replace('N','NGC') for image in glob.glob(os.path.join(otherDataDir,'temp_co','N*RD_fixed.fits'))]

for galaxy in degas_table[idx_dr1]:

    if galaxy['NAME'] in heracles_list:
        cubemask(os.path.join(otherDataDir,'heracles',
                              galaxy['NAME']+'_heracles_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)
    elif galaxy['NAME'] in bima_list:
        cubemask(os.path.join(otherDataDir,'bima_song',
                              galaxy['NAME']+'_bima_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)
    elif galaxy['NAME'] in ovro_list:
        cubemask(os.path.join(otherDataDir,'temp_co',
                              galaxy['NAME'].lower()+'.co.cmmsk_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)        
    elif galaxy['NAME'] in nro_list:
        cubemask(os.path.join(otherDataDir,'temp_co',
                              galaxy['NAME'].replace('NGC','N')+'RD_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)   
    elif galaxy['NAME'] == 'NGC4030':
        ## check on smoothing here -- this galaxy has a larger beam
        ## than 15arcsec (17.4").                
        cubemask(os.path.join(otherDataDir,'temp_co',
                              'NGC4030_12CO_RADEC.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)   
    elif galaxy['NAME'] == 'NGC4038':
        cubemask(os.path.join(otherDataDir,'temp_co',
                              'ngc4038_bigiel_carma_co_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)   
    else:
        print(galaxy['NAME']+" doesn't appear to have ancillary CO data.")
    



