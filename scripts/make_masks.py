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
peakCut = 5.0
lowCut = 3.0

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

    # for IC 0342 use Jialu's 12CO
    #if galaxy['NAME'] == 'IC0342':
        #idx = degas_table['NAME'] == 'IC0342'
        #galaxy = degas_table[idx][0]
        #cubemask(os.path.join(otherDataDir,'jialu',
        #                      'ic0342_12co_cube_Tmb_gauss15.fits'),
        #         galaxy['NAME']+'_mask.fits',
        #         outDir=maskDir,
        #         peakCut=peakCut,lowCut=lowCut,
        #         minBeamFrac=2.0,
        #         minNchan=3.0)
    ## use heracles first
    if galaxy['NAME'] in heracles_list:
        # for testing
        #idx = degas_table['NAME'] == 'NGC2903'
        #galaxy = degas_table[idx][0]
        cubemask(os.path.join(otherDataDir,'heracles',
                              galaxy['NAME']+'_heracles_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut,
                 minBeamFrac=1.5)
        ## use bima second
    elif galaxy['NAME'] in bima_list:
        cubemask(os.path.join(otherDataDir,'bima_song',
                              galaxy['NAME']+'_bima_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=1.5,minBeamFrac=0.5,minNchan=3.0,threeD=True)

        ## use ovro third
    elif galaxy['NAME'] in ovro_list:
        cubemask(os.path.join(otherDataDir,'temp_co',
                              galaxy['NAME'].lower()+'.co.cmmsk_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut)       
        
        ## use HERA data from Andreas next
    elif galaxy['NAME'] == 'NGC3631':

        cubemask(os.path.join(otherDataDir,'co_from_andreas',
                              galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut,
                 minBeamFrac=1.75)     

    elif galaxy['NAME'] == 'NGC4030':
        # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.
        cubemask(os.path.join(otherDataDir,'co_from_andreas',
                              galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=3.5,lowCut=2.0,
                 minBeamFrac=1.0,threeD=True,minNchan=3.0)    
        
    elif galaxy['NAME'] == 'NGC4501':

        cubemask(os.path.join(otherDataDir,'co_from_andreas',
                              galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut,
                 minBeamFrac=1.0,minNchan=3.0,
                 threeD=True)    


    ## use CARMA data from Frank
    elif galaxy['NAME'] == 'NGC4038':

        ## probably need to create a signal to noise image to clean this up.
        cubemask(os.path.join(otherDataDir,'temp_co',
                              'ngc4038_bigiel_carma_co_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut,
                 minNchan=3.0)   

    ## use Adam's cleaned up version of the NGC4501 cube. eliminated
    ## weirdness at edges.
    # elif galaxy['NAME'] == 'NGC4501':
    #     # The 12CO is from the NRO atlas of nearby galaxies and is
    #     # extremely noisy with edge effects. It needs some special
    #     # love and care to get a good mask.

    #     ## test for NGC4501
    #     #idx = degas_table['NAME'] == 'NGC4501'     
    #     #galaxy = degas_table[idx][0]  

    #     cubemask(os.path.join(otherDataDir,'extra_co_from_adam',
    #                           'ngc4501_nroatlas_co10_native.fits'),
    #              galaxy['NAME']+'_mask.fits',
    #              outDir=maskDir,
    #              peakCut = 3.5, # lower peak cut to get more "real" regions
    #              lowCut=3.0, # but keep same low threshold                 
    #              minBeamFrac= 2.0, # increase size of spurious region to prune
    #              minNchan=4.0,# minimum number of channels for signal
    #              skipChan = [0,1,2,3],#skip the noisy channels
    #              threeD=True) #make masks in 3d.
        
    # Last resort is NRO 45m data.
    #elif galaxy['NAME'] in nro_list:
    #    cubemask(os.path.join(otherDataDir,'temp_co',
    #                          galaxy['NAME'].replace('NGC','N')+'RD_fixed.fits'),
    #             galaxy['NAME']+'_mask.fits',
    #             outDir=maskDir,
    #             peakCut=peakCut,lowCut=lowCut)   
    
    else:
        print(galaxy['NAME']+" doesn't appear to have ancillary CO data.")
    



