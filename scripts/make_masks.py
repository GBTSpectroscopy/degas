# Purpose: make masks based on 12CO ancillary data for data reduction
# and analysis.

# Date          Programmer              Description of Changes
# ----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code

import os
from degas.masking import cubemask
#import degas
from astropy.table import Table, Column
import glob
import numpy as np

# set desired mask parameters
peakCut = 5.0
lowCut = 3.0

# set up the relevant directories
maskDir = os.environ['MASKDIR'] 
otherDataDir = os.environ['OTHERDATA']
scriptDir = os.environ['SCRIPTDIR']

# get list of galaxies in degas DR1
degas_table = Table.read(os.path.join(scriptDir,"degas_base.fits"))

# create a column for logging masks.
if 'MASK' in degas_table.colnames:
    degas_table.remove_column('MASK')
degas_table.add_column(Column(np.full_like(degas_table['NAME'],''),dtype='S15'),name='MASK')

#idx_dr1 = degas_table['DR1'] == 1
# for testing
idx_dr1 = degas_table['NAME'] == 'IC0342'

# Extract list of galaxies via fancy list comprehension

# heracles
heracles_list =  [os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'heracles','*gauss15.fits'))]

# bima song
bima_list = [ os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'bima_song','*gauss15_fixed.fits'))]

# everyHERACLES
extra_hera_list = [os.path.basename(image).split('_')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'co_from_andreas','*.cube.fits'))]

# OVRO
ovro_list = [ os.path.basename(image).split('.')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'temp_co','*co.cmmsk_gauss15_fixed.fits'))]

for galaxy in degas_table[idx_dr1]:

    #for IC 0342 use Jialu's 12CO
    if galaxy['NAME'] == 'IC0342':

        cubemask(os.path.join(otherDataDir,'jialu',
                              'ic0342_regrid_12co_cube_Tmb_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut,
                 threeD=True,
                 minBeamFrac=1.5,
                 minNchan=5.0)
    ## use heracles first
    elif galaxy['NAME'] in heracles_list:

        if galaxy['NAME'] == 'NGC0337':
            cubemask(os.path.join(otherDataDir,'heracles',
                                  galaxy['NAME']+'_heracles_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=3.5,lowCut=2.0,
                     minBeamFrac=1.5,
                     minNchan=3.0,
                     threeD=True)
        else:

            cubemask(os.path.join(otherDataDir,'heracles',
                                  galaxy['NAME']+'_heracles_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut,
                     minBeamFrac=1.5)

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'HERACLES'
   
    ## use HERA data from Andreas next (everyHERACLES)
    elif galaxy['NAME'] in extra_hera_list:
    
        if galaxy['NAME'] == 'NGC3631':

            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut,
                     skipChan=[5])     

        elif galaxy['NAME'] == 'NGC4030':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.
            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=3.5,lowCut=2.0,
                     minBeamFrac=1.0,threeD=True,minNchan=3.0)

        elif galaxy['NAME'] == 'NGC3147':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.
            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=3.5,lowCut=2.0,
                     minBeamFrac=1.5,threeD=True,minNchan=5.0)   

        elif galaxy['NAME'] == 'NGC4501':
            
            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut,
                     skipChan=[15])    
           
        elif galaxy['NAME'] == 'NGC4535':
            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut,
                     minBeamFrac=1.5)    
        else:
            cubemask(os.path.join(otherDataDir,'co_from_andreas',
                                  galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut)

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'everyHERACLES'

    ## use bima third
    elif galaxy['NAME'] in bima_list:
        
        # single pointing case
        if galaxy['NAME'] == 'NGC4414':
            cubemask(os.path.join(otherDataDir,'bima_song',
                                  galaxy['NAME']+'_bima_gauss15_fixed.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,
                     lowCut=lowCut)
        # multiple pointing cases
        else:
            cubemask(os.path.join(otherDataDir,'bima_song',
                                  galaxy['NAME']+'_bima_gauss15_fixed.fits'),
                     galaxy['NAME']+'_mask.fits',
                     outDir=maskDir,
                     peakCut=peakCut,
                     lowCut=2.0,
                     noise3D=True)
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'BIMASONG'

    ## use ovro next
    elif galaxy['NAME'] in ovro_list:
        cubemask(os.path.join(otherDataDir,'temp_co',
                              galaxy['NAME'].lower()+'.co.cmmsk_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut,lowCut=lowCut)       

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] == 'OVRO'

    ##
    elif galaxy['NAME'] == 'NGC4038':

        ## ALMA DATA FROM CHRIS WILSON
        cubemask(os.path.join(otherDataDir,'ngc4038_from_chris',
                              'ngc_4038_4039_7m_co10_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peakCut=peakCut, lowCut=lowCut)
        
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] == 'WILSON'

    
    else:
        print(galaxy['NAME']+" doesn't appear to have ancillary CO data.")
    
# write out degas data base table with the mask used.
degas_table.write(os.path.join(scriptDir,"degas_base.fits"),overwrite=True)
