# Purpose: make masks based on 12CO ancillary data for data reduction
# and analysis.

# Date          Programmer              Description of Changes
# ----------------------------------------------------------------------
# 4/14/2020     A.A. Kepley             Original Code
# 4/8/2021      A.A. Kepley             Updated to include phangs and new hera data


import os
from degas.masking import cubemask
from degas.products import makeMap
#import degas
from astropy.table import Table, Column
import glob
import numpy as np
from spectral_cube import SpectralCube

# calculate masks
domask = True
#domask = False

# calculate moments
generateMomentsInit = True
#generateMomentsInit = False

# set desired mask parameters
peakCut = 5.0
lowCut = 3.0

# set up the relevant directories
analysisDir = os.environ['ANALYSISDIR']
scriptDir = os.environ['SCRIPTDIR']

maskDir = os.path.join(analysisDir,'CO')
if not os.path.exists(maskDir):
    os.mkdir(maskDir)

otherDataDir = os.path.join(analysisDir,'ancillary_data')
#otherDataDir = os.environ['OTHERDATA']

# get list of galaxies in degas DR1
degas_table = Table.read(os.path.join(scriptDir,"degas_base.fits"))

# create a column for logging masks.
if 'MASK' in degas_table.colnames:
    degas_table.remove_column('MASK')
degas_table.add_column(Column(np.full_like(degas_table['NAME'],''),dtype='S25'),name='MASK')

idx_dr1 = degas_table['DR1'] == 1
#idx_dr1 = degas_table['NAME'] == 'IC0342'

# Extract list of galaxies via fancy list comprehension

# phangs
phangs_list = [os.path.basename(image).split('_')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'phangs','*10kms_gauss15.fits'))]

# heracles
heracles_list =  [os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'heracles','*gauss15_fixed_kms.fits'))]

# everyHeracles (from Adam)
extra_hera_adam_list = [os.path.basename(image).split('_')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'everyHeracles_fromadam_20210318','*10kms_gauss15.fits'))]

# everyHERACLES
extra_hera_list = [os.path.basename(image).split('_')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'co_from_andreas','*.cube.fits'))]

# bima song
bima_list = [ os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'bima_song','*gauss15_fixed.fits'))]

# OVRO
ovro_list = [ os.path.basename(image).split('.')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'temp_co','*co.cmmsk_gauss15_fixed.fits'))]

for galaxy in degas_table[idx_dr1]:

    generateMoments = generateMomentsInit

    #for IC 0342 use pipeline reduction of Jialu's 12CO
    if galaxy['NAME'] == 'IC0342':

        line='10'

        outName = galaxy['NAME']+'_12CO_mask.fits'

        #cubeFile = os.path.join(otherDataDir,'jialu',
        #                        'ic0342_regrid_12co_cube_Tmb_10kms_gauss15.fits')
        
        release = 'IR6p1'
        cubeFile = os.path.join(analysisDir,release,'IC0342_12CO_rebase7_smooth1.3_hanning1.fits')
        print('Warning: make sure you are using most current 12CO for IC0342. Currently using '+release+'.')
        if domask:
            cubemask(cubeFile,
                     outName,
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut,
                     threeD=True,
                     minBeamFrac=1.5,
                     minNchan=5.0)

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'GBT'
      
    ## use phangs first
    elif galaxy['NAME'] in phangs_list:
        
        line='21'

        outName = galaxy['NAME']+'_12CO_mask.fits'

        cubeFile = os.path.join(otherDataDir,'phangs',
                                galaxy['NAME'].lower() + '_12m+7m+tp_co21_10kms_gauss15.fits')

        if domask:
        ## works for NGC2903 and NGC3521 and NGC4569
            cubemask(cubeFile,
                     outName,
                     outDir=maskDir,
                     peakCut=peakCut, lowCut=lowCut,
                     minBeamFrac=2.0,
                     minNchan=3.0)
        
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'PHANGS'

    ## use heracles second
    elif galaxy['NAME'] in heracles_list:

        line='21'
        
        outName = galaxy['NAME']+'_12CO_mask.fits'

        cubeFile = os.path.join(otherDataDir,'heracles',
                                  galaxy['NAME']+'_heracles_gauss15_fixed_kms.fits')

        if galaxy['NAME'] == 'NGC0337':
            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=3.5,lowCut=2.0,
                         minBeamFrac=1.5,
                         minNchan=3.0,
                         threeD=True)
        else:

            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut,
                         minBeamFrac=1.5)

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'HERACLES'
   

    ## use extra HERA data from adam next 
    elif galaxy['NAME'] in extra_hera_adam_list:

        line='21'

        outName = galaxy['NAME']+'_12CO_mask.fits'
        
        cubeFile = os.path.join(otherDataDir,'everyHeracles_fromadam_20210318',
                                galaxy['NAME'].lower()+'_hera_co21_native_fixed_10kms_gauss15.fits')

        if galaxy['NAME'] == 'NGC3147':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.

            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=3.5,lowCut=2.0,
                         minBeamFrac=1.5,threeD=True,minNchan=5.0)   

        elif galaxy['NAME'] == 'NGC3631':
            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut,
                         skipChan=[5])     

        elif galaxy['NAME'] == 'NGC4030':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.

            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=3.5,lowCut=2.0,
                         minBeamFrac=1.0,threeD=True,minNchan=3.0)

        else:
            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut)
                
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'everyHERACLES'

    ## use HERA data from Andreas next (everyHERACLES)
    elif galaxy['NAME'] in extra_hera_list:

        line='21'

        outName = galaxy['NAME']+'_12CO_mask.fits'
        
        cubeFile = os.path.join(otherDataDir,'co_from_andreas',
                                galaxy['NAME'].lower()+'_hera_co21.cube_fixed_10kms_gauss15.fits')


        if galaxy['NAME'] == 'NGC3631':
            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut,
                         skipChan=[5])     

        elif galaxy['NAME'] == 'NGC4030':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.

            if domask:

                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=3.5,lowCut=2.0,
                         minBeamFrac=1.0,threeD=True,minNchan=3.0)

        elif galaxy['NAME'] == 'NGC3147':
            # this galaxy has a lot of fluffy low-level CO emission. By channel doesn't identify it as much as doing a three dimension cut for things that show up in multiple channels.

            if domask:
                
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=3.5,lowCut=2.0,
                         minBeamFrac=1.5,threeD=True,minNchan=5.0)   
                

        elif galaxy['NAME'] == 'NGC4501':
            
            if domask:

                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut,
                         skipChan=[15])    
           
        elif galaxy['NAME'] == 'NGC4535':

            if domask:

                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut,
                         minBeamFrac=1.5)    
        else:

            if domask:
            
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,lowCut=lowCut)
                
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'everyHERACLES_Andreas'

    ## use bima third
    elif galaxy['NAME'] in bima_list:

        line='10'

        outName = galaxy['NAME']+'_12CO_mask.fits'
        
        # single pointing case
        if galaxy['NAME'] == 'NGC4414':


            cubeFile = os.path.join(otherDataDir,'bima_song',
                                  galaxy['NAME']+'_bima_gauss15_fixed.fits')

            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,
                         lowCut=lowCut)

        # multiple pointing cases
        else:

            cubeFile = os.path.join(otherDataDir,'bima_song',
                                  galaxy['NAME']+'_bima_gauss15_fixed.fits')

            if domask:
                cubemask(cubeFile,
                         outName,
                         outDir=maskDir,
                         peakCut=peakCut,
                         lowCut=2.0,
                         noise3D=True)

        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'BIMASONG'

    ## use ovro next
    elif galaxy['NAME'] in ovro_list:

        line='10'

        outName = galaxy['NAME']+'_12CO_mask.fits'
        
        cubeFile = os.path.join(otherDataDir,'temp_co',
                              galaxy['NAME'].lower()+'.co.cmmsk_gauss15_fixed.fits')

        if domask:
            cubemask(cubeFile,
                     outName,
                     outDir=maskDir,
                     peakCut=peakCut,lowCut=lowCut)       
            
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'OVRO'

    ##
    elif galaxy['NAME'] == 'NGC4038':

        line='10'
        outName = galaxy['NAME']+'_12CO_mask.fits'
        
        cubeFile = os.path.join(otherDataDir,'ngc4038_from_chris',
                              'ngc_4038_4039_7m_co10_fixed_gauss15.fits')

        if domask:
            ## ALMA DATA FROM CHRIS WILSON
            cubemask(cubeFile,
                     outName,
                     outDir=maskDir,
                     peakCut=peakCut, lowCut=lowCut)
            
        degas_table['MASK'][degas_table['NAME'] == galaxy['NAME']] = 'WILSON'

    
    else:
        generateMoments=False
        print(galaxy['NAME']+" doesn't appear to have ancillary CO data.")

        
    if generateMoments:

        # copy cube over to 12CO directory. Potentially could just do
        # this via a copy command. Doing it via sectral cube so I can
        # add other features in and to make sure that the header is
        # sanitized.
        cube = SpectralCube.read(cubeFile)
        cube.write(os.path.join(maskDir,galaxy['NAME']+'_12CO'+line+'.fits'),overwrite=True)
        
        # Mom0
        makeMap(cubeFile,maskDir,
                #maskFile = os.path.join(maskDir,outName),
                baseName=galaxy['NAME']+'_12CO'+line,
                maptype='moment',order=0)
        
        # peakInt
        makeMap(cubeFile,maskDir,
                #maskFile = os.path.join(maskDir,outName),
                baseName=galaxy['NAME']+'_12CO'+line,
                maptype='peakIntensity')

        # moment 1
        makeMap(cubeFile, maskDir,
                maskFile = os.path.join(maskDir,outName),
                baseName=galaxy['NAME']+'_12CO',
                maptype='moment',order=1)

        # peak Vel 
        makeMap(cubeFile,maskDir,
                maskFile = os.path.join(maskDir,outName),
                baseName=galaxy['NAME']+'_12CO',
                maptype='peakVelocity')

        # 2D mask
        makeMap(os.path.join(maskDir,outName), maskDir,
                baseName = galaxy['NAME']+'_12CO',
                maptype='mask2D')

    
# write out degas data base table with the mask used.
degas_table.write(os.path.join(scriptDir,"degas_base.fits"),overwrite=True)
