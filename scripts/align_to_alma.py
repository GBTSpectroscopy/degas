from degas.products import makeMap
from degas.analysis_setup import regridData, smoothCube
from degas.masking import cubemask
from spectral_cube import SpectralCube

import glob
import os
import shutil
from astropy.table import Table

release = 'alma_hcn'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data',release)
coDir = os.path.join(os.environ['ANALYSISDIR'],'CO')
multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')

regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')

if not os.path.exists(regridDir):
    os.mkdir(regridDir)

hcnlist = glob.glob(os.path.join(releaseDir,'*_7m+tp_hcn_*_kms.fits'))

beam = 33.357764283516 # arcsec; EMPIRE beam
maxnchan = 131

for hcn in hcnlist:

    name = os.path.basename(hcn).split('_')[0].upper()

    print("** processing " + name + " HCN **")
    
    # process HCN -- just need to smooth here because taking as base.
    hcn_smooth = smoothCube(hcn, releaseDir, beam=beam)
    outfilename = name + "_HCN_7m+tp_33arcsec.fits"
    #outfilename = os.path.basename(hcn_smooth).replace('hcn','HCN').replace('ngc','NGC')
    shutil.copy(hcn_smooth,os.path.join(regridDir,outfilename))

    print("** processing " + name + " 12CO **")

    # process 12CO 
    co = glob.glob(os.path.join(coDir,name+'_12CO??.fits'))[0]
    if os.path.exists(co):
        co_smooth = smoothCube(co, releaseDir, beam=beam)
        co_regrid = regridData(hcn_smooth,co_smooth,regridDir)
    
    
        #calculate 12CO masks (used in stacking initial noise estimate)
        peakCut = 5.0
        lowCut = 3.0
        if name in ['IC0342']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut,lowCut=lowCut,
                     threeD=True,
                     minBeamFrac=1.5,
                     minNchan=5.0)

        elif name in ['NGC2903','NGC3521','NGC4321','NGC4535','NGC4569']:
            # use phangs parameters to create mask
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut, lowCut=lowCut,
                     minBeamFrac=2.0,
                     minNchan=3.0)

        elif name in ['NGC2146','NGC5055','NGC6946']:
            # use heracles parameters to create mask
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                      outDir=regridDir,
                      peakCut=peakCut,lowCut=lowCut,
                      minBeamFrac=1.5)

        elif name in ['NGC0337']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=3.5,lowCut=2.0,
                     minBeamFrac=1.5,
                     minNchan=3.0,
                    threeD=True)

        elif name in ['NGC3147']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=3.5,lowCut=2.0,
                     minBeamFrac=1.5,threeD=True,minNchan=5.0)   

        elif name in ['NGC3631']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut,lowCut=lowCut,
                     skipChan=[5])     

        elif name in ['NGC4030']:
            cubemask(co_regrid,
                    name+'_12CO_mask.fits',
                    outDir=regridDir,
                    peakCut=3.5,lowCut=2.0,
                    minBeamFrac=1.0,threeD=True,minNchan=3.0)

        elif name in ['NGC4501']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut,lowCut=lowCut,
                     skipChan=[15])    

        elif name in ['NGC4258']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut,lowCut=lowCut) 

        elif name in ['NGC4414']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut,
                     lowCut=lowCut)

        elif name in ['NGC4038']:
            cubemask(co_regrid,
                     name+'_12CO_mask.fits',
                     outDir=regridDir,
                     peakCut=peakCut, lowCut=lowCut)
        else:
            print("no 12CO mask created")
      
     
        # process CO products         
        maskFile = os.path.join(regridDir,name+'_12CO_mask.fits')
        makeMap(co_regrid,regridDir,maptype='peakIntensity')
        makeMap(co_regrid,regridDir,maptype='peakVelocity',maskFile=maskFile)
        makeMap(co_regrid,regridDir,maptype='moment',order=0)
        makeMap(co_regrid,regridDir,maptype='moment',order=1,maskFile=maskFile)


    print("** processing " + name + " SFR **")

    # process SFR 
    sfr = name +'_sfr_fuvw4_gauss15.fits' 
    inFile = os.path.join(multiDir,'data','sfr',sfr)
    if os.path.exists(inFile):   
        sfr_smooth = smoothCube(inFile, releaseDir,beam=beam)
        # WARNING: Could not parse unit MSUN/YR/KPC2 [spectral_cube.cube_utils]
        ## TODO -- CHECK TO MAKE SURE THIS IS OKAY.
        regridData(hcn_smooth,sfr_smooth,regridDir)

    sfr = name+'_sfr_nuvw4_gauss15.fits' 
    inFile = os.path.join(multiDir,'data','sfr',sfr)
    if os.path.exists(inFile):
        sfr_smooth = smoothCube(inFile, releaseDir,beam=beam)
        regridData(hcn_smooth,sfr_smooth,regridDir)

    sfr = name+'_sfr_fuvw3_gauss15.fits' 
    inFile = os.path.join(multiDir,'data','sfr',sfr)
    if os.path.exists(inFile):
        sfr_smooth = smoothCube(inFile, releaseDir, beam=beam)
        regridData(hcn_smooth,sfr_smooth,regridDir)

    sfr = name+'_sfr_nuvw3_gauss15.fits' 
    inFile = os.path.join(multiDir,'data','sfr',sfr)
    if os.path.exists(inFile):
        sfr_smooth = smoothCube(inFile, releaseDir, beam=beam)
        regridData(hcn_smooth,sfr_smooth,regridDir)

    # process the stellar mass data
    print("** processing " + name + " Mstar **")

    mstar = name + "_mstar_gauss15.fits"
   
    if os.path.exists(os.path.join(multiDir,'data','mstarirac1',mstar)):
        inFile = os.path.join(multiDir,'data','mstarirac1',mstar)
        mstar_smooth = smoothCube(inFile, releaseDir, beam=beam)
        regridData(hcn_smooth, mstar_smooth, regridDir)
    elif os.path.exists(os.path.join(multiDir,'data','mstarw1',mstar)):
        inFile = os.path.join(multiDir,'data','mstarw1',mstar)
        mstar_smooth = smoothCube(inFile, releaseDir, beam=beam)
        regridData(hcn_smooth, mstar_smooth, regridDir)
    else:
        print("No stellar mass map found!")

    # process the stellar mass data
    print("** processing " + name + " LTIR **")
    
    ltir = name + "_LTIR_gauss33.fits"
    inFile = os.path.join(multiDir,'data','LTIR_calc_33as',ltir)
    if os.path.exists(inFile):
        regridData(hcn_smooth, inFile,regridDir)
    else:
        print("No LTIR map found!")
                      

        
