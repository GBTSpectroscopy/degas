from degas.products import makeMap
from degas.analysis_setup import regridData, smoothCube
from spectral_cube import SpectralCube

import glob
import os
import shutil
from astropy.table import Table

release = 'empire'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','empire','EMPIRE_cubes')
coDir = os.path.join(os.environ['ANALYSISDIR'],'CO')
multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')

regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')

# TODO -- determine automatically rather than manually.
overlaplist = ['NGC2903','NGC4321','NGC5055','NGC6946']

if not os.path.exists(regridDir):
    os.mkdir(regridDir)

hcnlist = glob.glob(os.path.join(releaseDir,'*_hcn_*_fixed_kms.fits'))

for hcn in hcnlist:

    name = os.path.basename(hcn).split('_')[1].upper()

    if name in overlaplist:

        print("** processing " + name + " HCN **")

        # process HCN -- just need to smooth here because taking as base.
        hcn_cube = SpectralCube.read(hcn)
        beam = hcn_cube.beam.major.to('arcsec').value # assumes bmaj=bmin
        hcnOut = os.path.basename(hcn).replace('EMPIRE_','').replace(name.lower(),name).replace('hcn','HCN')
        shutil.copy(hcn, os.path.join(regridDir,hcnOut))

        print("** processing " + name + " HCO+ **")

        # process HCO+ -- regrid -- already smoothed
        hcop = hcn.replace('hcn','hcop')
        hcop_regrid = regridData(hcn,hcop,regridDir) 
        hcopOut = hcop_regrid.replace('EMPIRE_','').replace(name.lower(),name).replace('hcop','HCOp')
        os.rename(hcop_regrid,hcopOut)
        
        print("processing " + name + " 13CO")

        # process 13CO -- smooth and regrid
        thirteenCO = hcn.replace('hcn','13co').replace('33as','27as')
        if os.path.exists(thirteenCO):
            thirteenCO_smooth = smoothCube(thirteenCO,releaseDir,beam=beam)
            thirteenCO_regrid = regridData(hcn, thirteenCO_smooth, regridDir)
            thirteenCO_out = thirteenCO_regrid.replace('EMPIRE_','').replace(name.lower(),name).replace('13co','13CO').replace('27as','33as')
            os.rename(thirteenCO_regrid,thirteenCO_out)

        print("processing " + name + " C18O")

        # process C18O -- smooth and regrid
        c18o = hcn.replace('hcn','c18o').replace('33as','27as')
        if os.path.exists(c18o):
            c18o_smooth = smoothCube(c18o, releaseDir, beam=beam)
            c18o_regrid = regridData(hcn,c18o_smooth,regridDir)
            c18o_out = c18o_regrid.replace('EMPIRE_','').replace(name.lower(),name).replace('c18o','C18O').replace('27as','33as')
            os.rename(c18o_regrid, c18o_out)

        print("** processing " + name + " 12CO **")

        # process 12CO 
        co = glob.glob(os.path.join(coDir,name+'_12CO??.fits'))[0]
        if os.path.exists(co):
            co_smooth = smoothCube(co, releaseDir, beam=beam)
            co_regrid = regridData(hcn,co_smooth,regridDir)

        # process CO products         
        makeMap(co_regrid,regridDir,maptype='peakIntensity')
        makeMap(co_regrid,regridDir,maptype='peakVelocity')
        makeMap(co_regrid,regridDir,maptype='moment',order=0)
        makeMap(co_regrid,regridDir,maptype='moment',order=1)

        print("** processing " + name + " SFR **")

        # process SFR 
        sfr = name +'_sfr_fuvw4_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):   
            sfr_smooth = smoothCube(inFile, releaseDir,beam=beam)
            # WARNING: Could not parse unit MSUN/YR/KPC2 [spectral_cube.cube_utils]
            ## TODO -- CHECK TO MAKE SURE THIS IS OKAY.
            regridData(hcn,sfr_smooth,regridDir)

        sfr = name+'_sfr_nuvw4_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            sfr_smooth = smoothCube(inFile, releaseDir,beam=beam)
            regridData(hcn,inFile,regridDir)

        sfr = name+'_sfr_fuvw3_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            sfr_smooth = smoothCube(inFile, releaseDir, beam=beam)
            regridData(hcn,inFile,regridDir)

        sfr = name+'_sfr_nuvw3_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            sfr_smooth = smoothCube(inFile, releaseDir, beam=beam)
            regridData(hcn,inFile,regridDir)

        # process the stellar mass data
        print("** processing " + name + " Mstar **")

        mstar = name + "_mstar_gauss15.fits"
   
        if os.path.exists(os.path.join(multiDir,'data','mstarirac1',mstar)):
            inFile = os.path.join(multiDir,'data','mstarirac1',mstar)
            mstar_smooth = smoothCube(inFile, releaseDir, beam=beam)
            regridData(hcn, mstar_smooth, regridDir)
        elif os.path.exists(os.path.join(multiDir,'data','mstarw1',mstar)):
            inFile = os.path.join(multiDir,'data','mstarw1',mstar)
            mstar_smooth = smoothCube(inFile, releaseDir, beam=beam)
            regridData(hcn, mstar_smooth, regridDir)
        else:
            print("No stellar mass map found!")

        # process the stellar mass data
        print("** processing " + name + " LTIR **")
    
        ltir = name + "_LTIR_gauss15.fits"
        inFile = os.path.join(multiDir,'data','LTIR_calc',ltir)
        if os.path.exists(inFile):
            ltir_smooth = smoothCube(inFile, releaseDir, beam=beam)
            regridData(hcn, inFile,regridDir)
        else:
            print("No LTIR map found!")
                      

        
