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
        shutil.copy(hcn,os.path.join(regridDir,os.path.basename(hcn)))

        print("** processing " + name + " HCO+ **")

        # process HCO+ -- smooth and regrid
        hcop = hcn.replace('hcn','hcop')
        regridData(hcn,hcop,regridDir) 
    
        print("processing " + name + " 13CO")

        # process 13CO -- smooth and regrid
        thirteenCO = hcn.replace('hcn','13co').replace('33as','27as')
        if os.path.exists(thirteenCO):
            thirteenCO_smooth = smoothCube(thirteenCO,releaseDir,beam=beam)
            regridData(hcn, thirteenCO_smooth,regridDir)

        print("processing " + name + " C18O")

        # process C18O -- smooth and regrid
        c18o = hcn.replace('hcn','c18o').replace('33as','27as')
        if os.path.exists(c18o):
            c18o_smooth = smoothCube(c18o, releaseDir, beam=beam)
            regridData(hcn,c18o_smooth,regridDir)

        print("** processing " + name + " 12CO **")

        # process 12CO 
        co = glob.glob(os.path.join(coDir,name+'_12CO??.fits'))[0]
        if os.path.exists(co):
            co_smooth = smoothCube(co, releaseDir, beam=beam)
            regridData(hcn,co_smooth,regridDir)

        # process CO products 

        ##### START HERE ######
        
        ## TODO -- need to create mom1 and peakVelocity products with
        ## appropriate resolution? Or do? I could just regrid the
        ## higher resolution products?

        # for product in ['mask','mask2D','mom0','mom1','peakInt','peakVelocity']:
        #     if ( (product == 'mom0') or (product == 'peakInt')):
        #         coproduct = glob.glob(os.path.join(coDir,name+'_12CO??_'+product+'.fits'))[0]
        #     else:
        #         coproduct = os.path.join(coDir,name+'_12CO_'+product+'.fits')
        # if os.path.exists(coproduct):
        #     if ((product is 'mask') or (product is 'mask2D')):
        #         regridData(hcn_smooth,coproduct,regridDir,mask=True)
        #     else:
        #         regridData(hcn_smooth,coproduct,regridDir)

        print("** processing " + name + " SFR **")

        # process SFR 
        ## TODO -- smooth then regrid.
        sfr = name +'_sfr_fuvw4_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):   
            sfr_smooth = smoothCube(inFile, releaseDir,beam=beam)
            regridData(hcn_smooth,sfr_smooth,regridDir)

        sfr = name+'_sfr_nuvw4_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            regridData(hcn_smooth,inFile,regridDir)

        sfr = name+'_sfr_fuvw3_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            regridData(hcn_smooth,inFile,regridDir)

        sfr = name+'_sfr_nuvw3_gauss15.fits' 
        inFile = os.path.join(multiDir,'data','sfr',sfr)
        if os.path.exists(inFile):
            regridData(hcn_smooth,inFile,regridDir)

        # process the stellar mass data
        print("** processing " + name + " Mstar **")

        mstar = name + "_mstar_gauss15.fits"
   
        if os.path.exists(os.path.join(multiDir,'data','mstarirac1',mstar)):
            regridData(hcn_smooth, os.path.join(multiDir,'data','mstarirac1',mstar), regridDir)
        elif os.path.exists(os.path.join(multiDir,'data','mstarw1',mstar)):
            regridData(hcn_smooth, os.path.join(multiDir,'data','mstarw1',mstar), regridDir)
        else:
            print("No stellar mass map found!")

        # process the stellar mass data
        print("** processing " + name + " LTIR **")
    
        ltir = name + "_LTIR_gauss15.fits"
        if os.path.exists(os.path.join(multiDir,'data','LTIR_calc',ltir)):
            regridData(hcn_smooth, os.path.join(multiDir,'data','LTIR_calc',ltir),regridDir)
        else:
            print("No LTIR map found!")
                      

        
