from degas.analysis_setup import regridData, smoothCube

import glob
import os
import shutil
from astropy.table import Table

release = 'IR5'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],release)
coDir = os.path.join(os.environ['ANALYSISDIR'],'CO')
multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')

regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')

beam=15.0

if not os.path.exists(regridDir):
    os.mkdir(regridDir)

hcnlist = glob.glob(os.path.join(releaseDir,'*HCN*hanning1.fits'))

for hcn in hcnlist:

    name = os.path.basename(hcn).split('_')[0]

    print("** processing " + name + " HCN **")

    # process HCN -- just need to smooth here because taking as base.
    hcn_smooth = smoothCube(hcn, releaseDir, beam=beam)
    shutil.copy(hcn_smooth,os.path.join(regridDir,os.path.basename(hcn_smooth)))

    print("** processing " + name + " HCO+ **")

    # process HCO+ -- smooth and regrid
    hcop = hcn.replace('HCN','HCOp')
    hcop_smooth = smoothCube(hcop,releaseDir,beam=beam)
    regridData(hcn_smooth,hcop_smooth,regridDir) # HCN and HCO+ taken at same time, so not need to check for existance first.
    
    print("processing " + name + " 13CO")

    # process 13CO -- smooth and regrid
    thirteenCO = hcn.replace('HCN','13CO')
    if os.path.exists(thirteenCO):
        thirteenCO_smooth = smoothCube(thirteenCO,releaseDir,beam=beam)
        regridData(hcn_smooth, thirteenCO_smooth,regridDir)

    print("processing " + name + " C18O")

    # process C18O -- smooth and regrid
    c18o = hcn.replace('HCN','C18O')
    if os.path.exists(c18o):
        c18o_smooth = smoothCube(c18o, releaseDir, beam=beam)
        regridData(hcn_smooth,c18o_smooth,regridDir)

    print("** processing " + name + " 12CO **")

    # process 12CO -- regrid only -- already smoothed
    co = os.path.join(coDir,name+'_12CO.fits')
    if os.path.exists(co):
        regridData(hcn_smooth,co,regridDir)

    # process CO products -- regrid only -- already smoothed.
    for product in ['mask','mask2D','mom0','mom1','peakInt','peakVelocity']:
        coproduct = os.path.join(coDir,name+'_12CO_'+product+'.fits')
        if os.path.exists(coproduct):
            if ((product is 'mask') or (product is 'mask2D')):
                regridData(hcn_smooth,coproduct,regridDir,mask=True)
            else:
                regridData(hcn_smooth,coproduct,regridDir)

    print("** processing " + name + " SFR **")

    # process SFR -- regrid only -- already smoothed to 15arcsec
    sfr = name +'_sfr_fuvw4_gauss15.fits' 
    inFile = os.path.join(multiDir,'data','sfr',sfr)
    if os.path.exists(inFile):
        regridData(hcn_smooth,inFile,regridDir)

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

   
