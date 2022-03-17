from degas.analysis_setup import regridData, smoothCube
from spectral_cube import SpectralCube

import glob
import os
import shutil
from astropy.table import Table

release = 'IR6p1'
#release = 'IR6p0'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],release)
coDir = os.path.join(os.environ['ANALYSISDIR'],'CO')
multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')

regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid_nosmooth')

beam=15.0
maxnchan = 131

if not os.path.exists(regridDir):
    os.mkdir(regridDir)

hcnlist = glob.glob(os.path.join(releaseDir,'*HCN*hanning1.fits'))

for hcn in hcnlist:

    name = os.path.basename(hcn).split('_')[0]

    print("** processing " + name + " HCN **")

    # process HCN -- just need to smooth here because taking as base.
    
    # first cut off any extra channels to make all cubes the same size.
    hcnCube = SpectralCube.read(hcn)
    hcnOutCube = hcnCube[0:maxnchan] 
    hcnOut = hcn.replace('.fits','_maxnchan.fits')

    hcnOutCube.write(os.path.join(releaseDir,os.path.basename(hcnOut)),overwrite=True)
    hcnOut_regrid = os.path.join(regridDir,os.path.basename(hcnOut))
    shutil.copy(hcnOut,hcnOut_regrid)

    print("** processing " + name + " HCO+ **")

    # process HCO+ -- smooth and regrid
    hcop = hcn.replace('HCN','HCOp')
    regridData(hcnOut_regrid,hcop,regridDir) # HCN and HCO+ taken at same time, so not need to check for existance first.
    
    print("processing " + name + " 13CO")

    # process 13CO -- smooth and regrid
    thirteenCO = hcn.replace('HCN','13CO')
    if os.path.exists(thirteenCO):
        regridData(hcnOut_regrid, thirteenCO,regridDir)

    print("processing " + name + " C18O")

    # process C18O -- smooth and regrid
    c18o = hcn.replace('HCN','C18O')
    if os.path.exists(c18o):
        regridData(hcnOut_regrid,c18o,regridDir)

    print("** processing " + name + " 12CO **")

    # process 12CO -- regrid only -- already smoothed
    co = glob.glob(os.path.join(coDir,name+'_12CO??.fits'))[0]
    if os.path.exists(co):
        regridData(hcnOut_regrid,co,regridDir)

    # process CO products -- regrid only -- already smoothed.
    for product in ['mask','mask2D','mom0','mom1','peakInt','peakVelocity']:
        if ( (product == 'mom0') or (product == 'peakInt')):
            coproduct = glob.glob(os.path.join(coDir,name+'_12CO??_'+product+'.fits'))[0]
        else:
            coproduct = os.path.join(coDir,name+'_12CO_'+product+'.fits')
        if os.path.exists(coproduct):
            if ((product is 'mask') or (product is 'mask2D')):
                regridData(hcnOut_regrid,coproduct,regridDir,mask=True)
            else:
                regridData(hcnOut_regrid,coproduct,regridDir)

                      

        
