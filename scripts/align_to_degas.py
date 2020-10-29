from degas.analysis_setup import regridData, smoothCube

import glob
import os
import shutil

release = 'IR5'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],release)
coDir = os.path.join(os.environ['ANALYSISDIR'],'CO')
z0mgsDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','z0mgs')

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

    # process SFR -- regrid only -- already smoothed
    sfr = name+'_w4nuv_sfr.fits' # this name may change with final distribution of data to me by Adam.
    if os.path.exists(os.path.join(z0mgsDir,sfr)):
        regridData(hcn_smooth,os.path.join(z0mgsDir,sfr),regridDir)

    sfr = name+'_w4fuv_sfr.fits' # this name may change with final distribution of data to me by Adam.
    if os.path.exists(os.path.join(z0mgsDir,sfr)):
        regridData(hcn_smooth,os.path.join(z0mgsDir,sfr),regridDir)
    


    print("** processing " + name + " Mstar **")

    # process stellar mass -- regrid only -- already smoothed
    mstar = name+'_w1_stellarmass.fits' # this name may change with final distribution of data to me by Adam.
    if os.path.exists(os.path.join(z0mgsDir,mstar)):
        regridData(hcn_smooth,os.path.join(z0mgsDir,mstar),regridDir)

    w1mask = name+'_w1_gauss15_stars.fits'
    if os.path.exists(os.path.join(z0mgsDir,w1mask)):
         regridData(hcn_smooth,os.path.join(z0mgsDir,w1mask),regridDir, mask=True)

