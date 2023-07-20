from degas.analysis_setup import smoothCube
from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import os
import glob
import shutil

hiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','hi_from_jiayi')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
hi_db = Table.read(os.path.join(hiDir,'DEGAS-HI-summary.fits'))

for galaxy in degas:
    print("***Processing galaxy: " + galaxy['NAME'] + "***")
    
    if not hi_db['any'][hi_db['NAME'] == galaxy['NAME']]:
        print("No HI data available for this galaxy. Skipping.")
        continue
    
    if hi_db['THINGS'][hi_db['NAME'] == galaxy['NAME']]:
        print('galaxy in things')
        mom0_path = os.path.join(hiDir,'archival',galaxy['NAME'] + '_THINGS_21cm_strictmask_mom0.fits')
    elif hi_db['PHANGSVLA'][hi_db['NAME'] == galaxy['NAME']]:
        print('galaxy in phangsvla')
        mom0_path = os.path.join(hiDir,'PHANGSVLA',galaxy['NAME'] + '_PHANGSVLA_21cm_strictmask_mom0.fits')
    elif hi_db['VLAHERACLES'][hi_db['NAME'] == galaxy['NAME']]:
            print('galaxy in vlaheracles')
            mom0_path = os.path.join(hiDir,'archival',galaxy['NAME'] + '_VLAHERACLES_21cm_strictmask_mom0.fits')
    elif hi_db['HALOGAS'][hi_db['NAME'] == galaxy['NAME']]:
            print('galaxy in halogas')
            mom0_path = os.path.join(hiDir,'archival',galaxy['NAME'] + '_HALOGAS_21cm_strictmask_mom0.fits')
    elif hi_db['EVERYTHINGS'][hi_db['NAME'] == galaxy['NAME']]:
        print("Galaxy in everythings")
        mom0_path = os.path.join(hiDir,'EVERYTHINGS',galaxy['NAME'] + '_EVERYTHINGS_21cm_strictmask_mom0.fits')
    else:
        print('galaxy not found')

    hi_mom0 = fits.open(mom0_path)
    beam_arcsec = hi_mom0[0].header['BMAJ']*(206265.0/60.0)
    hi_mom0.close()

    print('beam: ' + str(beam_arcsec))
    analysis_beam = 15.0
    if beam_arcsec <= analysis_beam:
        outDir = os.path.join(hiDir,'15arcsec')
        print(outDir)
        if not os.path.exists(outDir):
            os.mkdir(outDir)

      
        if ((galaxy['NAME'] == 'NGC4569') or (galaxy['NAME'] == 'NGC3521')):            
            # beam here is 14.9arcsec and 14.85, which is close enough to 15 to not make a difference. Just copying over.
            shutil.copy(mom0_path,os.path.join(outDir,os.path.basename(mom0_path).replace('.fits','_smooth.fits')))
        else:
            smoothCube(mom0_path,outDir,
                       beam=analysis_beam)
            
    empire_beam = 33.357764283516
    if beam_arcsec <= empire_beam:
        outDir = os.path.join(hiDir,'33arcsec')
        print(outDir)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        smoothCube(mom0_path,outDir,
                   beam=empire_beam)



    


    
    
    
                            

