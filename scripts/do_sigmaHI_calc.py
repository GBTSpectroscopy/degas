from degas.analysis_setup import calc_sigmaHI
from astropy.table import Table, Column
import os
import glob


hiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','hi_from_jiayi')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
hi_db = Table.read(os.path.join(hiDir,'DEGAS-HI-summary.fits'))


outDir = os.path.join(hiDir,'sigmaHI_15arcsec')
if not os.path.exists(outDir):
    os.mkdir(outDir)

myfiles = glob.glob(os.path.join(hiDir,'15arcsec','*.fits'))
for fitsimage in myfiles:
    
    galaxy = os.path.basename(fitsimage).split('_')[0]
    incl = degas['INCL_DEG'][degas['NAME'] == galaxy]

    calc_sigmaHI(fitsimage,
                 outDir,
                 incl=incl)
    
outDir = os.path.join(hiDir,'sigmaHI_33arcsec')
if not os.path.exists(outDir):
    os.mkdir(outDir)

myfiles = glob.glob(os.path.join(hiDir,'33arcsec','*.fits'))
for fitsimage in myfiles:
    
    galaxy = os.path.basename(fitsimage).split('_')[0]
    incl = degas['INCL_DEG'][degas['NAME'] == galaxy]

    calc_sigmaHI(fitsimage,
                 outDir,
                 incl=incl)
    
