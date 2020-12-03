from degas.analysis_setup import  calc_stellar_mass
from astropy.table import Table
from astropy.io import fits
import os

multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')
iracDir = os.path.join(multiDir,'data','irac1','convolved_with_astropy')
MtoLDir = os.path.join(multiDir,'data','mtol')

outDir = os.path.join(multiDir,'data','mstarirac1')

if not os.path.exists(outDir):
    os.mkdir(outDir)

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

## should I just do dr1 galaxies here?
for galaxy in degas:
    print("***Processing galaxy: " + galaxy['NAME'] + "***")

    iracName = galaxy['NAME'].lower() + "_irac1_gauss15.fits"
    iracFile = os.path.join(iracDir,iracName)
    MtoLName = galaxy['NAME'].lower() + "_mtol_gauss15.fits"
    MtoLFile = os.path.join(MtoLDir,MtoLName)

    # figure out if there's 
    if os.path.exists(iracFile) and os.path.exists(MtoLFile):

        outFile = os.path.join(outDir,galaxy['NAME']+"_mstar_gauss15.fits")
        calc_stellar_mass(iracFile,MtoLFile,outFile)

    else:
        print("Both IRAC and MtoL images don't exist for galaxy")
        
