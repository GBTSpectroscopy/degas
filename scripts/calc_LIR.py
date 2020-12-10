from degas.analysis_setup import calc_LTIR
from astropy.table import Table
from astropy.io import fits
import os
import glob

multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')
TIRDir = os.path.join(multiDir,'data','TIR','convolved15arc')
z0mgsDir = os.path.join(multiDir,'data','interpolated_z0mgs')

outDir = os.path.join(multiDir,'data','LTIR_calc')

if not os.path.exists(outDir):
    os.mkdir(outDir)

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

## should I just do dr1 galaxies here?
for galaxy in degas:
    print("***Processing galaxy: " + galaxy['NAME'] + "***")
    
    # look for 24micron images        
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_lvl_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_lvl_release_mips24_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_bendomips_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_bendomips_release_mips24_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_sings_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_sings_release_mips24_gauss15.fits")
    else:
        b24 = None
        print("24micron image not found")

    # look for 70micron images
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs70_gauss15.fits")):
        b70 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs70_gauss15.fits")
    elif  os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs70_gauss15.fits")):
        b70 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs70_gauss15.fits")
    else:
        b70 = None
        print("70micron image not found")


    # look for 100micron image
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs100_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs100_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs100_gauss15.fits")
    else:
        b100 = None
        print("100micron image not found")

    #look for 160 micron image
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs160_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs160_gauss15.fits")
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs160_gauss15.fits")
    else:
        b160 = None
        print("160micron image not found")
   
    # look for w4 interpolated images
    if os.path.exists(os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w4_gauss15_interpol.fits")):
        w4 = os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w4_gauss15_interpol.fits")
    else:
        w4 = None
        print("w4 image not found.")    
     
    # look for w3 interpolated images
    # if os.path.exists(os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w3_gauss15_interpol.fits")):
    #     w3 = os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w3_gauss15_interpol.fits")
    # else:
    #     w3 = None
    #     print("w3 image not found.")
 
    outFile = os.path.join(outDir,galaxy['NAME']+"_LTIR_gauss15.fits")
    calc_LTIR(outFile,b24=b24,b70=b70, b100=b100,b160=b160, w4=w4)
