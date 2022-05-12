from degas.analysis_setup import calc_LTIR, colorCorrect_24micron
from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import os
import glob

multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')
TIRDir = os.path.join(multiDir,'data','TIR','convolved15arc')
z0mgsDir = os.path.join(multiDir,'data','interpolated_z0mgs')

outDir = os.path.join(multiDir,'data','LTIR_calc')

if not os.path.exists(outDir):
    os.mkdir(outDir)

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

# add in column for data source
IRBands = ['IR_24micron','IR_70micron','IR_100micron','IR_160micron']

for band in IRBands:
    if band in degas.colnames:
        degas.remove_column(band)
    degas.add_column(Column(np.full_like(degas['NAME'],''),dtype='S15'),name=band)

## should I just do dr1 galaxies here?
for galaxy in degas:
    print("***Processing galaxy: " + galaxy['NAME'] + "***")
    
    # look for 24micron images        
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_bendomips_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_bendomips_release_mips24_gauss15.fits")
        b24_src = 'BENDO'

    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_lvl_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_lvl_release_mips24_gauss15.fits")
        b24_src = 'LVL'

    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_sings_release_mips24_gauss15.fits")):
        b24 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_sings_release_mips24_gauss15.fits")
        b24_src = 'SINGS'

    # look for w4 interpolated images
    elif os.path.exists(os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w4_gauss15_interpol.fits")):
        b24 = os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w4_gauss15_interpol.fits")
        b24_src = 'W4'

    else:
        b24 = None
        b24_src = ''
        print("24micron image not found")

    degas['IR_24micron'][degas['NAME'] == galaxy['NAME']] = b24_src
    

    # look for 70micron images
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs70_gauss15.fits")):
        b70 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs70_gauss15.fits")
        b70_src = 'KINGFISH'

    elif  os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs70_gauss15.fits")):
        b70 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs70_gauss15.fits")
        b70_src = 'VNGS'

    else:
        b70 = None
        b70_src = ''
        print("70micron image not found")

    degas['IR_70micron'][degas['NAME'] == galaxy['NAME']] = b70_src


    # look for 100micron image
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs100_gauss15.fits")
        b100_src = 'KINGFISH'

    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs100_gauss15.fits")
        b100_src = 'HRS'

    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs100_gauss15.fits")):
        b100 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs100_gauss15.fits")
        b100_src = 'VNGS'
    else:
        b100 = None
        b100_src = ''
        print("100micron image not found")

    degas['IR_100micron'][degas['NAME'] == galaxy['NAME']] = b100_src

    #look for 160 micron image
    if os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_kingfish_pacs160_gauss15.fits")
        b160_src = 'KINGFISH'
        
    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_hrs_pacs160_gauss15.fits")
        b160_src = 'HRS'

    elif os.path.exists(os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs160_gauss15.fits")):
        b160 = os.path.join(TIRDir,galaxy['NAME'].lower()+"_vngs_pacs160_gauss15.fits")
        b160_src = 'VNGS'

    else:
        b160 = None
        b160_src = ''
        print("160micron image not found")
   
    degas['IR_160micron'][degas['NAME'] == galaxy['NAME']] = b160_src

  
    # look for w3 interpolated images
    # if os.path.exists(os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w3_gauss15_interpol.fits")):
    #     w3 = os.path.join(z0mgsDir,galaxy['NAME'].lower()+"_w3_gauss15_interpol.fits")
    # else:
    #     w3 = None
    #     print("w3 image not found.")

    ## color correct the images that have only 24micron, but do have
    ## lower resolution 70micron.
    outFile = os.path.join(outDir,galaxy['NAME']+"_LTIR_gauss15.fits")
    if b24 and not b70 and not b100 and not b160 and b24_src is 'LVL':

        print('LVL only -- color correcting LTIR')

        TIRDir_30arc = TIRDir.replace('convolved15arc','over15arc')

        b70 = os.path.join(TIRDir_30arc,galaxy['NAME'].lower()+'_lvl_release_mips70_gauss30.fits')
        if os.path.exists(b70):
            colorCorrect_24micron(b24, b70,outFile)

    # otherwise just calculate the LTIR at the given resolution
    else:
        calc_LTIR(outFile,b24=b24,b70=b70, b100=b100,b160=b160)

degas.write(os.path.join(os.environ['SCRIPTDIR'],"degas_base.fits"),overwrite=True)
