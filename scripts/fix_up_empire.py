from degas.analysis_setup import fixEMPIRE
from spectral_cube import SpectralCube

import glob
import os
import shutil
from astropy.table import Table

release = 'empire'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','empire','EMPIRE_cubes')

empirelist = glob.glob(os.path.join(releaseDir,"EMPIRE_*_*_??as.fits"))

for cube in empirelist:
    fixEMPIRE(cube)
