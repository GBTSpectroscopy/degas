from spectral_cube import SpectralCube
from degas.analysis_setup import fixALMAHCN

import glob
import os
import shutil
from astropy.table import Table

release = 'alma_hcn'

releaseDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','alma_hcn')

almalist = glob.glob(os.path.join(releaseDir, "*_7m+tp_hcn_pbcorr_trimmed_k.fits"))

for almacube in almalist:
    fixALMAHCN(almacube)
