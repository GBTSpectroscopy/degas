from degas.analysis_setup import smoothCube
from spectral_cube import SpectralCube

import glob
import os
import shutil
from astropy.table import Table

release = 'empire'
releaseDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','empire','EMPIRE_cubes')

multiDir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','multiwavelength')

z0mgsDir = os.path.join(multiDir,'data','interpolated_z0mgs')

TIRDir_15as = os.path.join(multiDir,'data','TIR','convolved15arc')
TIRDir_33as = os.path.join(multiDir,'data','TIR','convolved33arc')

hcnlist = glob.glob(os.path.join(releaseDir,'*_hcn_*_fixed_kms.fits'))
hcn = hcnlist[0]
hcn_cube = SpectralCube.read(hcn)
beam = hcn_cube.beam.major.to('arcsec').value

TIRlist = glob.glob(os.path.join(TIRDir_15as,"*.fits"))

if not os.path.exists(TIRDir_33as):
    os.mkdir(TIRDir_33as)

for TIR in TIRlist:
    smoothCube(TIR,TIRDir_33as, beam=beam)


w4list = glob.glob(os.path.join(z0mgsDir,'*_w4_gauss15_interpol.fits'))

for w4 in w4list:
    smoothCube(w4,z0mgsDir,beam=beam)

