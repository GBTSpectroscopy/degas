# Produce mom0/mom1 products
from degas.analysis_products import make_line_products
from astropy.table import Table
import os
import glob
import astropy.units as u

release = 'IR6p1'
regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
maskDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
degas_dr1 = degas[degas['DR1'] == 1]
#degas_dr1 = degas[degas['NAME'] == 'NGC2903']

line_list = ['HCN','HCOp','13CO','C18O']

noise_kwargs = {'do_map':True,'do_spec':True,'spec_box':5}

outDir = os.path.join(os.environ['ANALYSISDIR'],'moments_'+release)
if not os.path.exists(outDir):
    os.mkdir(outDir)

for galaxy in degas_dr1:
    for line in line_list:
        if (galaxy['NAME'] == 'NGC6946') & ((line == '13CO') | (line == 'C18O')):
            continue
        else:
            make_line_products(galaxy=galaxy['NAME'], 
                               line=line, 
                               inDir=regridDir, 
                               maskDir=maskDir, 
                               outDir=outDir, 
                               noise_kwargs=noise_kwargs,
                               line_width=10*u.km/u.s, 
                               snr_cut=[3.0,5.0],
                               vel_tolerance=50.0) 



