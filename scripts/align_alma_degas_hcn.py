from degas.products import makeMap
from degas.analysis_setup import regridData, smoothCube
from degas.masking import cubemask
from spectral_cube import SpectralCube
from degas.analysis_products import calc_products
import astropy.units as u

import glob
import os
import shutil
from astropy.table import Table

alma_dir = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','alma_hcn')
degas_dir = os.path.join(os.environ['ANALYSISDIR'],'IR6p1')
out_dir = os.path.join(os.environ['ANALYSISDIR'],'alma_degas_regrid')

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

gal_list = ['NGC2903','NGC4321']

noise_kwargs = {'do_map':True,'do_spec':True,'spec_box':5}

for gal in gal_list:
    alma_hcn = os.path.join(alma_dir,gal.lower()+'_7m+tp_hcn_pbcorr_trimmed_k_kms.fits')
    degas_hcn = os.path.join(degas_dir,gal+'_HCN_rebase7_smooth1.3_hanning1_maxnchan.fits')

    # get ALMA hcn beam
    alma_cube = SpectralCube.read(alma_hcn)
    beam = alma_cube.beam.major.to('arcsec').value

    # smooth and regrid DEGAS to same coordinate system
    degas_smooth = smoothCube(degas_hcn,out_dir,beam=beam)
    degas_regrid = regridData(alma_hcn,degas_smooth,out_dir)

    
    # copy alma hcn file to same directory
    shutil.copy(alma_hcn,out_dir)

    ## create mask and make moment maps
    if gal in ['NGC4321','NGC2903']:
        cubemask(alma_hcn,
                 gal+'_HCN_mask.fits',
                 outDir = out_dir,
                 peakCut = 5.0,
                 lowCut = 3.0)

        # make mask
        maskFile  = os.path.join(out_dir,gal+'_HCN_mask.fits')        

        rms_file = degas_regrid.replace('.fits','_rms.fits')
        calc_products(degas_regrid,
                      maskFile,
                      rms_file,
                      out_dir = out_dir,
                      mom1_template_file=None,
                      noise_kwargs = noise_kwargs,
                      line_width=10*u.km/u.s, 
                      snr_cut=[4.0,4.0],
                      vel_tolerance=30.0)
                      
        
        rms_file = alma_hcn.replace('.fits','_rms.fits')
        calc_products(alma_hcn,
                      maskFile,
                      rms_file,
                      out_dir = out_dir,
                      mom1_template_file=None,
                      noise_kwargs = noise_kwargs,
                      line_width=10*u.km/u.s, 
                      snr_cut=[4.0,4.0],
                      vel_tolerance=30.0)
                      
        
