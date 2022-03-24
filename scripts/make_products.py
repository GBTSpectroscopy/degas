#this is the driver script for making moment maps and associated error maps##
#the working functions are in /degas/degas/analysis_products.py
import sys
sys.path.append('/lustre/cv/users/ysong/degas/degas/degas')
from analysis_products import make_line_products
import os
import glob
#test
line_name='HCOp'
#in_dir='/lustre/cv/users/ysong/degas/IR6p1_regrid/' #smoothed
in_dir='/lustre/cv/users/ysong/degas/IR6p1_regrid_nosmooth/' #native
mask_dir='/lustre/cv/users/ysong/degas/IR6p1_regrid/'
#out_dir='/lustre/cv/users/ysong/degas/test_product/smoothed/'
out_dir='/lustre/cv/users/ysong/degas/test_product/native/' #native
noise_kwargs = {'do_map':True,'do_spec':True,'spec_box':5}

filelist=glob.glob(in_dir+'*'+line_name+'*.fits')
for file in filelist:
    name=os.path.basename(file).split('_')[0]
    make_line_products(galaxy=name, line=line_name, inDir=in_dir, maskDir=mask_dir, outDir=out_dir, noise_kwargs=noise_kwargs)


