from astropy.table import Table
import os
from degas.analysis_fit import fit_trend,make_fit_table

release = 'IR6p1'

# R21 simple

stack_dir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)

stack = Table.read(os.path.join(stack_dir,'stack_IR6p1_mom1_pruned.fits'))

myresult = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'],outfile=os.path.join(stack_dir,'stack_fits.fits'),plotDir=os.path.join(stack_dir,'fit_plots'))

myresult_noNGC4414 = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'], outfile=os.path.join(stack_dir,'stack_fits_noNGC4414.fits'),plotDir=os.path.join(stack_dir,'fit_plots_noNGC4414'),exclude_gal=['NGC4414'])

## R21 spatial

stack_dir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_spatialR21')

stack = Table.read(os.path.join(stack_dir,'stack_IR6p1_spatialR21_mom1_pruned.fits'))

myresult = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'],outfile=os.path.join(stack_dir,'stack_fits.fits'),plotDir=os.path.join(stack_dir,'fit_plots'))

myresult_noNGC4414 = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'], outfile=os.path.join(stack_dir,'stack_fits_noNGC4414.fits'),plotDir=os.path.join(stack_dir,'fit_plots_noNGC4414'),exclude_gal=['NGC4414'])

# 33 arcsec

stack_dir = '/lustre/cv/users/akepley/degas/stack_IR6p1_33as/'

stack = Table.read(os.path.join(stack_dir,'stack_IR6p1_33as_mom1_pruned.fits'))

myresult = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'],outfile=os.path.join(stack_dir,'stack_fits.fits'),plotDir=os.path.join(stack_dir,'fit_plots'))

myresult_noNGC4414 = make_fit_table(stack,columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN','ratio_HCOp_CO','ratio_ltir_mean_HCOp'],outfile=os.path.join(stack_dir,'stack_fits_noNGC4414.fits'),plotDir=os.path.join(stack_dir,'fit_plots'),exclude_gal=['NGC4414'])
