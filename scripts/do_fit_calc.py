from astropy.table import Table
import os
from degas.analysis_fit import fit_trend,make_fit_table

release = 'IR6p1'

stack_dir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)

stack = Table.read(os.path.join(stack_dir,'stack_IR6p1_mom1_pruned.fits'))

myresult = make_fit_table(stack,outfile=os.path.join(stack_dir,'stack_fits.fits'),plotDir=os.path.join(stack_dir,'fit_plots'))

myresult_noNGC4414 = make_fit_table(stack,outfile=os.path.join(stack_dir,'stack_fits_noNGC4414.fits'),plotDir=os.path.join(stack_dir,'fit_plots_noNGC4414'),exclude_gal=['NGC4414'])

#test2 = make_fit_table(stack,outfile=os.path.join(stack_dir,'stack_fits_noNGC4414.fits'),exclude_gal=['NGC4414'])

