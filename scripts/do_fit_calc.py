from astropy.table import Table
import os
from degas.analysis_fit import fit_trend

release = 'IR6p1'

stack_dir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)

stack = Table.read(os.path.join(stack_dir,'stack_IR6p1_mom1_pruned.fits'))

test = fit_trend(stack,bin_type='r25',stack_col='ratio_HCN_CO',exclude_gal='NGC4414')

test2 = fit_trend(stack,bin_type='mstar',stack_col='ratio_HCN_CO')
