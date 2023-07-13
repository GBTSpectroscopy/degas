from astropy.table import Table
from degas.analysis_plot import ratio_vs_ratio_plots
import os

degas = Table.read('/lustre/cv/users/akepley/degas/code/degas/scripts/degas_base.fits')

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends'

stack = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_IR6p1_spatialR21_mom1_pruned.fits')
stack_fit =  Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits.fits')

ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_13CO_CO', ydata='ratio_HCN_CO',pltname=os.path.join(plotdir,'ratio_13CO_CO_ratio_HCN_CO_mstar'))
                     
ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_13CO_CO', ydata='ratio_HCN_13CO', pltname=os.path.join(plotdir,'ratio_13CO_CO_ratio_HCN_13CO_mstar'))
                     
ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_13CO_CO', ydata='ratio_HCOp_HCN',pltname=os.path.join(plotdir,'ratio_13CO_CO_ratio_HCOp_HCN_mstar'))
   
ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_13CO_CO', ydata='ratio_13CO_C18O',pltname=os.path.join(plotdir,'ratio_13CO_CO_ratio_13CO_C18O_mstar'))
   
ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_HCN_CO', ydata='ratio_13CO_C18O',pltname=os.path.join(plotdir,'ratio_HCN_CO_ratio_13CO_C18O_mstar'))
   

ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='ratio_HCN_CO', ydata='ratio_HCN_13CO', pltname=os.path.join(plotdir,'ratio_HCN_CO_ratio_HCN_13CO_mstar'))

ratio_vs_ratio_plots(degas, stack, bin_type='mstar', xdata='int_intensity_sum_CO', ydata='ratio_13CO_CO',pltname=os.path.join(plotdir,'ratio_int_intensity_sum_CO_ratio_13CO_CO'))
