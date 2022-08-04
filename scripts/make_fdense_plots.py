from astropy.table import Table
from degas.analysis_plot import ratio_plots
import os

degas = Table.read('/lustre/cv/users/akepley/degas/code/degas/scripts/degas_base.fits')

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends'

# IR6p1 mom1
stack = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_IR6p1_mom1_pruned.fits')

ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='ICO',pltname=os.path.join(plotdir,"fdense_vs_ICO_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')

ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1"))
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='mstar',pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1"))
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='ICO',pltname=os.path.join(plotdir,"sfe_dense_vs_ICO_IR6p1_mom1"))

ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='ICO',pltname=os.path.join(plotdir,"fdense_vs_ICO_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)

ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_hcop"),add_empire=False)
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='mstar',pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_hcop"),add_empire=False)
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='ICO',pltname=os.path.join(plotdir,"sfe_dense_vs_ICO_IR6p1_mom1_hcop"),add_empire=False)



# IR6p1 mom1 spatial R21
stack2 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_IR6p1_spatialR21_mom1_pruned.fits')

ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_spatialR21"),pltlabel=r'spatial R$_{21}$')
ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_spatialR21"),pltlabel=r'spatial R$_{21}$')
ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='ICO',pltname=os.path.join(plotdir,"fdense_vs_ICO_IR6p1_mom1_spatialR21"),pltlabel=r'spatial R$_{21}$')

ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_hcop_spatialR21"),pltlabel=r'spatial R$_{21}$',add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_hcop_spatialR21"),pltlabel=r'spatial R$_{21}$',add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='ICO',pltname=os.path.join(plotdir,"fdense_vs_ICO_IR6p1_mom1_hcop_spatialR21"),pltlabel=r'spatial R$_{21}$',add_empire=False)
