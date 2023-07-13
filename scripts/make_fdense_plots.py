from astropy.table import Table
from degas.analysis_plot import ratio_plots
import os

degas = Table.read('/lustre/cv/users/akepley/degas/code/degas/scripts/degas_base.fits')

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends'

# IR6p1 mom1
stack = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_IR6p1_mom1_pruned.fits')

ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_HCN_CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')

ratio_plots(degas, stack, ydata='ratio_HCN_13CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_13CO_vs_r25_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCN_13CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_13CO_vs_mstar_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCN_13CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_13CO_vs_molgas_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)

ratio_plots(degas, stack, ydata='ratio_HCN_C18O',xdata='r25',pltname=os.path.join(plotdir,"fdense_C18O_vs_r25_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCN_C18O',xdata='mstar',pltname=os.path.join(plotdir,"fdense_C18O_vs_mstar_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCN_C18O',xdata='molgas',pltname=os.path.join(plotdir,"fdense_C18O_vs_molgas_IR6p1_mom1"),pltlabel=r'constant R$_{21}$',add_empire=False)

ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='mstar',pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1"),pltlabel=r'constant R$_{21}$')


ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='r25', pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_CO',xdata='molgas', pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)

ratio_plots(degas, stack, ydata='ratio_HCOp_13CO',xdata='r25', pltname=os.path.join(plotdir,"fdense_13CO_vs_r25_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_13CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_13CO_vs_mstar_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_13CO',xdata='molgas', pltname=os.path.join(plotdir,"fdense_13CO_vs_molgas_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)

ratio_plots(degas, stack, ydata='ratio_HCOp_C18O',xdata='r25', pltname=os.path.join(plotdir,"fdense_C18O_vs_r25_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_C18O',xdata='mstar',pltname=os.path.join(plotdir,"fdense_C18O_vs_mstar_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)
ratio_plots(degas, stack, ydata='ratio_HCOp_C18O',xdata='molgas', pltname=os.path.join(plotdir,"fdense_C18O_vs_molgas_IR6p1_mom1_hcop"),pltlabel=r'constant R$_{21}$',add_empire=False)


ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_hcop"),add_empire=False,pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='mstar', pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_hcop"),add_empire=False,pltlabel=r'constant R$_{21}$')
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCOp',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1_hcop"),add_empire=False,pltlabel=r'constant R$_{21}$')



# IR6p1 mom1 spatial R21
stack2 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_IR6p1_spatialR21_mom1_pruned.fits')
stack2_fit =  Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits.fits')
#stack2_fit_noNGC4414 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits_noNGC4414.fits')

ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='r25', pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_spatialR21"), stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_spatialR21"), stack_fit=stack2_fit,xlim=(1.5,4.0),ylim=(-3,0))
ratio_plots(degas, stack2, ydata='ratio_HCN_CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1_spatialR21"), stack_fit=stack2_fit)

ratio_plots(degas, stack2, ydata='ratio_HCN_13CO',xdata='r25', pltname=os.path.join(plotdir,"fdense_13CO_vs_r25_IR6p1_mom1_spatialR21") ,add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCN_13CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_13CO_vs_mstar_IR6p1_mom1_spatialR21"), xlim=(1.5,4.0),ylim=(-1.5,-1.5+3.0), add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCN_13CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_13CO_vs_molgas_IR6p1_mom1_spatialR21"),  add_empire=False)

ratio_plots(degas, stack2, ydata='ratio_HCN_C18O',xdata='r25', pltname=os.path.join(plotdir,"fdense_C18O_vs_r25_IR6p1_mom1_spatialR21"), add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCN_C18O',xdata='mstar',pltname=os.path.join(plotdir,"fdense_C18O_vs_mstar_IR6p1_mom1_spatialR21"), add_empire=False)
ratio_plots(degas, stack2, ydata='ratio_HCN_C18O',xdata='molgas',pltname=os.path.join(plotdir,"fdense_C18O_vs_molgas_IR6p1_mom1_spatialR21"), add_empire=False)

ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCN',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_spatialR21"),add_empire=True, stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCN',xdata='mstar',pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_spatialR21"),add_empire=True,stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCN',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1_spatialR21"),add_empire=True,stack_fit=stack2_fit)


ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_HCOp_CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)

ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCOp',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCOp',xdata='mstar',pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)
ratio_plots(degas, stack2, ydata='ratio_ltir_mean_HCOp',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1_hcop_spatialR21"),add_empire=False, stack_fit=stack2_fit)


# IR6p1 mom1 33as
stack3 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_33as/stack_IR6p1_33as_mom1_pruned.fits')
#stack3_fit = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_33as/stack_fits.fits')
#stack3_fit_noNGC4414 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_33as/stack_fits_noNGC4414.fits')

ratio_plots(degas, stack3, ydata='ratio_HCN_CO',xdata='r25',pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_33as"))
ratio_plots(degas, stack3, ydata='ratio_HCN_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_33as"))
ratio_plots(degas, stack3, ydata='ratio_HCN_CO',xdata='molgas',pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1_33as"))

ratio_plots(degas, stack3, ydata='ratio_ltir_mean_HCN',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_33as"))
ratio_plots(degas, stack3, ydata='ratio_ltir_mean_HCN',xdata='mstar', pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_33as"))
ratio_plots(degas, stack, ydata='ratio_ltir_mean_HCN',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1_33as"))

ratio_plots(degas, stack3, ydata='ratio_HCOp_CO',xdata='r25', pltname=os.path.join(plotdir,"fdense_vs_r25_IR6p1_mom1_hcop_33as"), add_empire=False)
ratio_plots(degas, stack3, ydata='ratio_HCOp_CO',xdata='mstar',pltname=os.path.join(plotdir,"fdense_vs_mstar_IR6p1_mom1_hcop_33as"), add_empire=False)
ratio_plots(degas, stack3, ydata='ratio_HCOp_CO',xdata='molgas', pltname=os.path.join(plotdir,"fdense_vs_molgas_IR6p1_mom1_hcop_33as"), add_empire=False)

ratio_plots(degas, stack3, ydata='ratio_ltir_mean_HCOp',xdata='r25',pltname=os.path.join(plotdir,"sfe_dense_vs_r25_IR6p1_mom1_hcop_33as"),add_empire=False)
ratio_plots(degas, stack3, ydata='ratio_ltir_mean_HCOp',xdata='mstar', pltname=os.path.join(plotdir,"sfe_dense_vs_mstar_IR6p1_mom1_hcop_33as"),add_empire=False)
ratio_plots(degas, stack3, ydata='ratio_ltir_mean_HCOp',xdata='molgas',pltname=os.path.join(plotdir,"sfe_dense_vs_molgas_IR6p1_mom1_hcop_33as"),add_empire=False)


