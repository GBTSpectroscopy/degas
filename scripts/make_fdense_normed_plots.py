import os
from astropy.table import Table
from degas.analysis_plot import ratio_plots

release =  'IR6p1'
scriptDir = os.environ['SCRIPTDIR']
stackDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_spatialR21')

degas_db = Table.read(os.path.join(scriptDir,'degas_base.fits'))
stack_norm = Table.read(os.path.join(stackDir,'stack_IR6p1_spatialR21_mom1_pruned_norm.fits'))

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends_normalized'

stack_norm_fit = Table.read(os.path.join(stackDir,'stack_norm_fits.fits'))

### HCN/CO

ratio_plots(degas_db, stack_norm, ydata='ratio_HCN_CO_norm',
            xdata = 'r25',
            add_empire=False,
            norm=True,
            stack_fit = stack_norm_fit,
            pltname=os.path.join(plotdir, 'fdense_vs_r25_IR6p1_mom1_norm'))

ratio_plots(degas_db, stack_norm, ydata='ratio_HCN_CO_norm',
            xdata = 'mstar',
            add_empire=False,
            norm=True,
            stack_fit = stack_norm_fit,
            pltname=os.path.join(plotdir, 'fdense_vs_mstar_IR6p1_mom1_norm'))

ratio_plots(degas_db, stack_norm, ydata='ratio_HCN_CO_norm',
            xdata = 'molgas',
            add_empire=False,
            norm=True,
            stack_fit = stack_norm_fit,
            pltname=os.path.join(plotdir, 'fdense_vs_molgas_IR6p1_mom1_norm'))

### LTIR/HCN

ratio_plots(degas_db, stack_norm, ydata='ratio_ltir_mean_HCN_norm',
            xdata = 'r25',
            add_empire=False,
            norm=True,
            pltname=os.path.join(plotdir, 'sfe_dense_vs_r25_IR6p1_mom1_norm'))

ratio_plots(degas_db, stack_norm, ydata='ratio_ltir_mean_HCN_norm',
            xdata = 'mstar',
            add_empire=False,
            norm=True,
            stack_fit = stack_norm_fit,
            pltname=os.path.join(plotdir, 'sfe_dense_vs_mstar_IR6p1_mom1_norm'))

ratio_plots(degas_db, stack_norm, ydata='ratio_ltir_mean_HCN_norm',
            xdata = 'molgas',
            add_empire=False,
            norm=True,
            stack_fit = stack_norm_fit,
            pltname=os.path.join(plotdir, 'sfe_dense_vs_molgas_IR6p1_mom1_norm'))

