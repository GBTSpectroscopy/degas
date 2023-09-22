import os
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from degas.analysis_plot import plot_stack, plot_trends

matplotlib.use('cairo')

#release = 'IR6p1'
#release = 'IR6p1_spatialR21'
#stacktype = 'mom1'
#stacktype = 'peakVelocity'

for release in ['IR6p1','IR6p1_spatialR21']:
    for stacktype in ['mom1','peakVelocity']:
        print("Creating plots for " + release + " " + stacktype)

        if stacktype == 'mom1':
            stackDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)
        elif stacktype == 'peakVelocity':
            stackDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_'+stacktype)
        else:
            print(stacktype+" unknown. Assuming directory name of stack_release.")
            stackDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)

        # setup information sources
        stack = Table.read(os.path.join(stackDir,'stack_'+release+'_'+stacktype+'.fits'))
        stack_pruned = Table.read(os.path.join(stackDir,'stack_'+release+'_'+stacktype+'_pruned.fits'))

        degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

        # Plot individual stacks

        plot_dir =  os.path.join(stackDir,'stack_plots')

        binlist = ['radius','r25','mstar','ICO', 'molgas','PDE']
        #binlist = ['radius']

        for bin_type in binlist:

            plot_stack(stack, bin_type,
                       plot_dir,
                       degas_db, release='DR1',fwhm_factor=5.0)

        # plot trends

        plot_dir =  os.path.join(stackDir,'stack_trends')

        binlist = ['radius','r25','mstar','ICO','molgas', 'PDE']

        mystyledict = {'int_intensity_sum_CO': 
                       {'marker':'o','color':'orange', 'label':'CO'}, 
                       'int_intensity_sum_HCN': 
                       {'marker':'^','color':'green', 'label':'HCN'}, 
                       'int_intensity_sum_HCOp': 
                       {'marker':'s','color':'blue', 'label':'HCO+'}, 
                       'int_intensity_sum_13CO': 
                       {'marker':'>','color':'red', 'label':'13CO'}, 
                       'int_intensity_sum_C18O': 
                       {'marker':'D','color': 'magenta', 'label':'C18O'},
                       'ratio_HCN_CO':
                       {'marker':'^','color':'green', 'label':r'HCN/$^{12}$CO'},
                        'ratio_HCOp_CO':
                       {'marker':'s','color':'blue', 'label':r'HCO+/$^{12}$CO'},
                       'ratio_ltir_mean_HCN':
                       {'marker':'^','color':'green', 'label':r'L$_{TIR}$/HCN (L$_\odot$ (K km s$^{-1}$ pc$^2$)$^{-1}$'},
                       'ratio_ltir_mean_HCOp':
                       {'marker':'^','color':'blue', 'label':r'L$_{TIR}$/HCO+ (L$_\odot$ (K km s$^{-1}$ pc$^2$)$^{-1}$'},
                       'ratio_HCOp_HCN':
                       {'marker':'s','color':'blue', 'label':r'HCO+/HCN'},
                       'ratio_13CO_CO':
                       {'marker':'>','color':'red', 'label':r'$^{13}$CO/$^{12}$CO'},
                       'ratio_C18O_CO':
                       {'marker':'D','color': 'magenta', 'label':r'$C^{18}$O/$^{12}$CO'},
                       'molgas_mean': 
                       {'marker':'o', 'color':'orange','label':r'$\Sigma_{mol}$ (M$_\odot/pc^2$)'},
                       'mstar_mean':
                       {'marker':'o','color':'red','label':r'$\Sigma_*$ (M$_\odot/pc^2$)'},
                       'sfr_mean':
                       {'marker':'o','color':'blue','label':r'$\Sigma_{SFR}$ (M$_\odot/yr/pc^2$)'}}


        for bin_type in binlist:

            if bin_type in ['mstar', 'ICO','molgas','PDE']:
                xlog = True
            else:
                xlog = False

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['int_intensity_sum_CO',
                         'int_intensity_sum_HCN',
                         'int_intensity_sum_HCOp',
                         'int_intensity_sum_13CO',
                         'int_intensity_sum_C18O'],
                        bin_type+'_lines',
                    #factordict= {'int_intensity_sum_CO':30},
                        styledict = mystyledict,
                        yaxislabel = r'Stacked Integrated Intensity (K km s$^{-1}$)',
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['molgas_mean','mstar_mean','sfr_mean'],
                        bin_type + '_other',
                        factordict = {'sfr_mean': 1/1e10},
                        yaxislabel = r'Variable',
                        styledict = mystyledict,
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ltir_mean','ltir_total'],
                        bin_type + '_ltir',
                        yaxislabel = r'Variable',
                        #styledict = mystyledict,
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ratio_ltir_mean_HCN','ratio_ltir_mean_HCOp'],
                        bin_type + '_sfedense',
                        yaxislabel = r'L$_{TIR}$/HCN',
                        styledict = mystyledict,
                        xlog=xlog)         

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ratio_HCN_CO','ratio_HCOp_CO'],
                        bin_type+'_fdense',
                        styledict = mystyledict,
                        yaxislabel = r'Dense gas fraction',
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ratio_HCN_CO','ratio_HCOp_CO','ratio_13CO_CO','ratio_C18O_CO'],
                        bin_type+'_coratios',
                        styledict = mystyledict,
                        yaxislabel = r'Dense gas ratios',
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ratio_HCOp_HCN'],
                        bin_type+'_hcop_hcn',
                        yaxislabel = r'Ratios',
                        styledict = mystyledict,
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack, bin_type, plot_dir, degas_db,
                        ['ratio_13CO_C18O'],
                        'radius_13CO_C18O',
                        yaxislabel = r'Ratios',
                        styledict = mystyledict,
                        ylog=False,
                        xlog=xlog)

        # plot individual stacks after pruning
        plot_dir =  os.path.join(stackDir,'stack_plots_pruned')


        binlist = ['radius','r25','mstar','ICO', 'molgas','PDE']

        for bin_type in binlist:

            plot_stack(stack_pruned, bin_type,
                       plot_dir,
                       degas_db, release='DR1', fwhm_factor=5.0)

        # plot trends
        plot_dir =  os.path.join(stackDir,'stack_trends_pruned')

        binlist = ['radius','r25','mstar','ICO','molgas','PDE']

        for bin_type in binlist:

            if bin_type in ['mstar', 'ICO','molgas','PDE']:
                xlog = True
            else:
                xlog = False

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['int_intensity_sum_CO',
                         'int_intensity_sum_HCN',
                         'int_intensity_sum_HCOp',
                         'int_intensity_sum_13CO',
                         'int_intensity_sum_C18O'],
                        bin_type+'_lines',
                    #factordict= {'int_intensity_sum_CO':30},
                        styledict = mystyledict,
                        yaxislabel = r'Stacked Integrated Intensity (K km s$^{-1}$)',
                        xlog=xlog)


            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['molgas_mean','mstar_mean','sfr_mean'],
                        bin_type + '_other',
                        factordict = {'sfr_mean': 1/1e10},
                        styledict = mystyledict,
                        yaxislabel = r'Variable',
                        xlog=xlog)

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ltir_mean','ltir_total'],
                        bin_type + '_ltir',
                        yaxislabel = r'Variable',
                        #styledict = mystyledict,
                        xlog=xlog)        

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ratio_ltir_mean_HCN','ratio_ltir_mean_HCOp'],
                        bin_type + '_sfedense',
                        yaxislabel = r'L$_{TIR}$/HCN',
                        styledict = mystyledict,
                        xlog=xlog)         

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ratio_HCN_CO','ratio_HCOp_CO'],
                        bin_type+'_fdense',
                        styledict = mystyledict,
                        yaxislabel = r'Dense gas fraction',
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ratio_HCN_CO','ratio_HCOp_CO','ratio_13CO_CO','ratio_C18O_CO'],
                        bin_type+'_coratios',
                        styledict = mystyledict,
                        yaxislabel = r'Dense gas ratios',
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ratio_HCOp_HCN'],
                        bin_type+'_hcop_hcn',
                        yaxislabel = r'Ratios',
                        styledict = mystyledict,
                        ylog=False,
                        xlog=xlog)

            plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                        ['ratio_13CO_C18O'],
                        'radius_13CO_C18O',
                        yaxislabel = r'Ratios',
                        styledict = mystyledict,
                        ylog=False,
                        xlog=xlog)
