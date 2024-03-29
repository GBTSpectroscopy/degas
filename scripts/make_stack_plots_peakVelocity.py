import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

from degas.analysis_plot import plot_stack, plot_trends

release = 'IR6p1'
stacktype = 'peakVelocity'

# setup information sources
stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_'+stacktype,b'stack_'+release+'_'+stacktype+'.fits'))
stack_pruned = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_'+stacktype,'stack_'+release+'_'+stacktype+'_pruned.fits'))
degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

# Plot individual stacks

plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'st,'stack_plots')

binlist = ['radius','r25','mstar','ICO']
#binlist = ['radius']

for bin_type in binlist:

    plot_stack(stack, bin_type,
               plot_dir,
               degas_db, release='DR1')
    
# plot trends

plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_'+release,'stack_trends')

binlist = ['radius','r25','mstar','ICO']

mystyledict = {'int_intensity_sum_CO': 
               {'marker':'o','color':'orange'}, 
               'int_intensity_sum_HCN': 
               {'marker':'^','color':'green'}, 
               'int_intensity_sum_HCOp': 
               {'marker':'s','color':'blue'}, 
               'int_intensity_sum_13CO': 
               {'marker':'>','color':'red'}, 
               'int_intensity_sum_C18O': 
               {'marker':'D','color': 'magenta'},
               'ratio_HCN_CO':
               {'marker':'^','color':'green'},
               'ratio_HCOp_CO':
               {'marker':'s','color':'blue'},
               'ratio_13CO_CO':
               {'marker':'>','color':'red'},
               'ratio_C18O_CO':
               {'marker':'D','color': 'magenta'}}

for bin_type in binlist:

    if bin_type in ['mstar', 'ICO']:
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
                ['comass_mean','mstar_mean','sfr_mean'],
                bin_type + '_other',
                factordict = {'sfr_mean': 1/1e10},
                yaxislabel = r'Variable',
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
                ylog=False,
                xlog=xlog)

    plot_trends(stack, bin_type, plot_dir, degas_db,
                ['ratio_13CO_C18O'],
                'radius_13CO_C18O',
                yaxislabel = r'Ratios',
                ylog=False,
                xlog=xlog)

# plot individual stacks after pruning

plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_'+release,'stack_plots_pruned')


binlist = ['radius','r25','mstar','ICO']

for bin_type in binlist:

    plot_stack(stack_pruned, bin_type,
               plot_dir,
               degas_db, release='DR1')

# plot trends

plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_'+release,'stack_trends_pruned')

binlist = ['radius','r25','mstar','ICO']

for bin_type in binlist:

    if bin_type in ['mstar', 'ICO']:
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
                ['comass_mean','mstar_mean','sfr_mean'],
                bin_type + '_other',
                factordict = {'sfr_mean': 1/1e10},
                yaxislabel = r'Variable',
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
                ylog=False,
                xlog=xlog)

    plot_trends(stack_pruned, bin_type, plot_dir, degas_db,
                ['ratio_13CO_C18O'],
                'radius_13CO_C18O',
                yaxislabel = r'Ratios',
                ylog=False,
                xlog=xlog)
