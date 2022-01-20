import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import ipdb

## TODO?

## Could I just define some default plot settings here and reference below?


def plot_stack(stack, bin_type, plot_dir, degas_db, release='DR1'):
    '''
    general purpose routine to plot stacks.

    stack: astropy Table with stack

    bin_type: type of bin to plot

    plot_dir: plotting directory

    degas_db: degas galaxy property database in astropy Table format

    Date        Programmer      Description of Code
    ----------------------------------------------------------------------
    12/21/2021  A.A. Kepley     Original Code

    '''


    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    binlist = np.unique(stack['bin_type'])

    if bin_type not in binlist:
        print(bin_type + " not found in stack.")
        return

    # plot setup
    style = {'CO': {'marker':'o','color':'orange','name':'CO'}, 
             'HCN': {'marker':'^','color':'green','name':'HCN'}, 
             'HCOp': {'marker':'s','color':'blue','name':'HCO+'}, 
             '13CO': {'marker':'*','color':'red', 'name':'13CO'}, 
             'C18O': {'marker':'D','color': 'magenta','name':'C18O'}}
    
    # only look at dr1 galaxies
    gallist = degas_db[release] == 1

    for galaxy in degas_db[gallist]:

        plt.close()

        idx = ( (stack['galaxy'] == galaxy['NAME']) \
                & (stack['bin_type'] == bin_type))
    
        nprofile = np.sum(idx)
        
        ncols = 3
        nrows = int(np.ceil(nprofile / ncols))
    
        fig, myax = plt.subplots(nrows=nrows, ncols=ncols,
                                 sharex='all', sharey='row',
                                 figsize=(8,8))

        fig.subplots_adjust(hspace=0.1,wspace=0.1)

        flatax = myax.flatten()

        for i in range(0,nprofile):
        
            if galaxy['NAME'] == 'NGC6946':
                linelist = ['CO', 'HCN','HCOp']
            else:
                linelist = ['CO', 'HCN','HCOp','13CO','C18O']

            for line in linelist:

                if line == 'CO':
                    factor = 30
                    mylabel = style[line]['name'] + "/"+str(factor)
                else:
                    factor = 1.0
                    mylabel =  style[line]['name']
        
                flatax[i].plot(stack[idx]['spectral_axis'][i],
                               stack[idx]['stack_profile_'+line][i] / factor,
                               color=style[line]['color'],
                               linestyle='-',
                               marker=None,
                               label=mylabel,
                               linewidth=1)
            
                flatax[i].axhline(0,linestyle=':',linewidth=1,color='gray')
                bin_label = 'bin = {:5.2f}' .format(stack[idx]['bin_mean'][i])+ ' ' + str(stack[idx]['bin_unit'][0])
                flatax[i].text(0.1, 0.9, bin_label,
                               transform=flatax[i].transAxes,fontsize='x-small')
                flatax[i].legend(loc='upper right',fontsize='x-small')

            
        for ax in flatax:
            ax.set(xlabel=r'Velocity (km s$^{-1}$)',
                   ylabel=r'Stacked Mean Intensity (K)')
    
        for ax in flatax:
            ax.label_outer()

        plt.show()
        plt.savefig(os.path.join(plot_dir,galaxy['NAME']+'_'+bin_type+'_stacks.pdf'))
        plt.savefig(os.path.join(plot_dir,galaxy['NAME']+'_'+bin_type+'_stacks.png'))


def plot_trends(stack, bin_type, plot_dir, degas_db, 
                plot_quant,
                plot_name, 
                styledict = None,
                factordict = None,
                release='DR1',
                ylog = True,
                xlog = False,
                yaxislabel = 'Quantity'):
    '''
    Plot radial trends for given stacks

    stack: astropy Table with stack

    bin_type: type of bin to plot

    plot_dir: plotting directory. will be create if doesn't exist

    degas_db: degas galaxy property database in astropy Table format

    plot_quant: quantities to plot. Much be columns in stack

    plot_name: name for plot file. format will be <galaxy name>_<plot_name>.pdf or *.png

    styledict: dictionary to style plots

    factordict: dictionary of factors to divide data by
    
    release: select data from a particular release
    
    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    1/6/2022    A.A. Kepley     Original Code

    '''

    # create plot directory if doesn't already exist
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # check to make sure bin_type exists in stack.
    binlist = np.unique(stack['bin_type'])

    if bin_type not in binlist:
        print(bin_type + " not found in stack.")
        return

    #only look at dr1 galaxies
    gallist = degas_db[release] == 1

    # check to make sure plot quantities exist in stacks.
    for quant in plot_quant:
        if quant not in stack.columns:
            print("Unknown quantity " + quant + "\n")
            return

    # setup default plot colors and markers
    nquant = len(plot_quant)
    if not styledict:
        markers = ['o','v','^','s','*','D'] # 6 items
        markerlist = np.tile(markers,int(np.ceil(nquant/len(markers))))
        markerlist = markerlist[0:nquant]

        colors = ['royalblue','forestgreen','darkorange','royalblue','crimson','rebeccapurple','darkcyan','darkmagenta']
        colorlist = np.tile(colors,int(np.ceil(nquant/len(colors))))
        colorlist = colorlist[0:nquant]

        styledict = {}

        for (quant,marker,color) in zip(plot_quant, markerlist,colorlist):
            styledict[quant] = {'marker': marker, 'color': color}

    # setup label information 
    ## TO DO: move this up top so read in automatically?
    labeldict = {'int_intensity_sum_CO': 'CO',
                 'int_intensity_sum_HCN': 'HCN',
                 'int_intensity_sum_HCOp': 'HCO+',
                 'int_intensity_sum_13CO': '13CO',
                 'int_intensity_sum_C18O': 'C18O',
                 'comass_mean': r'$\Sigma_{CO}$ (M$_\odot/pc^2$)',
                 'mstar_mean': r'$\Sigma_*$ (M$_\odot/pc^2$)',
                 'sfr_mean': r'$\Sigma_{SFR}$ (M$_\odot/yr/pc^2$)',
                 'ratio_HCN_CO': r'HCN/$^{12}$CO',
                 'ratio_HCOp_CO': r'HCO+/$^{12}$CO',
                 'ratio_13CO_CO': r'$^{13}$CO/$^{12}$CO',
                 'ratio_C18O_CO': r'$C^{18}$O/$^{12}$CO',
                 'ratio_HCOp_HCN': r'HCO+/HCN',
                 'ratio_13CO_C18O': r'$^{13}$CO/$C^{18}$O'}

    for galaxy in degas_db[gallist]:

        plt.close()
        
        # select data
        idx = ( (stack['galaxy'] == galaxy['NAME'] ) & (stack['bin_type'] == bin_type))

        # convert x-axis if needed
        if bin_type == 'radius':
            # convert to kpc -- assuming all bins have same unit
            xvals = (stack[idx]['bin_mean'] * u.Unit(stack[idx]['bin_unit'][0])) * galaxy['DIST_MPC'] * u.Mpc
            xvals = xvals.to('kpc',equivalencies=u.dimensionless_angles())
            xvals = xvals.data
            xlabel = r"Galactocentric Radius (kpc)"
        elif bin_type == 'r25':
            xvals = stack[idx]['bin_mean']
            xlabel = r'r/r$_{25}$'
        elif bin_type == 'mstar':
            xvals = stack[idx]['bin_mean']
            xlabel = r'$\Sigma_{*}$ (' + u.Unit(stack[idx]['bin_unit'][0]).to_string('latex_inline') + ')'
        elif bin_type == 'ICO':
            xvals = stack[idx]['bin_mean']
            xlabel = r'$I_{CO}$ (' + u.Unit(stack[idx]['bin_unit'][0]).to_string('latex_inline') + ')'
        else:
            # keep the units as is for the data
            xvals = stack[idx]['bin_mean'] 
            xlabel = stack[idx]['bin_type'][0] + ' (' +  stack[idx]['bin_unit'][0] + ')'

        # plot each quantity
        for quant in plot_quant:

            # get factor
            if factordict:            
                if quant in factordict:
                    factor = factordict[quant]
                    
                    if factor >= 1:
                        mylabel = labeldict[quant] + '/' + str(factor)
                    else:
                        mylabel = labeldict[quant] + "*" + "{:5.1e}".format(1/factor)
                else:
                    factor = 1.0
                    mylabel = labeldict[quant]                
            else:
                factor = 1.0
                mylabel = labeldict[quant]
                    

            # get values to plot
            yvals = stack[idx][quant].data / factor
            
            if quant+ '_err' in stack.columns:
                # get errors
                yerrs = stack[idx][quant+'_err'].data / factor
            else:
                yerrs = None

            # get upper and lower limits
            if quant + '_uplim' in stack.columns:
                uplims = stack[idx][quant + '_uplim']
                yerrs[uplims] = yvals[uplims] * 0.3
            else:
                uplims = False

            if quant + '_lolim' in stack.columns:
                lolims = stack[idx][quant + '_lolim']
                yerrs[lolims] = yvals[lolims] * 0.3
            else:
                lolims = False

            # make plot 

            # TODO -- I want empty symbols for limits and filled symbols
            # for measurements.
            plt.errorbar(xvals, yvals,
                         yerr = yerrs,
                         uplims = uplims,
                         lolims = lolims,
                         marker = styledict[quant]['marker'], 
                         color = styledict[quant]['color'], 
                         linestyle = '--',
                         label = mylabel) 

            
        # add info to plot
        if ylog:
            plt.yscale('log')
        
        if xlog:
            plt.xscale('log')

        plt.legend() ### want to get rid of error bars here....
        plt.title(galaxy['NAME'])
        plt.xlabel(xlabel)
        plt.ylabel(yaxislabel)
        plt.show()
            
        # save plot
        plt.savefig(os.path.join(plot_dir,galaxy['NAME'] + '_' + plot_name + '.png'))
        plt.savefig(os.path.join(plot_dir,galaxy['NAME'] + '_' + plot_name + '.pdf'))

        plt.close()






