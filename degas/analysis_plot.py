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

    styledict: dictionary to style plots (includes labels)

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
        markers = ['o','v','^','s','>','D'] # 6 items
        markerlist = np.tile(markers,int(np.ceil(nquant/len(markers))))
        markerlist = markerlist[0:nquant]

        colors = ['royalblue','forestgreen','darkorange','royalblue','crimson','rebeccapurple','darkcyan','darkmagenta']
        colorlist = np.tile(colors,int(np.ceil(nquant/len(colors))))
        colorlist = colorlist[0:nquant]

        styledict = {}
        for (quant,marker,color) in zip(plot_quant, markerlist,colorlist):
            styledict[quant] = {'marker': marker, 'color': color, 'label':quant}


  
    for galaxy in degas_db[gallist]:

        plt.close()
        
        # select data
        idx = ( (stack['galaxy'] == galaxy['NAME'] ) & (stack['bin_type'] == bin_type))

        # convert x-axis if needed
        if bin_type == 'radius':
            # convert to kpc -- assuming all bins have same unit
            xvals = (stack[idx]['bin_mean'] * u.Unit(stack[idx]['bin_unit'][0])) * galaxy['DIST_MPC'] * u.Mpc
            xvals = xvals.to('kpc',equivalencies=u.dimensionless_angles())
            xvals = xvals.value
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
                
            # set to defaults if values don't exist in dictionary
            if quant in styledict.keys():
                if 'label' in styledict[quant].keys():
                    mylabel = styledict[quant]['label']
                else:
                    mylabel = quant

                if 'marker' in styledict[quant].keys(): 
                    mymarker = styledict[quant]['marker']
                else:
                    mymarker = 'o'
                    
                if 'color' in styledict[quant].keys():
                    mycolor = styledict[quant]['color']
                else:
                    mycolor = 'green'

            else:
                mylabel = quant
                mymarker = 'o'
                mycolor = 'green'
                

            # get factor
            if factordict:            
                if quant in factordict:
                    factor = factordict[quant]
                    
                    if factor >= 1:
                        mylabel = mylabel + '/' + str(factor)
                    else:
                        mylabel = mylabel + "*" + "{:5.1e}".format(1/factor)
                else:
                    factor = 1.0
            else:
                factor = 1.0

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


            if ((type(lolims) == bool ) & (type(uplims) == bool)):
                lims = False
            elif (type(lolims) == bool):
                lims = uplims
            elif (type(uplims) == bool):
                lims = lolims
            else:
                lims = uplims | lolims

            # make plot
            plt.errorbar(xvals, yvals,
                         yerr = yerrs,
                         uplims = uplims,
                         lolims = lolims,
                         marker = mymarker, 
                         markerfacecolor = 'none',
                         markeredgecolor = mycolor, 
                         color = mycolor,
                         linestyle = '--')

            if np.any(lims):
                plt.scatter(xvals[~lims], yvals[~lims],
                            marker = mymarker,
                            color = mycolor,
                            label = mylabel)
            else:
                plt.scatter(xvals, yvals,
                            marker = mymarker,
                            color = mycolor,
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

#----------------------------------------------------------------------

def compare_stacks_line(stack_list, stack_name_list, 
                        empire_db = None, ### NEED TO FIX UP OLD SCRIPTS TO POINT TO DB.
                        outdir = None,
                        line='HCN',
                        galaxy='NGC2903',
                        dist=None,
                        ylim=None,
                        xlim=None):
    
    nquant = len(stack_list)
    markers = ['o','^','s','*','D','>'] # 5 items
    markerlist = np.tile(markers,int(np.ceil(nquant/len(markers))))
    markerlist = markerlist[0:nquant]
    

    # compare results of different stacks
    for (stack,stack_name,marker) in zip(stack_list, stack_name_list,markerlist):
        idx = (stack['galaxy'] == galaxy) & (stack['bin_type'] == 'radius')

        if dist:
            binval = dist.to('kpc') * stack[idx]['bin_upper']/206265.0
            #binval = dist.to('kpc') * stack[idx]['bin_lower']/206265.0
        else:
            binval =  stack[idx]['bin_upper']

        binval = binval.value

        if line in ['HCN','HCOp','13CO','C18O','CO']:
            quant = 'int_intensity_sum_'+line
            quant_err = 'int_intensity_sum_'+line+'_err'
        elif line in stack.columns:
            quant = line
            if quant+'_err' in stack.columns:
                quant_err = quant+'_err'
            else:
                quant_err = ''
        else:
            print("The quantity " + line + " not in data base. Returning.\n")
            return

        val = stack[idx][quant].value
        if quant_err:
            val_err = stack[idx][quant_err].value
            plt.errorbar(binval, val, yerr=val_err, marker=marker,label=stack_name)
        else:
            plt.plot(binval, val, marker=marker,label=stack_name)

    # add empire
    if empire_db:
        if line in ['HCN','HCOp']:
            idx = empire_db['ID'] == galaxy.lower()
            plt.errorbar(empire_db[idx]['Rad'],empire_db[idx]['I'+line.replace('HCOp','HCO+')],
                         yerr = empire_db[idx]['e_I'+line.replace('HCOp','HCO+')],
                         marker='+',label='EMPIRE (pub)',color='purple')
        else:
            print("Value " + line + " not in empire database.\n")

    if ylim:
        plt.ylim(ylim)
        

    if xlim:
        plt.xlim(xlim)

    plt.yscale('log')
    plt.legend()
    plt.title(galaxy + ' - ' + line.replace('HCOp','HCO+'),size='large')
    plt.xlabel('radius (kpc)',size='medium')
    plt.ylabel('I (K km/s)',size='medium')

    if outdir == None:
        outDir = './'
    
    plt.savefig(os.path.join(outdir,galaxy + '_' + line + '_stack_compare.png'))
    plt.savefig(os.path.join(outdir,galaxy + '_' + line + '_stack_compare.pdf'))

    plt.close()

#----------------------------------------------------------------------



def ratio_plots(degas, stack, ydata='ratio_HCN_CO', xdata='r25', 
                add_empire=True, 
                include_uplims=False,
                release='DR1',
                pltname=None,
                pltlabel=None):

    '''

    Purpose: general purpose plotting routine for sfe_dense vs. x and
    fdense vs. x 

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    7/21/2022   A.A. Kepley     Original Code

    '''

    # select release galaxies. TODO -- offer all option?
    print("Only plotting release " + release + "\n")
    release = degas[release] == 1
    nrelease = np.sum(release)

    # set up x labels
    if xdata == 'r25':
        xlabel = r'R/R$_{25}$'
    elif xdata == 'mstar':
        xlabel = r'$\Sigma_{*}$ [M$_\odot$ pc$^{-2}$]'
    elif xdata == 'ICO':
        xlabel = r'$\Sigma_{mol}$ [M$_\odot$ pc$^{-2}$]'
    else:
        print("no valid xlabel value using xdata name: " + xdata)
        xlabel = xdata
    
    # set up y labels
    if ydata == 'ratio_HCN_CO':
        ylabel = r'HCN-to-CO'
    elif ydata == 'ratio_ltir_mean_HCN':
        ylabel = r'L$_{TIR}$ / HCN'
    elif ydata == 'ratio_HCOp_CO':
        ylabel = r'HCO+-to-CO'
    elif ydata == 'ratio_ltir_mean_HCOp':
        ylabel = r'L$_{TIR}$ / HCO+'
    else:
        ylabel = ydata

  
    # setup plot style
    markers = ['o','v','^','s','>','D'] # 6 items
    colors = ['royalblue','forestgreen','darkorange','royalblue','crimson','rebeccapurple','darkcyan','darkmagenta']

    markerlist = np.tile(markers,int(np.ceil(nrelease/len(markers))))
    markerlist = markerlist[0:nrelease]

    colorlist = np.tile(colors,int(np.ceil(nrelease/len(colors))))
    colorlist = colorlist[0:nrelease]
    
    # set up plot
    fig = plt.figure(figsize=(8,6),facecolor='white',edgecolor='white')
    fig.subplots_adjust(left=0.15,right=0.78,bottom=0.15, top=0.9)

    ax = fig.add_subplot(1,1,1)
   
    # for each dr1 galaxy, show radial trends for each line.
    for (galaxy,color,marker) in zip(degas[release],colorlist,markerlist):
        idx = ( (stack['galaxy'] == galaxy['NAME']) \
                & (stack['bin_type'] == xdata))

        if xdata == 'ICO':
            # apply alpha_co and correct for inclination
            alpha_co = float(stack.meta['ALPHA_CO'].split(' ')[0])
            factor = alpha_co * np.cos(np.radians(galaxy['INCL_DEG']))
        else:
            factor = 1.0

        xvals = stack[idx]['bin_mean'] * factor



        yvals = stack[idx][ydata] 
        yvals_err = stack[idx][ydata+'_err'] 

        if ydata+'_uplim' in stack.columns:
            lims = stack[idx][ydata+'_uplim']
        elif ydata + '_lolim' in stack.columns:
            lims = stack[idx][ydata+'_lolim']
            yvals_err[lims] = yvals[lims] * 0.3 

        ax.errorbar(xvals[~lims], yvals[~lims],
                    yerr = yvals_err[~lims],
                    marker = marker,
                    linestyle= '--',
                    color=color,
                    label=galaxy['NAME'],zorder=10)

    ax.set_yscale('log')
    if xdata != 'r25':
        ax.set_xscale('log')
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylabel(ylabel,fontsize=20)
    ax.tick_params(axis='both',labelsize=16)
    
    if pltlabel:
        ax.set_title(pltlabel)

    # empire fit parameters
    if add_empire:
        xmin, xmax = ax.get_xlim()
        if ydata == 'ratio_HCN_CO':
            if xdata == 'r25':
                xvals_empire = np.linspace(xmin, xmax)
                yvals_empire = 10**(-1.5-0.8*xvals_empire)
            elif xdata == 'mstar':
                xvals_empire = np.logspace(np.log10(xmin),np.log10(xmax))
                yvals_empire = 10**(-2.7+0.4*np.log10(xvals_empire))
            elif xdata == 'ICO':
                xvals_empire = np.logspace(np.log10(xmin),np.log10(xmax))
                yvals_empire = 10**(-2.4 + 0.5*np.log10(xvals_empire))
            else:
                print("No empire equation available for " + xdata + ' vs. ' + ydata)
        elif ydata == 'ratio_ltir_mean_HCN':
            if xdata == 'r25': 
                xvals_empire = np.linspace(xmin,xmax)
                yvals_empire = 10**(2.8+0.6*xvals_empire)
            elif xdata == 'mstar':
                xvals_empire = np.logspace(np.log10(xmin),np.log10(xmax))
                yvals_empire = 10**(4.0-0.4*np.log10(xvals_empire))
            elif xdata == 'ICO':
                xvals_empire = np.logspace(np.log10(xmin),np.log10(xmax))
                yvals_empire = 10**(3.5-0.4*np.log10(xvals_empire))
            else:
                print("No empire equation available for " + xdata + ' vs. ' + ydata)
        else:
            print("Not a validate y data set for EMPIRE: " + ydata)
                                
        if xdata != 'r25':
            ax.loglog(xvals_empire,yvals_empire,color='gray',linestyle=':',linewidth=5, label='EMPIRE',zorder=1)
        else:
            ax.semilogy(xvals_empire,yvals_empire,color='gray',linestyle=':',linewidth=5, label='EMPIRE',zorder=1)


    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # remove errorbars
    new_handles = [h[0] for h in handles[1:-1]]
    new_handles.insert(0,handles[0])


    mylegend = ax.legend(new_handles, labels, loc='upper left',
                         bbox_to_anchor=(1.0,1.0),handlelength=3, 
                         borderpad=1.2)    

    fig.show()

    if pltname:
        # save as both png and pdf
        fig.savefig(pltname+'.png')
        fig.savefig(pltname+'.pdf')


def plot_correlation_coeffs(results, 
                            corr_quant='ratio_HCN_CO',
                            pltname=None):

    '''
    plot the correlation coefficients for each galaxy and the overall
    sample for correlation bins and 

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    8/4/2022    A.A. Kepley     Original code

    '''
    corr_bins = ['r25','mstar','ICO']
 
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,8),
                                      edgecolor='white')
    fig.subplots_adjust(hspace=0.2, bottom=0.12,left=0.15,top=0.95)
    

    for (mybin,myax) in zip(corr_bins,[ax1,ax2,ax3]):

        idx_all = (results['galaxy'] == 'all') & (results['corr_bin'] == mybin) & (results['quant'] == corr_quant)
        idx_gal = (results['galaxy'] != 'all') & (results['corr_bin'] == mybin) & (results['quant'] == corr_quant)

        myax.set_axisbelow(True)
        myax.grid(axis='x')

        # plot correlation coefficients
        myax.scatter(results['galaxy'][idx_gal],results['corr_val'][idx_gal],
                     facecolor='gray',s=20,edgecolor='None',marker='o',
                     zorder=1, label="n < 3")

        # mark those with >3 values
        large_nval = results[idx_gal]['nval'] >= 3
        myax.scatter(results[idx_gal][large_nval]['galaxy'], 
                     results[idx_gal][large_nval]['corr_val'],
                     marker='o',facecolor="orange",edgecolor='None',
                     s=30,
                     label="n >= 3",
                     zorder=1)

        # mark the significant values
        good_pval = results[idx_gal]['p_val'] < 0.05
        myax.scatter(results[idx_gal][good_pval]['galaxy'], 
                     results[idx_gal][good_pval]['corr_val'],
                     marker='o',facecolor="None",edgecolor='black',
                     s=40,
                     label='p<0.05',
                     zorder=1)
        
        # plot the average value
        myax.axhline(results[idx_all]['corr_val'],color='gray',linewidth=3,
                     zorder=2)
        myax.axhline(0,color='gray',linewidth=1)

        myax.set_ylim(-1.2,1.2)
        myax.set_ylabel(r'$\rho$')

        if mybin == 'r25':
            mybinstr = r"R/R$_{25}$"
        elif mybin == 'mstar':
            #mybinstr = r"$\Sigma_*$ [M$_\odot$ pc$^{-2}$]"
            mybinstr = r"$\Sigma_*$"
        elif mybin == 'ICO':
            #mybinstr = r"I$_{CO}$ [K km s$^{-1}$]"
            mybinstr = r"I$_{CO}$"
        else:
            print("Bin " + mybin + " not recognized. Using bin value as text.")
        
        if corr_quant == 'ratio_ltir_mean_HCN':
            #corr_quant_str = r"L$_{TIR}$/HCN [L$_\odot$ pc$^{-2}$ (K km s$^{-1}$)$^{-1}$]"
            corr_quant_str = r"L$_{TIR}$/HCN"
        elif corr_quant == 'ratio_HCN_CO':
            corr_quant_str = r"HCN-to-CO"
        else:
            print("Correlation Quantity " + corr_quant + " not recognized. Using corr_quant value as text")

        myax.set_title(corr_quant_str + ' - ' + mybinstr) 

    if corr_quant == 'ratio_HCN_CO':
        myloc = 'lower right'
    elif corr_quant == 'ratio_ltir_mean_HCN':
        myloc = 'upper right'

    ax3.legend(loc=myloc)
    ax3.tick_params(axis='x',which='both',labelrotation=90)

    if pltname:
        fig.savefig(pltname+'.png')
        fig.savefig(pltname+'.pdf')
        
        
