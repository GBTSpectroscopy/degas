import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np


def plot_stack(stack, bin_type, plot_dir, degas_db, release='DR1'):
    '''
    general purpose routine to plot stacks.

    stack: astropy Table with stack

    bin_type: type of bin to plot

    plot_dir: plotting directory

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
             'HCN': {'marker':'o','color':'green','name':'HCN'}, 
             'HCOp': {'marker':'o','color':'blue','name':'HCO+'}, 
             '13CO': {'marker':'o','color':'red', 'name':'13CO'}, 
             'C18O': {'marker':'o','color': 'magenta','name':'C18O'}}
    
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

    
