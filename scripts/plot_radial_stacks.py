import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

# setup information sources
stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_test','stack_IR6p1_mom1.fits'))
plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots','radial_stacks')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

if not os.path.exists(plotDir):
    os.mkdir(plotDir)

# plot setup
style = {'CO': {'marker':'o','color':'orange','name':'CO'}, 
         'HCN': {'marker':'o','color':'green','name':'HCN'}, 
         'HCOp': {'marker':'o','color':'blue','name':'HCO+'}, 
         '13CO': {'marker':'o','color':'red', 'name':'13CO'}, 
         'C18O': {'marker':'o','color': 'magenta','name':'C18O'}}

# only look at dr1 galaxies
dr1 = degas['DR1'] == 1

# For each dr1 galaxy, show radial trends for each line.
for galaxy in degas[dr1]:

    plt.close()

    idx = ( (stack['galaxy'] == galaxy['NAME']) \
            & (stack['bin_type'] == 'radius'))
    
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
            else:
                factor = 1.0
        
            flatax[i].plot(stack[idx]['spectral_axis'][i],
                           stack[idx]['stack_profile_'+line][i] / factor,
                           color=style[line]['color'],
                           linestyle='-',
                           marker=None,
                           label=style[line]['name'],
                           linewidth=1)
            
            flatax[i].axhline(0,linestyle=':',linewidth=1,color='gray')
            radius_label = 'r = {:5.2f}' .format(stack[idx]['bin_mean'][i])+' arcsec'
            flatax[i].text(0.1, 0.9, radius_label,
                           transform=flatax[i].transAxes,fontsize='x-small')
            flatax[i].legend(loc='upper right',fontsize='x-small')

            
            

    for ax in flatax:
        ax.set(xlabel=r'Velocity (km s$^{-1}$)',
               ylabel=r'Stacked Mean Intensity (K)')
    
    for ax in flatax:
        ax.label_outer()

    


    plt.show()
    plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_radial_stacks.pdf'))
