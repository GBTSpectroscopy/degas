import os
from astropy.table import Table
import matplotlib.pyplot as plt

# setup information sources
stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_test','test_mom1.fits'))
plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots','radial_plots')

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

# for each dr1 galaxy, show radial trends for each line.
for galaxy in degas[dr1]:

    plt.close()

    
    for line in ['CO','HCN','HCOp','13CO','C18O']:

        if (galaxy['NAME'] == 'NGC6946') & ((line == '13CO') | (line == 'C18O')):
            continue

        idx = ( (stack['galaxy'] == galaxy['NAME']) \
                & (stack['bin_type'] == 'radius') \
                & (stack[line+'_stack_sum'] > 0))
        
        # get radius in kpc -- radius stacks in arcsec.
        radius = (stack[idx]['bin_mean'] * galaxy['DIST_MPC'] * 1e6 / 206265.0) / 1e3
        # determine whether upper or lower limits
        uplims = stack[idx][line+'_stack_noise']*3.0 > stack[idx][line+'_stack_sum']
   
        stack_sum = stack[idx][line+'_stack_sum']
        stack_sum[uplims] = stack[idx][line+'_stack_noise'][uplims]*3.0

        yerr = stack[idx][line+'_stack_noise']
        yerr[uplims] = stack_sum[uplims] * 0.3


        plt.errorbar(radius, stack_sum,
                     yerr = yerr,
                     uplims = uplims,
                     marker = style[line]['marker'],
                     color = style[line]['color'], 
                     linestyle = '--', 
                     label = style[line]['name'])


    plt.yscale('log')
   
    plt.legend()
    plt.title(galaxy['NAME'])
    plt.xlabel(r"Galactocentric Radius (kpc)") 
    plt.ylabel(r"Stacked Integrated Intensity (K km s$^{-1}$)") 
    plt.show()

    plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_radial_lines.pdf'))
