import matplotlib.pyplot as plt
import os
from astropy.table import Table

# setup information sources
degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
stack = Table.read('/lustre/cv/users/akepley/degas/stack_test/stack_IR6p0_mom1.fits')


plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots','fdense_plots')

if not os.path.exists(plotDir):
    os.mkdir(plotDir)

# only look at dr1 galaxies
dr1 = degas['DR1'] == 1

# for each dr1 galaxy, show radial trends for each line.
for galaxy in degas[dr1]:
    
    plt.close()

    idx = ( (stack['galaxy'] == galaxy['NAME']) \
            & (stack['bin_type'] == 'stellarmass'))

    mstar = stack[idx]['bin_mean']

    lolims = stack[idx]['ratio_ltir_mean_HCN_lolim']

    sfe_dense = stack[idx]['ratio_ltir_mean_HCN']
    sfe_dense_err = stack[idx]['ratio_ltir_mean_HCN_err']
    sfe_dense_err[lolims] = sfe_dense[lolims] * 0.3 

    plt.errorbar(mstar, sfe_dense,
                 yerr = sfe_dense_err,
                 uplims = lolims,
                 marker = 'o',
                 markerfacecolor='none',
                 markeredgecolor = 'blue',
                 linestyle= '--')


    plt.scatter(mstar[~lolims], sfe_dense[~lolims],
                marker='o',
                color='blue')

    
    plt.yscale('log')
    plt.xscale('log')

    plt.title(galaxy['NAME'])
    plt.xlabel(r'$\Sigma_{*}$ (M$_\odot$ pc$^{-2}$)')
    plt.ylabel(r'log$_{10}$ (IR-to-HCN)')

    plt.show()
    plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_sfe_dense_vs_mstar.pdf'))    
    plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_sfe_dense_vs_mstar.png'))    

    
    
