import matplotlib.pyplot as plt
import os
from astropy.table import Table
import numpy as np

# setup information sources
degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
stack = Table.read('/lustre/cv/users/akepley/degas/stack_test/stack_IR6p0_mom1.fits')

plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots','fdense_plots')

if not os.path.exists(plotDir):
    os.mkdir(plotDir)

# only look at dr1 galaxies
dr1 = degas['DR1'] == 1

ndr1 = np.sum(dr1)

# setup plot style
markers = ['o','v','^','s','*','D'] # 6 items
colors = ['royalblue','forestgreen','darkorange','royalblue','crimson','rebeccapurple','darkcyan','darkmagenta']

ndr1 = np.sum(dr1)
markerlist = np.tile(markers,int(np.ceil(ndr1/len(markers))))
markerlist = markerlist[0:ndr1]

colorlist = np.tile(colors,int(np.ceil(ndr1/len(colors))))
colorlist = colorlist[0:ndr1]

# set up plot
fig = plt.figure(figsize=(8,6),facecolor='white',edgecolor='white')
fig.subplots_adjust(left=0.1,right=0.8,bottom=0.1, top=0.9)

ax = fig.add_subplot(1,1,1)

# for each dr1 galaxy, show radial trends for each line.
for (galaxy,color,marker) in zip(degas[dr1],colorlist,markerlist):
    idx = ( (stack['galaxy'] == galaxy['NAME']) \
            & (stack['bin_type'] == 'intensity'))

    # unit of molmass is solMass/pc^2
    alpha_co = float(stack.meta['ALPHA_CO'].split()[0])
    alpha_co_units = ' '.join(stack.meta['ALPHA_CO'].split()[1:])
    molmass = stack[idx]['bin_mean'] * alpha_co

    lolims = stack[idx]['ratio_ltir_mean_HCN_lolim']

    sfe_dense = stack[idx]['ratio_ltir_mean_HCN']
    sfe_dense_err = stack[idx]['ratio_ltir_mean_HCN_err']
    sfe_dense_err[lolims] = sfe_dense[lolims] * 0.3 

    ax.errorbar(molmass, sfe_dense,
                 yerr = sfe_dense_err,
                 uplims = lolims,
                 marker = marker,
                 markerfacecolor='none',
                 markeredgecolor=color,
                 linestyle= '--',
                 color=color)

    ax.scatter(molmass[~lolims], sfe_dense[~lolims],
                marker=marker,
                color=color,
                label=galaxy['NAME'])

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.legend(loc='upper left',bbox_to_anchor=(1.0,1.0))

    ax.set_xlabel(r'$\Sigma_{mol}$ (M$_{\odot}$ pc$^{-2}$)')
    ax.set_ylabel(r'log$_{10}$ (IR-to-HCN)')

    fig.show()

fig.savefig(os.path.join(plotDir,'sfe_dense_vs_molmass_combined.pdf'))    
fig.savefig(os.path.join(plotDir,'sfe_dense_vs_molmass_combined.png'))    

plt.close()

    
