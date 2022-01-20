import os
from astropy.table import Table
import matplotlib.pyplot as plt

# setup information sources
stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_test','stack_IR6p1_mom1.fits'))
plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots','radial_plots')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

if not os.path.exists(plotDir):
    os.mkdir(plotDir)

# plot setup
style = {'CO': {'marker':'o','color':'orange','name':'CO'}, 
         'HCN': {'marker':'o','color':'green','name':'HCN'}, 
         'HCOp': {'marker':'o','color':'blue','name':'HCO+'}, 
         '13CO': {'marker':'o','color':'red', 'name':'13CO'}, 
         'C18O': {'marker':'o','color': 'magenta','name':'C18O'},
         'comass_mean': {'marker':'o','color':'orange','name':r'$\Sigma_{CO}$ (M$_\odot/pc^2$)'},
         'mstar_mean': {'marker':'o','color':'green','name':r'$\Sigma_*$ (M$_\odot/pc^2$)'},
         'sfr_mean':{'marker':'o','color':'blue','name':r'$\Sigma_{SFR}$ (M$_\odot/yr/pc^2$)'},
         'ratio_HCN_CO': {'marker':'o','color':'green','name':'HCN/$^{12}$CO'},
         'ratio_HCOp_CO': {'marker':'o','color':'blue','name':'HCO+/$^{12}$CO'}, 
         'ratio_13CO_CO': {'marker':'o','color':'red','name':'$^{13}$CO/$^{12}$CO'},
         'ratio_C18O_CO': {'marker':'o','color':'orange','name':'$C^{18}$O/$^{12}$CO'},
         'ratio_HCOp_HCN': {'marker':'o','color':'green','name':'HCO+/HCN'}}
         

# only look at dr1 galaxies
dr1 = degas['DR1'] == 1

# turn on different plotting
doline = True
doother = True
docoratio = True
doratio = True

# for each dr1 galaxy, show radial trends for each line.
for galaxy in degas[dr1]:

    plt.close()

    idx = ( (stack['galaxy'] == galaxy['NAME']) \
            & (stack['bin_type'] == 'radius'))

    # get radius in kpc -- radius stacks in arcsec.
    radius = (stack[idx]['bin_mean'] * galaxy['DIST_MPC'] * 1e6 / 206265.0) / 1e3

    if doline:
    
        for line in ['CO','HCN','HCOp','13CO','C18O']:

            if (galaxy['NAME'] == 'NGC6946') & ((line == '13CO') | (line == 'C18O')):
                continue

            # determine whether upper or lower limits
            uplims = stack[idx]['int_intensity_sum_uplim_'+line]

            int_intensity = stack[idx]['int_intensity_sum_'+line]

            yerr = stack[idx]['int_intensity_sum_err_'+line]
            yerr[uplims] = int_intensity[uplims] * 0.3


            plt.errorbar(radius, int_intensity,
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

        plt.close()

    if doother:

        for otherdata in ['comass_mean', 'mstar_mean',  'sfr_mean']:

            if otherdata == 'sfr_mean':
                factor = 1e10
            else:
                factor = 1.0

            plt.plot(radius, stack[idx][otherdata]*factor,
                     marker = style[otherdata]['marker'],
                     color = style[otherdata]['color'],
                     linestyle = '--',
                     label = style[otherdata]['name'])
           
        plt.yscale('log')
        
        plt.legend()
        plt.title(galaxy['NAME'])
        plt.xlabel(r"Galactocentric Radius (kpc)") 
        plt.ylabel(r"Surface Density")
        plt.show()

        plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_radial_other.pdf'))

        plt.close()


    if docoratio:

        for coratio in ['ratio_HCN_CO','ratio_HCOp_CO','ratio_13CO_CO','ratio_C18O_CO']:
            uplims = stack[idx][coratio+'_uplim']
            
            ratio = stack[idx][coratio]

            yerr = stack[idx][coratio+'_err']
            yerr[uplims] = ratio[uplims] * 0.3


            plt.errorbar(radius, ratio,
                         yerr = yerr,
                         uplims = uplims,
                         marker = style[coratio]['marker'],
                         color = style[coratio]['color'],
                         linestyle = '--',
                         label = style[coratio]['name'])

        plt.yscale('log')
        plt.legend()
        plt.title(galaxy['NAME'])
        plt.xlabel(r"Galactocentric Radius (kpc)") 
        plt.ylabel(r"Ratio")
        plt.show()
        
        plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_radial_coratio.pdf'))
        
        plt.close()
        
    if doratio:
        
        for myratio in ['ratio_HCOp_HCN']:
            uplims = stack[idx][myratio+'_uplim']
            lolims = stack[idx][myratio+'_lolim']

            ratio = stack[idx][myratio]

            yerr = stack[idx][myratio+'_err']
            yerr[uplims] = ratio[uplims] *0.3
            yerr[lolims] = ratio[lolims] * 0.3

            plt.errorbar(radius, ratio,
                         yerr = yerr,
                         uplims = uplims,
                         lolims = lolims,
                         marker = style[myratio]['marker'],
                         color = style[myratio]['color'],
                         linestyle = '--',
                         label = style[myratio]['name'])

        plt.legend()
        plt.title(galaxy['NAME'])
        plt.xlabel(r"Galactocentric Radius (kpc)") 
        plt.ylabel(r"Ratio")
        plt.show()

        
        plt.savefig(os.path.join(plotDir,galaxy['NAME']+'_radial_ratio.pdf'))
        
        plt.close()
