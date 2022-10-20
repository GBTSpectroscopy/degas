from astropy.table import Table
import os
from degas.analysis_plot import plot_fit_coeff_hist

stack_fit = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_fits.fits')
stack_fit_noNGC4414 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_fits_noNGC4414.fits')

stack_fit_pergal = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_fits_pergal.fits')

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends'

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_HCN_CO',pltname=os.path.join(plotdir,'ratio_HCN_CO_fit_hist'),pltlabel=r'constant R$_{21}$')


plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_ltir_mean_HCN',pltname=os.path.join(plotdir,'ratio_ltir_mean_HCN_fit_hist'),pltlabel=r'constant R$_{21}$')

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_HCOp_CO',pltname=os.path.join(plotdir,'ratio_HCOp_CO_fit_hist'),pltlabel=r'constant R$_{21}$')

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_ltir_mean_HCOp',pltname=os.path.join(plotdir,'ratio_ltir_mean_HCOp_fit_hist'),pltlabel=r'constant R$_{21}$')

### R21 spatial fits here.

stack_fit = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits.fits')
stack_fit_noNGC4414 = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits_noNGC4414.fits')

stack_fit_pergal = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21/stack_fits_pergal.fits')

plotdir = '/lustre/cv/users/akepley/degas/sfe_fdense_trends'

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_HCN_CO',pltname=os.path.join(plotdir,'ratio_HCN_CO_fit_hist_spatialR21'),pltlabel=r'spatial R$_{21}$')

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_ltir_mean_HCN',pltname=os.path.join(plotdir,'ratio_ltir_mean_HCN_fit_hist_spatialR21'),pltlabel=r'spatial R$_{21}$')

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_HCOp_CO',pltname=os.path.join(plotdir,'ratio_HCOp_CO_fit_hist_spatialR21'),pltlabel=r'spatial R$_{21}$')

plot_fit_coeff_hist(stack_fit_pergal,stack_fit,alt_r25=stack_fit_noNGC4414,fit_quant='ratio_ltir_mean_HCOp',pltname=os.path.join(plotdir,'ratio_ltir_mean_HCOp_fit_hist_spatialR21'),pltlabel=r'spatial R$_{21}$')

