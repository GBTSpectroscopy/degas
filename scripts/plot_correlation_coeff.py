import os
from astropy.table import Table
from degas.analysis_plot import plot_correlation_coeffs

## simple R21

coeffs = Table.read(os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','correlation_coeffs_IR6p1.fits'))

plot_correlation_coeffs(coeffs,corr_quant='ratio_HCN_CO',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','fdense_correlation_coeffs'),corr_type='kendall')

plot_correlation_coeffs(coeffs,corr_quant='ratio_ltir_mean_HCN',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','sfe_dense_correlation_coeffs'),corr_type='kendall')

# spatial R21

coeffs2 = Table.read(os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','correlation_coeffs_IR6p1_spatialR21.fits'))

plot_correlation_coeffs(coeffs2,corr_quant='ratio_HCN_CO',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','fdense_correlation_coeffs_spatialR21'),corr_type='kendall')

plot_correlation_coeffs(coeffs2,corr_quant='ratio_ltir_mean_HCN',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','sfe_dense_correlation_coeffs_spatialR21'),corr_type='kendall')

## 33as, simple R21
coeffs3 = Table.read(os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','correlation_coeffs_IR6p1_33as.fits'))

plot_correlation_coeffs(coeffs3,corr_quant='ratio_HCN_CO',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','fdense_correlation_coeffs_33as'),corr_type='kendall')

plot_correlation_coeffs(coeffs3,corr_quant='ratio_ltir_mean_HCN',pltname=os.path.join(os.environ['ANALYSISDIR'],'sfe_fdense_trends','sfe_dense_correlation_coeffs_33as'),corr_type='kendall')
