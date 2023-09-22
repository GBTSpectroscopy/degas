import os
from astropy.table import Table
from degas.analysis_fit import calc_correlation_coeffs

degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

## simple R21

stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1','stack_IR6p1_mom1_pruned.fits'))

myresult = calc_correlation_coeffs(stack,degas_db,
                                   corr_bins = ['r25','mstar','molgas','PDE'],
                                   corr_quants=['ratio_HCN_CO','ratio_ltir_mean_HCN',
                                                'ratio_HCN_13CO'], 
                                   outfile=os.path.join(os.environ['ANALYSISDIR'],'analysis','sfe_fdense_trends','correlation_coeffs_IR6p1.fits'))

## spatial R21

stack2 = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1_spatialR21','stack_IR6p1_spatialR21_mom1_pruned.fits'))

myresult = calc_correlation_coeffs(stack2,degas_db,
                                   corr_bins = ['r25','mstar','molgas','PDE'],
                                   corr_quants=['ratio_HCN_CO','ratio_ltir_mean_HCN',
                                                'ratio_HCN_13CO'], 
                                   outfile=os.path.join(os.environ['ANALYSISDIR'],'analysis','sfe_fdense_trends','correlation_coeffs_IR6p1_spatialR21.fits'))


## 33as, simple r21
stack3 = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1_33as','stack_IR6p1_33as_mom1_pruned.fits'))

myresult = calc_correlation_coeffs(stack3,degas_db,outfile=os.path.join(os.environ['ANALYSISDIR'],'analysis','sfe_fdense_trends','correlation_coeffs_IR6p1_33as.fits'))

## spatial R21 -- normed

stack4 = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1_spatialR21','stack_IR6p1_spatialR21_mom1_pruned_norm.fits'))

myresult = calc_correlation_coeffs(stack4,degas_db,norm=True,
                                   corr_quants = ['ratio_HCN_CO_norm',
                                                  'ratio_ltir_mean_HCN_norm'],
                                   outfile=os.path.join(os.environ['ANALYSISDIR'],'analysis','sfe_fdense_trends_normalized','correlation_coeffs_IR6p1_norm.fits'))
