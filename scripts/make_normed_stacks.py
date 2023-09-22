import os
from astropy.table import Table
from degas.analysis_stack import normalizeStacks

release =  'IR6p1'

scriptDir = os.environ['SCRIPTDIR']
stackDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_spatialR21')

degas_db = Table.read(os.path.join(scriptDir,'degas_base.fits'))
stack = Table.read(os.path.join(stackDir,'stack_IR6p1_spatialR21_mom1_pruned.fits'))
outfile = os.path.join(stackDir,'stack_IR6p1_spatialR21_mom1_pruned_norm.fits')

stack_norm = normalizeStacks(stack,degas_db,outfile=outfile)

# -----------------

# # plotting  code to check results
# idx = (test['galaxy'] == 'NGC4321') & (test['bin_type'] == 'mstar')
# #plt.loglog(test[idx]['bin_mean'],test[idx]['ratio_HCN_CO'])
# plt.loglog(test[idx]['bin_mean_norm'],test[idx]['ratio_HCN_CO_norm'])

# idx = (test['galaxy'] == 'NGC2903') & (test['bin_type'] == 'mstar')
# plt.loglog(test[idx]['bin_mean_norm'],test[idx]['ratio_HCN_CO_norm'])


# idx = (test['galaxy'] == 'NGC4321') & (test['bin_type'] == 'mstar')
# #plt.loglog(test[idx]['bin_mean'],test[idx]['ratio_ltir_mean_HCN'])
# plt.loglog(test[idx]['bin_mean_norm'],test[idx]['ratio_ltir_mean_HCN_norm'])
# idx = (test['galaxy'] == 'NGC2903') & (test['bin_type'] == 'mstar')
# plt.loglog(test[idx]['bin_mean_norm'],test[idx]['ratio_ltir_mean_HCN_norm'])


