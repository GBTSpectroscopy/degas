import os
from astropy.table import Table
import glob
from degas.analysis_plot import plot_galaxy_overview
import pickle
import csv
import itertools


mom_dir = os.path.join(os.environ['ANALYSISDIR'], 'moments_IR6p1')
regrid_dir = os.path.join(os.environ['ANALYSISDIR'], 'IR6p1_regrid')
out_dir = os.path.join(os.environ['ANALYSISDIR'], 'overview_plots_IR6p1')

levels_file = os.path.join(out_dir,'levels.pkl')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
degas_dr1 = degas[degas['DR1'] == 1]
#degas_dr1 = degas[degas['NAME'] == 'NGC4321']
#degas_dr1 = degas[degas['NAME'] == 'NGC2903']



if not os.path.exists(out_dir):
    os.mkdir(out_dir)

base_value_dict = {}

for galaxy in degas_dr1:
    hcn_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCN_*_mom0_masked.fits'))[0]
    hcn_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCN_*_mom0.fits'))[0]

    hcn_err_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCN_*_emom0.fits'))[0]

    hcop_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCOp_*_mom0_masked.fits'))[0]
    hcop_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCOp_*_mom0.fits'))[0]

    hcop_err_file = glob.glob(os.path.join(mom_dir,galaxy['NAME']+'_HCOp_*_emom0.fits'))[0]

    # DO I WANT TO DO NATIVE FILE OR TO SHOW RESCALED WITH RADIUS?
    co_file = glob.glob(os.path.join(regrid_dir,galaxy['NAME']+'_12CO10_*_sigmaSFR_mom0*.fits'))
    
    ## TODO: TEST THAT THIS WORKS WITH NATIVE 12CO
    if not co_file:
         co_file = glob.glob(os.path.join(regrid_dir,galaxy['NAME']+'_12CO10_*mom0*.fits'))
    co_file = co_file[0]

    ir_file = glob.glob(os.path.join(regrid_dir,galaxy['NAME']+'_LTIR_gauss15_regrid.fits'))[0]

    # set up dictionary for level values
    base_value_dict[galaxy['NAME']] = {}                      

    ## HCN plots
    tmpbase = plot_galaxy_overview(galaxy,
                         hcn_file, hcn_err_file,
                         co_file,
                         plot_title = '$^{12}$CO(1-0)\nHCN',
                         out_file = galaxy['NAME'] + '_overview_HCN_CO.pdf',
                         out_dir = out_dir)                        
    base_value_dict[galaxy['NAME']]['HCN'] = tmpbase
    
    plot_galaxy_overview(galaxy,
                         hcn_file, hcn_err_file,
                         ir_file,
                         plot_title='TIR\nHCN',
                         out_file = galaxy['NAME'] + '_overview_HCN_IR.pdf',
                         out_dir = out_dir)

    ## HCOp plots
    tmpbase = plot_galaxy_overview(galaxy,
                         hcop_file, hcop_err_file,
                         co_file,
                         plot_title = '$^{12}$CO(1-0)\nHCO+',
                         out_file = galaxy['NAME'] + '_overview_HCOp_CO.pdf',
                         out_dir = out_dir)

    base_value_dict[galaxy['NAME']]['HCOp'] = tmpbase

    plot_galaxy_overview(galaxy,
                         hcop_file, hcop_err_file,
                         ir_file,
                         plot_title='TIR\nHCO+',
                         out_file = galaxy['NAME'] + '_overview_HCOp_IR.pdf',
                         out_dir = out_dir)


fields = ['galaxy','HCN','HCOp']

csvfile = levels_file.replace('.pkl','.csv')
with open(csvfile,"w",newline='') as f:
    w = csv.DictWriter(f,fields)
    w.writeheader()
    for k in base_value_dict:
        w.writerow({field: base_value_dict[k].get(field) or k for field in fields})

# save the levels to a pickle file. Might be best to convert to csv,
# but pickle works for now.
pickle.dump(base_value_dict,open(levels_file,'wb'))

