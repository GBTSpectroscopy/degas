import os
from astropy.table import Table
import glob
from degas.analysis_plot import plot_galaxy_overview

mom_dir = os.path.join(os.environ['ANALYSISDIR'], 'moments_IR6p1')
regrid_dir = os.path.join(os.environ['ANALYSISDIR'], 'IR6p1_regrid')
out_dir = os.path.join(os.environ['ANALYSISDIR'], 'overview_plots_IR6p1')

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
degas_dr1 = degas[degas['DR1'] == 1]
#degas_dr1 = degas[degas['NAME'] == 'NGC4321']
#degas_dr1 = degas[degas['NAME'] == 'NGC2903']

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


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
                                      
    ## HCN plots
    plot_galaxy_overview(galaxy,
                         hcn_file, hcn_err_file,
                         co_file,
                         out_file = galaxy['NAME'] + '_overview_HCN_CO.pdf',
                         out_dir = out_dir)                        
    
    plot_galaxy_overview(galaxy,
                         hcn_file, hcn_err_file,
                         ir_file,
                         out_file = galaxy['NAME'] + '_overview_HCN_IR.pdf',
                         out_dir = out_dir)

    ## HCOp plots
    plot_galaxy_overview(galaxy,
                         hcop_file, hcop_err_file,
                         co_file,
                         out_file = galaxy['NAME'] + '_overview_HCOp_CO.pdf',
                         out_dir = out_dir)

    plot_galaxy_overview(galaxy,
                         hcop_file, hcop_err_file,
                         ir_file,
                         out_file = galaxy['NAME'] + '_overview_HCOp_IR.pdf',
                         out_dir = out_dir)
