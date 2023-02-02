import os
from astropy.table import Table
from degas.analysis_plot import plot_moments

indir = '/lustre/cv/users/akepley/degas/moments_IR6p1'
outdir = '/lustre/cv/users/akepley/degas/moments_IR6p1/figures/'

degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
degas_dr1 = degas[degas['DR1'] == 1]
#degas_dr1 = degas[degas['NAME'] == 'NGC2903']


gal_list = degas_dr1['NAME'] #import list of DEGAS targets

line_list = ['13CO', 'C18O', 'HCN', 'HCOp']

if not os.path.exists(outdir):
    os.mkdir(outdir)

for gal in gal_list:
    for line in line_list:
        plot_moments(gal, line, indir=indir, outdir=outdir,
                     moments=[0], masked=True)

        plot_moments(gal, line, indir=indir, outdir=outdir,
                     moments=[1], masked=False)
