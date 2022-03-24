from degas.analysis_plot import compare_stacks_line
from astropy.table import Table
import astropy.units as u

# setup information sources
degas = Table.read('/lustre/cv/users/akepley/degas/code/degas/scripts/degas_base.fits')
empire_db = Table.read('/lustre/cv/users/akepley/degas/ancillary_data/empire/apjab2b95t10_mrt.txt',format='ascii.cds')

# set up stacks
stack = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1/stack_IR6p1_mom1_pruned.fits')
stack_empire = Table.read('/lustre/cv/users/akepley/degas/stack_empire/stack_empire_mom1.fits')
stack_33as = Table.read('/lustre/cv/users/akepley/degas/stack_IR6p1_33as/stack_IR6p1_33as_mom1_pruned.fits')

stack_list = [stack,  stack_33as, stack_empire]
stack_name_list = ['DEGAS', 'DEGAS (33 arcsec)',  'EMPIRE (DEGAS CODE)']

overlaplist = ['NGC2903','NGC4321','NGC5055','NGC6946']
linelist = ['HCN','HCOp','13CO','C18O']


for galaxy in overlaplist:

    dist = degas[degas['NAME'] == galaxy]['DIST_MPC'][0] * u.Mpc

    for line in linelist:

        compare_stacks_line(stack_list, stack_name_list, empire_db,
                            outdir = '/lustre/cv/users/akepley/degas/empire_degas_comp',
                            line = line,
                            galaxy = galaxy,
                            dist = dist)



