import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

from degas.analysis_plot import plot_stack

# setup information sources
stack = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1','stack_IR6p1_mom1.fits'))
stack_pruned = Table.read(os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1','stack_IR6p1_mom1_pruned.fits'))
degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1','stack_plots')

# Plot various stacks
plot_stack(stack, 'radius',
           plot_dir,
           degas_db, release='DR1')

plot_stack(stack, 'r25',
           plot_dir,
           degas_db, release='DR1')

plot_stack(stack, 'mstar',
           plot_dir,
           degas_db, release='DR1')

plot_stack(stack, 'ICO',
           plot_dir,
           degas_db, release='DR1')


plot_dir =  os.path.join(os.environ['ANALYSISDIR'],'stack_IR6p1','stack_plots_pruned')

plot_stack(stack_pruned, 'mstar',
           plot_dir,
           degas_db, release='DR1')

plot_stack(stack_pruned, 'ICO',
           plot_dir,
           degas_db, release='DR1')
