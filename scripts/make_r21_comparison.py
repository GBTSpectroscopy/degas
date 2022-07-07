import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

from degas.analysis_plot import compare_stacks_line

stackDir = '/lustre/cv/users/akepley/degas/stack_IR6p1'
r21_simple = Table.read(os.path.join(stackDir,"stack_IR6p1_mom1_pruned.fits"))

stackDir = '/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21'
r21_spatial = Table.read(os.path.join(stackDir,"stack_IR6p1_spatialR21_mom1_pruned.fits"))

degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

plot_dir = os.path.join(os.environ['ANALYSISDIR'],'R21_simple_spatial_comparison') 

dr1 = degas_db['DR1'] == 1

for galaxy in degas_db[dr1]:
    compare_stacks_line([r21_simple,r21_spatial],['R21_simple','R21_spatial'],
                        line='CO',galaxy=galaxy['NAME'],
                        outdir=plot_dir)
