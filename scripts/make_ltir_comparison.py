import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

from degas.analysis_plot import compare_stacks_line

stackDir = '/lustre/cv/users/akepley/degas/stack_IR6p1'
ltir_all = Table.read(os.path.join(stackDir,"stack_IR6p1_mom1_pruned.fits"))

stackDir = '/lustre/cv/users/akepley/degas/stack_IR6p1_L24micron'
ltir_24 = Table.read(os.path.join(stackDir,"stack_IR6p1_L24micron_mom1_pruned.fits"))

degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

plot_dir = os.path.join(os.environ['ANALYSISDIR'],'ltir_comparison') 

dr1 = degas_db['DR1'] == 1

for galaxy in degas_db[dr1]:
    compare_stacks_line([ltir_all,ltir_24],['LTIR_ALL','LTIR_24micron'],
                        line='ltir_mean',galaxy=galaxy['NAME'],
                        outdir=plot_dir)
