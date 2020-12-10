from astropy.table import Table
import matplotlib.pyplot as plt
import os

# read in data files
degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
stack = Table.read('/users/akepley/cv_lustre/degas/stack_test/test_mom1.fits')

dr1 = degas[degas['DR1'] == 1]

mymarkers = ['o','o','o','o','o','o','o','o','o','o','s','s','s','s','s','s','s']


for (galaxy,marker) in zip(dr1,mymarkers):
    idx = ((stack['galaxy'] == galaxy['NAME']) & (stack['bin_type'] == 'radius'))
    plt.loglog(stack[idx]['HCN_stack_sum'], stack[idx]['ltir_mean'],
               marker=marker,linestyle=':',
               label=galaxy['NAME'])


plt.xlabel(r'$L_{HCN}$ (K km s$^{-1}$ pc$^{2}$)') ### CHECK ON UNITS
plt.ylabel(r'$S_{TIR}$ (L$_{\odot}$)') #### CHECK ON UNITS
plt.grid(linestyle=":",color='gray')


plt.legend()
           
           
plt.show()    
