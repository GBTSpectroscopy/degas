import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from importlib import reload
from scipy.ndimage import gaussian_filter

#----------------------------------------------------------------------
#                           read in files
#----------------------------------------------------------------------
databaseDir = os.path.join(os.environ['ANALYSISDIR'],'database')
scriptDir = os.environ['SCRIPTDIR']
z0mgs = Table.read(os.path.join(databaseDir,'apjsab3925t4_mrt.txt'),format='cds')
degas = Table.read(os.path.join(scriptDir,'degas_base.fits'),format='fits')

#----------------------------------------------------------------------
# Plot z0mgs sample to show star-forming main sequence
#----------------------------------------------------------------------


idx = degas['DR1'] == 1
idxout = degas['DR1'] == 0

plt.plot(degas[idx]['LOGMSTAR'],
         degas[idx]['LOGSFR']-degas[idx]['LOGMSTAR'], marker='o',
         markerfacecolor='orange',label='DEGAS DR1',
         linestyle='None',markeredgecolor='black',
         markersize='8',zorder=10)


idx = (((degas['LOGSFR'] - degas['LOGMSTAR']) > -9.75) & 
       (degas['DR1'] == 1))

for i in np.arange(len(degas['NAME'][idx])):
    plt.annotate(degas['NAME'][idx][i], (degas['LOGMSTAR'][idx][i],(degas['LOGSFR'][idx]-degas['LOGMSTAR'][idx])[i]),zorder=3, textcoords='offset points',xytext=(7,-7))

idx = z0mgs['D'] < 50.0 # only plot things less than 50Mpc.

plt.plot(z0mgs['logM*'][idx],
         z0mgs['logSFR'][idx] - z0mgs['logM*'][idx], 
         marker='.',
         markerfacecolor='gray',
         alpha=0.5,
         linestyle='None',markeredgecolor='None',zorder=4,
         label='Galaxies D < 50Mpc\n(Leroy+ 2019)')


plt.xlim(9.5,12.25)
plt.ylim(-11,-8.75)

massRange = np.logspace(9.5,12.0,100)

logSSFR = -0.32 * np.log10(massRange/10**10) - 10.17 

plt.plot(np.log10(massRange),logSSFR,
         marker='None',
         linestyle='-',
         color='black',
         linewidth='2',
         label='SF Main Sequence\n(Leroy+ 2019)',zorder=5)

plt.plot(np.log10(massRange),logSSFR-0.36,
         marker='None',
         linestyle='--',
         color='black',
         linewidth='2',zorder=5)

plt.plot(np.log10(massRange),logSSFR+0.36,
         marker='None',
         linestyle='--',
         color='black',
         linewidth='2',zorder=5)


plt.legend(loc='upper right')

plt.xlabel(r'$\log$ M$_*$ [$M_\odot$]')
plt.ylabel(r'$\log$ SFR/M$_*$ [yr$^{-1}$]')

plt.savefig(os.path.join(os.environ['ANALYSISDIR'],'plots','sample_degas_dr1_only_v2.pdf'),pad_inches=0.2,bbox_inches='tight')

plt.plot(degas[idxout]['LOGMSTAR'],
         degas[idxout]['LOGSFR']-degas[idxout]['LOGMSTAR'],
         marker='o',linestyle='None',
         markeredgecolor='black',markerfacecolor='blue',alpha=0.5,
         label='DEGAS Future',zorder=3,markersize='8')

plt.legend(loc='upper right')

plt.savefig(os.path.join(os.environ['ANALYSISDIR'],'plots','sample_degas_all_v2.pdf'),pad_inches=0.2,bbox_inches='tight')

for i in np.arange(len(degas['NAME'])):
    plt.annotate(degas['NAME'][i], (degas['LOGMSTAR'][i],(degas['LOGSFR']-degas['LOGMSTAR'])[i]),zorder=3)


plt.savefig(os.path.join(os.environ['ANALYSISDIR'],'plots','sample_degas_all_withnames_v2.pdf'), pad_inches=0.2,bbox_inches='tight')

plt.clf()

plt.close()
