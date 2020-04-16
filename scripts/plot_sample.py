import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from importlib import reload
from scipy.ndimage import gaussian_filter

#----------------------------------------------------------------------
#                           read in files
#----------------------------------------------------------------------
databaseDir = '/lustre/cv/users/akepley/degas/database'
z0mgs = Table.read(os.path.join(databaseDir,'apjsab3925t4_mrt.txt'),format='cds')
degas = Table.read(os.path.join(databaseDir,'degas_base.fits'),format='fits')

#----------------------------------------------------------------------
# Plot z0mgs sample to show star-forming main sequence
#----------------------------------------------------------------------

# Using histogram2d to get the point density
z0mgsHist, xedges, yedges = np.histogram2d(z0mgs['logM*'],
                                           z0mgs['logSFR']-z0mgs['logM*'],
                                           density=True,
                                           bins=(np.linspace(8,12,20),
                                                 np.linspace(-13,-7.5,20)))

# manipulating to get points and values
z0mgsHist=z0mgsHist.T
xcenters = (xedges[:-1] + xedges[1:])/2.0
ycenters = (yedges[:-1] + yedges[1:])/2.0

#z0mgsHistGauss = gaussian_filter(z0mgsHist,sigma=0.25)

# plotting contours
plt.contour(xcenters,ycenters,
             z0mgsHist,levels=np.linspace(0.1,0.50,10),zorder=1)

# THE CONTOURS HERE ARE KIND OF CHUNKY. THERE SHOULD BE A WAY I CAN SMOOTH -- THE APPROACH MIGHT BE TO USE SCIPY.INTERPOLATE.GRIDDATA OR scipy.ndimage.gaussian_filter


idx = degas['DR1'] == 1
idxout = degas['DR1'] == 0
plt.scatter(degas[idx]['LOGMSTAR'],degas[idx]['LOGSFR']-degas[idx]['LOGMSTAR'],edgecolor='black',facecolor='orange',label='DEGAS DR1',zorder=3)


plt.xlim(9,11.5)
plt.ylim(-12,-8.5)
plt.legend(loc='upper left')

plt.xlabel(r'$\log$ M$_*$ [$M_\odot$]')
plt.ylabel(r'$\log$ SFR/M$_*$ [yr$^{-1}$]')

plt.savefig('sample_degas_dr1_only.pdf')

plt.scatter(degas[idxout]['LOGMSTAR'],degas[idxout]['LOGSFR']-degas[idxout]['LOGMSTAR'],edgecolor='black',facecolor='blue',alpha=0.5,label='DEGAS Future',zorder=2)

plt.legend(loc='upper left')

plt.savefig('sample_degas_all.pdf')

for i in np.arange(len(degas['NAME'])):
    plt.annotate(degas['NAME'][i], (degas['LOGMSTAR'][i],(degas['LOGSFR']-degas['LOGMSTAR'])[i]))


plt.savefig('sample_degas_all_withnames.pdf')

plt.close()
