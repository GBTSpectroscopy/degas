from spectral_cube import SpectralCube
degas = SpectralCube.read('NGC4321_HCN_rebase3_smooth1.3_hanning1.fits')
degas
empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
ls
from astropy.io import fits
empire_hdu = fits.open('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
run santitize
run /home/erosolow/sanitize
empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
run /home/erosolow/sanitize
empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
empire
import astropy.units as u
degas_con = degas.convolve_to(33*u.arcsec)
?degas.convolve_to
from radio_beam import Beam
degas_con = degas.convolve_to(Beam(major = 33*u.arcsec, minor=33*u.arcsec, pa=0*u.deg))
empire_reproj = empire.reproject(degas_con.header)
empire_reproj
degas_con
plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:])
import matplotlib.pyplot as plt
plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:])
plt.show()
plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:],logstretch=True)
?plt.hexbin
plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:],bins='log')
plt.show()
from corner import corner
corner(empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel())
?corner
corner(np.c_[empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel()])
import numpy as np
corner(np.c_[empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel()])
corner(np.c_[empire_reproj.filled_data[:].ravel().value, degas_con.filled_data[:].ravel().value])
emp = =empire_reproj.filled_data[:].value
emp  =empire_reproj.filled_data[:].value
deg = degas_con.filled_data[:].value
idx = np.isfinite(emp) * np.isfinite(deg)
corner(np.c_[emp[idx], deg[idx]])
plt.show()
?corner
corner(np.c_[emp[idx], deg[idx]])
plt.set_xrange([-0.006, 0.012])
plt.set_xlim([-0.006, 0.012])
f = corner(np.c_[emp[idx], deg[idx]])
f.axes
f.axes[2].set_xlim([-0.006, 0.02])
f.axes[2].set_ylim([-0.006, 0.02])
plt.show()
f = corner(np.c_[emp[idx], deg[idx]/0.3],bins=100)
f.axes[2].set_xlim([-0.006, 0.02])
f.axes[2].set_ylim([-0.006, 0.02])
f.axes[0].set_xlim([-0.006, 0.02])
f.axes[3].set_xlim([-0.006, 0.02])
f.axes[2].set_xlabel('EMPIRE')
f.axes[2].set_ylabel('DEGAS')
f.set_size_inches(6,6)
f.savefig('DEGAS_vs_EMPIRE_NGC4321_HCN.pdf')

