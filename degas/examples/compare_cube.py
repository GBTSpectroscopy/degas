from spectral_cube import SpectralCube
import glob
from astropy.convolution import Box1DKernel
import astropy.units as u
from corner import corner
import numpy as np
import matplotlib.pyplot as plt
empdir = '/mnt/space/erosolow/surveys/empire/'
degasdir = '/mnt/space/erosolow/surveys/DEGAS/'
maskdir = '/mnt/space/erosolow/surveys/DEGAS/masks/'

applymask = True

gals = ['ngc2903','ngc4321','ngc5055','ngc6946']
# gals = gals[-1:]
for g in gals:
    for species in ['HCN','HCOp','13CO','C18O']:
#        degas = SpectralCube.read(degasdir + g.upper() + '/images/{0}_{1}_rebase7_smooth1.3_hanning1.fits'.format(g.upper(), species))
        try:
            degas = SpectralCube.read(degasdir + '/IR6p0/{0}_{1}_rebase7_smooth1.3_hanning1.fits'.format(g.upper(), species))
        except:
            continue
        fl = glob.glob(empdir + 'empire_{0}_{1}_*.fits'.format(g, species.lower()))
        empire = SpectralCube.read(fl[0])
        dv_ratio = ((degas.spectral_axis[1]-degas.spectral_axis[0]) / (empire.spectral_axis[1] - empire.spectral_axis[0])).to(u.dimensionless_unscaled).value
        if dv_ratio > 1:
            kernel = Box1DKernel(dv_ratio)
            empire = empire.spectral_smooth(kernel)
            empire = empire.spectral_interpolate(degas.spectral_axis)
        degas = degas.convolve_to(empire.beam)
        empire = empire.reproject(degas.header)
        emp = empire.filled_data[:].value
        deg = degas.filled_data[:].value
        p999 = np.nanpercentile(deg, 99.9) * 2
        p001 = np.nanpercentile(deg, 0.1) * 1.5

        if applymask:
            mask = SpectralCube.read(maskdir
                                     + '{0}_12CO_mask.fits'.format(g.upper()))
            mask.allow_huge_operations=True
            mask = mask.spectral_interpolate(degas.spectral_axis)
            mask = mask.reproject(degas.header, order='nearest-neighbor')
            mask = mask.filled_data[:].value
            mask = mask > 0
            idx = np.isfinite(emp) * np.isfinite(deg) * (mask)
        else:
            idx = np.isfinite(emp) * np.isfinite(deg)
            
        topfrac = np.logical_or((emp > np.nanpercentile(emp, 99)),
                                (deg > np.nanpercentile(deg, 99)))
        medrat = np.nanmedian(deg[topfrac] / emp[topfrac])
        f = corner(np.c_[emp[idx], deg[idx]], bins=100)
        f.axes[2].set_xlim([p001, p999])
        f.axes[2].set_ylim([p001, p999])
        f.axes[0].set_xlim([p001, p999])
        f.axes[3].set_xlim([p001, p999])
        f.axes[2].set_xlabel('EMPIRE')
        f.axes[2].set_ylabel('DEGAS')
        f.axes[2].plot([p001,p999], [p001,p999], color='r',
                       linewidth=3, alpha=0.4)
        f.axes[2].plot([p001, p999], [p001 * medrat,
                                      p999 * medrat],
                       color='b',linewidth=3, linestyle='--', alpha=0.4)
        f.text(0.6, 0.76, 
               '{0} {1}'.format(g.upper(), species), transform=plt.gcf().transFigure)
        f.text(0.6, 0.68, 'Median DEGAS/EMPIRE: {0:4.2f}'.format(medrat))
        f.set_size_inches(6,6)
        f.savefig('DEGAS_vs_EMPIRE_{0}_{1}.pdf'.format(g, species))
        degas.write('/mnt/space/erosolow/degas_{0}_{1}.fits'.format(g,species), overwrite=True)
        empire.write('/mnt/space/erosolow/empire_{0}_{1}.fits'.format(g,species), overwrite=True)
        # import pdb; pdb.set_trace()
        
# degas
# empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
# ls
# empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
# run /home/erosolow/sanitize
# empire = SpectralCube.read('../../../empire/EMPIRE_ngc4321_hcn_33as.fits')
# empire
# import astropy.units as u
# degas_con = degas.convolve_to(33*u.arcsec)
# ?degas.convolve_to
# from radio_beam import Beam
# degas_con = degas.convolve_to(Beam(major = 33*u.arcsec, minor=33*u.arcsec, pa=0*u.deg))
# empire_reproj = empire.reproject(degas_con.header)
# empire_reproj
# degas_con
# plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:])
# import matplotlib.pyplot as plt
# plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:])
# plt.show()
# plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:],logstretch=True)
# ?plt.hexbin
# plt.hexbin(empire_reproj.filled_data[:],degas_con.filled_data[:],bins='log')
# plt.show()
# from corner import corner
# corner(empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel())
# ?corner
# corner(np.c_[empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel()])
# import numpy as np
# corner(np.c_[empire_reproj.filled_data[:].ravel(), degas_con.filled_data[:].ravel()])
# corner(np.c_[empire_reproj.filled_data[:].ravel().value, degas_con.filled_data[:].ravel().value])
# emp = =empire_reproj.filled_data[:].value
# emp  =empire_reproj.filled_data[:].value
# deg = degas_con.filled_data[:].value
# idx = np.isfinite(emp) * np.isfinite(deg)
# corner(np.c_[emp[idx], deg[idx]])
# plt.show()
# ?corner
# corner(np.c_[emp[idx], deg[idx]])
# plt.set_xrange([-0.006, 0.012])
# plt.set_xlim([-0.006, 0.012])
# f = corner(np.c_[emp[idx], deg[idx]])
# f.axes
# f.axes[2].set_xlim([-0.006, 0.02])
# f.axes[2].set_ylim([-0.006, 0.02])
# plt.show()
# f = corner(np.c_[emp[idx], deg[idx]/0.3],bins=100)
# f.axes[2].set_xlim([-0.006, 0.02])
# f.axes[2].set_ylim([-0.006, 0.02])
# f.axes[0].set_xlim([-0.006, 0.02])
# f.axes[3].set_xlim([-0.006, 0.02])
# f.axes[2].set_xlabel('EMPIRE')
# f.axes[2].set_ylabel('DEGAS')
# f.set_size_inches(6,6)
# f.savefig('DEGAS_vs_EMPIRE_NGC4321_HCN.pdf')
