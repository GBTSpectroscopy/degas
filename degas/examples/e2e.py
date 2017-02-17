# For calibration
import degas

# for Gridding
import glob
import gbtpipe

# For plotting
from spectral_cube import SpectralCube
import aplpy
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl


# Assuming we have a directory containing all the files.
degas.calscans('TGBT15A_901_34.raw.vegas', start=82,
               stop=105, 
               refscans=[80],
               outdir='testdir')

# Find all the files that we want to use in our map.
flist = glob.glob('testdir/*fits')

# Trim 100 channels from the end of each spectrum 
# (which is 1024 channels long, in total)
edgetrim = 100
# Use a 64 channel window for baselining on each end of the spectrum.
basebuff = 64
# Note that we also use a few channels in the middle.

# Grid them dataz.
gbtpipe.griddata(filelist=flist,
                 startChannel=edgetrim,
                 endChannel=1024-edgetrim,
                 baselineRegion = [slice(edgetrim,
                                         edgetrim+basebuff,1),
                                   slice(448,576,1),
                                   slice(1024-edgetrim-basebuff,
                                         1024-basebuff,1)],
                 flagRMS=True,
                 flagRipple=True)

# Now, we have a calibrated cube.  Let's make pictures.

mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'FreeSerif'

hcn = SpectralCube.read('NGC1068_HCN.fits')
hcn = hcn[:,35:75,35:75]
hcn = hcn.with_spectral_unit(u.km/u.s,
                             velocity_convention='radio',
                             rest_value=89.18852 * u.GHz)
hcnslab = hcn.spectral_slab(800*u.km/u.s,1400*u.km/u.s)
# This writes out an HCN only cube with the corect axes.
hcnslab.write('ngc1068_hcn.fits', overwrite=True)

hcop = SpectralCube.read('NGC1068_HCN.fits')
hcop = hcop[:,35:75,35:75]
hcop = hcop.with_spectral_unit(u.km/u.s,
                               velocity_convention='radio',
                               rest_value=88.63394 * u.GHz)
hcopslab = hcop.spectral_slab(800*u.km/u.s,1400*u.km/u.s)
# This writes out an HCO+ only cube with the corect axes.
hcopslab.write('ngc1068_hcop.fits',overwrite=True)

# Let's make Moments!
hcop_mom0 = hcopslab.moment(0,0)
hcn_mom0 = hcnslab.moment(0,0)

# And plot moments!
fig = aplpy.FITSFigure(hcn_mom0.hdu,figsize=(4.5,4.0))
fig.show_colorscale(vmin=-1,cmap='inferno')
fig.show_contour(data=hcop_mom0.hdu, levels=[3,6,9,12,15],colors='white')
fig.add_colorbar()
fig.colorbar.set_axis_label_text('Integrated Intensity (K km/s)')
fig.add_beam()
fig.savefig('ngc1068_hcnhcop.png')
