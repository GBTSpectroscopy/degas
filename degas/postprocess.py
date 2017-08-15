import glob
import gbtpipe
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.convolution import Kernel1D
import scipy.ndimage as nd

def edgetrim(cube, wtsFile=None, weightCut=None):
    wts = SpectralCube.read(wtsFile)
    wtvals = np.squeeze(wts.filled_data[:])
    if weightCut:
        mask = wtvals > weightCut
    else:
        mask = wtvals > wtvals.max() / 5 # 5 is ad hoc
        radius = 3
        elt = np.zeros((2*radius + 1, 2*radius+1))
        xx, yy = np.indices(elt.shape)
        # create a disk structuring element
        elt = ((xx - radius)**2 + (yy - radius)**2)<=radius
        mask = nd.binary_closing(mask, structure=elt)
    cubemask = np.ones(cube.shape) * mask
    cube = cube.with_mask(cubemask)
    cube = cube.minimal_subcube()
    return


def cleansplit(filename, galaxy=None,
               Vwindow = 650 * u.km / u.s,
               Vgalaxy = 300 * u.km / u.s,
               blorder=3, HanningLoops=2,
               edgeMask=True,
               weightCut=None):
    Cube = SpectralCube.read(filename)
    CatalogFile = get_pkg_data_filename('./data/dense_survey.cat',
                                        package='degas')
                                    
    Catalog = Table.read(CatalogFile, format='ascii')
    if galaxy is None:
        RABound, DecBound = Cube.world_extrema
        match = np.zeros_like(Catalog, dtype=np.bool)
        for index, row in enumerate(Catalog):
            galcoord = SkyCoord(row['RA'],
                                row['DEC'],
                                unit=(u.hourangle, u.deg))
            if (galcoord.ra < RABound[1] and
                galcoord.ra > RABound[0] and
                galcoord.dec < DecBound[1] and
                galcoord.dec > DecBound[0]):
                match[index] = True
        MatchRow = Catalog[match]
        Galaxy = MatchRow['NAME'].data[0]
        print "Catalog Match with " + Galaxy
    V0 = MatchRow['CATVEL'].data[0] * u.km / u.s

    # Check spectral setups.  If > 100 GHz, assume it's the 13CO/C180
    if Cube.spectral_axis.max() > 100 * u.GHz:
        print "Assuming 13CO/C18O spectral setup..."
        CEighteenO = Cube.with_spectral_unit(u.km / u.s,
                                             velocity_convention='radio',
                                             rest_value=109.78217 * u.GHz)
        ThirteenCO = Cube.with_spectral_unit(u.km / u.s,
                                      velocity_convention='radio',
                                      rest_value=110.20135 * u.GHz)
        CubeList = (CEighteenO, ThirteenCO)
        LineList = ('C18O','13CO')
    else:
        print "Assuming HCN/HCO+ spectral setup..."
        HCN = Cube.with_spectral_unit(u.km / u.s,
                                      velocity_convention='radio',
                                      rest_value=88.631847 * u.GHz)
        HCOp = Cube.with_spectral_unit(u.km / u.s,
                                      velocity_convention='radio',
                                        rest_value=89.188518 * u.GHz)
        CubeList = (HCN, HCOp)
        LineList = ('HCN','HCOp')

    for ThisCube, ThisLine in zip(CubeList, LineList):
        ThisCube =ThisCube.spectral_slab(V0 - Vwindow, V0 + Vwindow)
        ThisCube.write(Galaxy + '_' + ThisLine + '.fits', overwrite=True)
        StartChan = ThisCube.closest_spectral_channel(V0 - Vgalaxy)
        EndChan = ThisCube.closest_spectral_channel(V0 + Vgalaxy)
        
        # Rebaseline
        gbtpipe.Baseline.rebaseline(Galaxy + '_' + ThisLine + '.fits',
                                    baselineRegion=[slice(0,StartChan,1),
                                                    slice(EndChan,
                                                          ThisCube.shape[0],1)],
                                    blorder=3)
        ThisCube = SpectralCube.read(Galaxy + '_' + ThisLine +
                                     '_rebase{0}'.format(blorder) + '.fits')
        # Smooth
        Kern = Kernel1D(array=np.array([0.5, 1.0, 0.5]))
        for i in range(HanningLoops):
            ThisCube.spectral_smooth(Kern)
            ThisCube = ThisCube[::2,:,:]
        
        # Trim edges.
        if edgeMask:
            edgetrim(ThisCube, 
                     wtsfile = filename.replace('.fits','_wts.fits'),
                     weightCut=weightCut)

        ThisCube.write(Galaxy + '_' + ThisLine +
                       '_rebase{0}'.format(blorder) +
                       '_smooth{0}.fits'.format(HanningLoops),
                       overwrite=True)

