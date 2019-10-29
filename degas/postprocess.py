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
from astropy.io import fits
from radio_beam import Beam
from gbtpipe.Preprocess import buildMaskLookup
from gbtpipe.Baseline import robustBaseline
import warnings

def circletrim(cube, wtsFile, x0, y0, weightCut=0.2):
    wtvals = np.squeeze(fits.getdata(wtsFile))
    badmask = ~(wtvals > (weightCut * wtvals.max()))
    yy, xx = np.indices(wtvals.shape)
    dist = (yy - y0)**2 + (xx - x0)**2
    mindist = np.min(dist[badmask])
    mask = dist < mindist
    cubemask = (np.ones(cube.shape) * mask).astype(np.bool)
    cube = cube.with_mask(cubemask)
    cube = cube.minimal_subcube()
    return(cube)

def edgetrim(cube, wtsFile, weightCut=None):
    """
    This trims off the edges of the cubes based on where the weights
    are lower than required to get a good signal.
    """

    wtvals = np.squeeze(fits.getdata(wtsFile))

    if weightCut:
        mask = wtvals > weightCut
    else:
        mask = wtvals > wtvals.max() / 5 # 5 is ad hoc
        radius = 11
        elt = np.zeros((2*radius + 1, 2*radius+1))
        xx, yy = np.indices(elt.shape)
        # create a disk structuring element
        elt = ((xx - radius)**2 + (yy - radius)**2)<=radius
        mask = nd.binary_closing(mask, structure=elt)
    cubemask = (np.ones(cube.shape) * mask).astype(np.bool)
    cube = cube.with_mask(cubemask)
    cube = cube.minimal_subcube()
    return(cube)


def cleansplit(filename, galaxy=None,
               Vwindow = 650 * u.km / u.s,
               Vgalaxy = 300 * u.km / u.s,
               blorder=3, HanningLoops=0,
               maskfile=None,
               circleMask=True,
               edgeMask=False,
               weightCut=0.2,
               spectralSetup=None,
               spatialSmooth=1.0):
    """
    Takes a raw DEGAS cube and produces individual cubes for each
    spectral line.
    
    Paramters
    ---------
    filename : str
        The file to split.
    
    Keywords
    --------
    galaxy : Galaxy object
        Currently unused
    Vwindow : astropy.Quantity
        Width of the window in velocity units
    Vgalaxy : astropy.Quantity
        Line of sight velocity of the galaxy centre
    blorder : int
        Baseline order
    HanningLoops : int
        Number of times to smooth and resample the data
    edgeMask : bool
        Determine whether to apply an edgeMask
    weightCut : float
        Minimum weight value to include in the data
    spatialSmooth : float
        Factor to increase the (linear) beam size by in a convolution.
    spectralSetup : str
        String to determine how we set up the spectrum
        'hcn_hcop' -- split based on HCN/HCO+ setup
        '13co_c18o' -- split based on 13CO/C18O setup
        '12co' -- don't split; assume single line
    """

    Cube = SpectralCube.read(filename)
    CatalogFile = get_pkg_data_filename('./data/dense_survey.cat',
                                        package='degas')
    Catalog = Table.read(CatalogFile, format='ascii')

    # Find which galaxy in our catalog corresponds to the object we
    # are mapping
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
        galcoord = SkyCoord(MatchRow['RA'],
                            MatchRow['DEC'],
                            unit=(u.hourangle, u.deg))
        Galaxy = MatchRow['NAME'].data[0]
        print("Catalog Match with " + Galaxy)
        V0 = MatchRow['CATVEL'].data[0] * u.km / u.s

    # Check spectral setups.  Use the max frequencies present to
    # determine which spectral setup we used if not specifed.
    if spectralSetup is None:
        if (Cube.spectral_axis.max() > 105 * u.GHz and
            Cube.spectral_axis.max() < 113 * u.GHz):
            warnings.warn("assuming 13CO/C18O spectral setup")
            spectralSetup = '13CO_C18O'
            filestr = '13co_c18o'
        if (Cube.spectral_axis.max() > 82 * u.GHz and
            Cube.spectral_axis.max() < 90 * u.GHz):
            warnings.warn("assuming HCN/HCO+ spectral setup")
            spectralSetup = 'HCN_HCO+'
            filestr = 'hcn_hcop'
        if (Cube.spectral_axis.max() > 113 * u.GHz):
            warnings.warn("assuming 12CO spectral setup")
            spectralSetup = '12CO'
            filestr = '12co'

    if spectralSetup == '13CO_C18O': 
        CEighteenO = Cube.with_spectral_unit(u.km / u.s,
                                             velocity_convention='radio',
                                             rest_value=109.78217 * u.GHz)
        ThirteenCO = Cube.with_spectral_unit(u.km / u.s,
                                             velocity_convention='radio',
                                             rest_value=110.20135 * u.GHz)
        CubeList = (CEighteenO, ThirteenCO)
        LineList = ('C18O','13CO')

    elif spectralSetup == 'HCN_HCO+':
        HCN = Cube.with_spectral_unit(u.km / u.s,
                                      velocity_convention='radio',
                                      rest_value=88.631847 * u.GHz)
        HCOp = Cube.with_spectral_unit(u.km / u.s,
                                       velocity_convention='radio',
                                       rest_value=89.188518 * u.GHz)
        CubeList = (HCN, HCOp)
        LineList = ('HCN','HCOp')

    elif spectralSetup == '12CO':
        TwelveCO = Cube.with_spectral_unit(u.km / u.s,
                                           velocity_convention='radio',
                                           rest_value=115.27120180 * u.GHz)
        CubeList = (TwelveCO,)
        LineList = ('12CO',)

    for ThisCube, ThisLine in zip(CubeList, LineList):
        if circleMask:
            x0, y0, _ = ThisCube.wcs.wcs_world2pix(galcoord.ra,
                                                   galcoord.dec,
                                                   0, 0)
            ThisCube = circletrim(ThisCube,
                                  filename.replace('.fits','_wts.fits'),
                                  x0, y0,
                                  weightCut=weightCut)
        if edgeMask:
            ThisCube = edgetrim(ThisCube, 
                                filename.replace('.fits','_wts.fits'),
                                weightCut=weightCut)

        # Trim each cube to the specified velocity range
        ThisCube =ThisCube.spectral_slab(V0 - Vwindow, V0 + Vwindow)
        ThisCube.write(Galaxy + '_' + ThisLine + '.fits', overwrite=True)
        StartChan = ThisCube.closest_spectral_channel(V0 - Vgalaxy)
        EndChan = ThisCube.closest_spectral_channel(V0 + Vgalaxy)

        if maskfile is not None:
            maskLookup = buildMaskLookup(maskfile)
            shp = ThisCube.shape
            TmpCube = ThisCube.with_spectral_unit(u.Hz)
            spaxis = TmpCube.spectral_axis
            spaxis = spaxis.value
            data = ThisCube.filled_data[:].value
            for y in np.arange(shp[1]):
                for x in np.arange(shp[2]):
                    spectrum = data[:, y, x]
                    if np.any(np.isnan(spectrum)):
                        continue
                    coords = ThisCube.world[:, y, x]
                    mask = maskLookup(coords[2].value,
                                      coords[1].value,
                                      spaxis)
                    spectrum = robustBaseline(spectrum, blorder=blorder,
                                              baselineIndex=~mask)
                    data[:, y, x] = spectrum
            ThisCube = SpectralCube(data, Cube.wcs, header=Cube.header)
            ThisCube.write(Galaxy
                           + '_' + ThisLine
                           + '_rebase{0}.fits'.format(blorder), overwrite=True)
        else:
            gbtpipe.Baseline.rebaseline(Galaxy + '_' + ThisLine + '.fits',
                                        baselineRegion=[slice(0,StartChan,1),
                                                        slice(EndChan,
                                                              ThisCube.shape[0],1)],
                                        blorder=blorder)
        ThisCube = SpectralCube.read(Galaxy + '_' + ThisLine +
                                     '_rebase{0}'.format(blorder) + '.fits')
        # Smooth
        Kern = Kernel1D(array=np.array([0.5, 1.0, 0.5]))
        for i in range(HanningLoops):
            ThisCube.spectral_smooth(Kern)
            ThisCube = ThisCube[::2,:,:]
            
        # Spatial Smooth
        if spatialSmooth > 1.0:
            newBeam = Beam(major=ThisCube.beam.major * spatialSmooth,
                           minor=ThisCube.beam.minor * spatialSmooth)
            ThisCube.convolve_to(newBeam)
            smoothstr = '_smooth{0}'.format(spatialSmooth)
        else:
            smoothstr = ''

        # Final Writeout
        ThisCube.write(Galaxy + '_' + ThisLine +
                       '_rebase{0}'.format(blorder) + smoothstr +
                       '_hanning{0}.fits'.format(HanningLoops),
                       overwrite=True)

