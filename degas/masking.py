from spectral_cube import SpectralCube
import astropy.units as u
from scipy.ndimage import map_coordinates, binary_dilation
import numpy as np
import copy
import astropy.wcs as wcs
from astropy.io import fits
import os
import glob

mad_to_std_fac = 1.482602218505602

def mad_zero_centered(data, mask=None):

    # Taken from PHANGS pipeline. I believe what this does is use the
    # data less than zero to estimate the mad under the assumption
    # that it's symmetric. Then it uses that to cut out what is
    # significant and estimate the noise based on the remaining
    # data. Basically it's a way to estimate the noise in a signal
    # filled FOV. It assumes no absorption (or not significant
    # absorption). For cases where there is little emission should be
    # very similar or the same as a regular mad.

    import scipy.stats as ss
    from astropy.stats import mad_std
    
    # ignore FOV masks
    where_data_valid = np.logical_and(~np.isnan(data), np.isfinite(data))

    if mask is None:
        sig_false = ss.norm.isf(0.5 / data.size) # greater than 0.5 probability that you get a false detection?
        data_lt_zero = np.less(data, 0,
                               where=where_data_valid,
                               out=np.full(where_data_valid.shape,
                                           False, dtype=bool))
        mad1 = mad_to_std_fac * np.abs(np.median(data[data_lt_zero])) #not quite sure I understand the math here.
     
        data_lt_mad1 = np.less(data, sig_false * mad1,
                               where=where_data_valid,
                               out=np.full(where_data_valid.shape,
                                           False, dtype=bool))
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data_lt_mad1])))
    else:
        nData = mask.sum()
        if nData == 0:
            logger.info('No data in mask. Returning NaN, which will now break things')
            return(np.nan)
        sig_false = ss.norm.isf(0.5 / nData)
        data_lt_zero = np.logical_and(np.less(data, 0, 
                                      where=where_data_valid,
                                      out=np.full(where_data_valid.shape,
                                                  False, dtype=bool)), mask)
        mad1 = mad_to_std_fac * np.abs(np.median(data[data_lt_zero]))
        data_lt_mad1 = np.logical_and(np.less(data, (sig_false * mad1),
                                      where=where_data_valid,
                                      out=np.full(where_data_valid.shape,
                                                  False, dtype=bool)), mask)
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data_lt_mad1])))
    return(mad2)

def noise_cube(data, mask=None, box=None, spec_box=None,
               nThresh=30, iterations=1,
               bandpass_smooth_window=None,
               bandpass_smooth_order=3):
    """

    # Taken from PHANGS pipeline

    Makes an empirical estimate of the noise in a cube assuming that it 
    is normally distributed.
    
    Parameters:
    -----------
    
    data : np.array
        Array of data (floats)
    
    Keywords:
    ---------
    
    mask : np.bool
        Boolean array with False indicating where data can be 
        used in the noise estimate. (i.e., True is Signal)
    
    box : int
        Spatial size of the box over which noise is calculated (correlation 
        scale).  Default: no box
    
    spec_box : int
        Spectral size of the box overwhich the noise is calculated.  Default:
        no box
    
    nThresh : int
        Minimum number of data to be used in a noise estimate.
    
    iterations : int
        Number of times to iterate the noise solution to force Gaussian 
        statistics.  Default: no iterations.
    
    bandpass_smooth_window : int
        Number of channels used in bandpass smoothing kernel.  Defaults to 
        nChan / 4 where nChan number of channels.  Set to zero to suppress 
        smoothing. Uses Savitzky-Golay smoothing
        
    bandpass_smooth_order : int
        Polynomial order used in smoothing kernel.  Defaults to 3.
    
    """

    from astropy.convolution import convolve, Gaussian2DKernel
    from scipy.signal import savgol_coeffs
    
    noisemask = np.isfinite(data)
    if mask is not None:
        noisemask[mask] = False
    step = 1
    boxr = step // 2 # floor division
    if box is not None:
        step = np.floor(box/2.5).astype(np.int)
        boxr = int(box // 2)
    if spec_box is not None:
        spec_step = np.floor(spec_box / 2).astype(np.int)
        boxv = int(spec_box // 2)
    else:
        boxv = 0

    noise_cube_out = np.ones_like(data)

    if bandpass_smooth_window is None:
        bandpass_smooth_window = 2 * (data.shape[0] // 8) + 1
        
    for i in np.arange(iterations):
        noise_map = np.zeros(data.shape[1:]) + np.nan
        noise_spec = np.zeros(data.shape[0]) + np.nan
        xx = np.arange(data.shape[2])
        yy = np.arange(data.shape[1])
        zz = np.arange(data.shape[0])
        for x in xx[boxr::step]:
            for y in yy[boxr::step]:
                spec = data[:, (y-boxr):(y+boxr+1),
                            (x-boxr):(x+boxr+1)]
                spec_mask = noisemask[:, (y-boxr):(y+boxr+1),
                                      (x-boxr):(x+boxr+1)]
                if np.sum(spec_mask) > nThresh:
                    noise_map[y, x] = mad_zero_centered(spec, mask=spec_mask)

        if boxr > 0:
            data_footprint = np.any(np.isfinite(data), axis=0)
            kernel = Gaussian2DKernel(box / np.sqrt(8 * np.log(2)))
            wt_map = np.isfinite(noise_map).astype(np.float)
            noise_map[np.isnan(noise_map)] = 0.0
            noise_map = convolve(noise_map, kernel)
            wt_map = convolve(wt_map, kernel)
            noise_map /= wt_map
            noise_map[~data_footprint] = np.nan

        for z in zz:
            lowz = np.clip(z - boxv, 0, data.shape[0])
            hiz = np.clip(z + boxv + 1, 0, data.shape[0])
            plane = data[lowz:hiz, :, :] / noise_map[np.newaxis, :, :]
            plane_mask = noisemask[lowz:hiz, :, :]
            noise_spec[z] = mad_zero_centered(plane, mask=plane_mask)
        # Smooth spectral shape
        if bandpass_smooth_window > 0:
            kernel = savgol_coeffs(int(bandpass_smooth_window),
                                   int(bandpass_smooth_order))
            noise_spec = convolve(noise_spec, kernel, boundary='extend')
        noise_spec /= np.nanmedian(noise_spec)
        noise_cube = np.ones_like(data)
        noise_cube *= (noise_map[np.newaxis, :]
                       * noise_spec[:, np.newaxis, np.newaxis])
        if iterations == 1:
            return(noise_cube)
        else:
            data = data / noise_cube
            noise_cube_out *= noise_cube
    return(noise_cube_out)




def cubemask(infile,outfile,
             peakCut=5.0,
             lowCut=3.0, 
             minBeamFrac = 1.0,
             minNchan = 1.0,
             skipChan = None,
             threeD = False,
             noise3D = False,
             outDir='./mask'):
    """ Creates a mask for a cube
    
    Parameters
    ----------

    infile : str
          input file

    output file : str
          output file

    lowCut : float
          Minimum signal-to-noise to expand the initial mask down to

    peakCut : float
          Minimum signal-to-noise for initial mask

    """

    from astrodendro import Dendrogram
    from astrodendro import pruning as p
    from matplotlib import pyplot as plt
    from scipy import ndimage

    # read in the data.
    cube = SpectralCube.read(infile)

    # set up cube and initialize mask
    cubedata = cube.unmasked_data[:,:,:].value
    mask = np.zeros(cube.shape, dtype=bool)

    # get info about cube
    nchan = cube.spectral_axis.size #number of channels

    # calculate noise
    if noise3D:
        # inputs based on PHANGS defaults.
        noise = noise_cube(cubedata, box=40, spec_box=5,bandpass_smooth_order=2,iterations=3)
        #SpectralCube(noise, cube.wcs, beam=cube.beam).write(os.path.join(outDir,outfile).replace('.fits','_noise.fits'),format='fits',overwrite=True)
    else: 
        noise = cube.mad_std().value # one value for whole cube. already scaled.

    sncube = cubedata / noise  # operate on the S/N cube to make everything easier

    #SpectralCube(sncube, cube.wcs, beam=cube.beam).write(os.path.join(outDir,outfile).replace('.fits','_sncube.fits'),format='fits',overwrite=True)

    if threeD:
        # compute dendrogram with a min threshold (input), 1sigma contrast, 
        # 1 beamsize as lower limit and a peak value lower limit(input)
        # all distinct regions will need to have a peak value > peakCut*sigma, 
        # or else it will be merged
        d = Dendrogram.compute(sncube,
                               min_value = lowCut,
                               min_delta = 1.0,
                               min_npix = minBeamFrac * cube.pixels_per_beam * minNchan,
                               is_independent = p.min_peak(peakCut))
        
        for t in d.trunk:
            mask = mask | t.get_mask()
    
    else:
        for chan in np.arange(nchan):
            d = Dendrogram.compute(sncube[chan,:,:],
                                   min_value = lowCut,
                                   min_delta = 1.0,
                                   min_npix = minBeamFrac * cube.pixels_per_beam,
                                   is_independent = p.min_peak(peakCut))
        
            for t in d.trunk:
                mask[chan,:,:] = mask[chan,:,:] | t.get_mask()  
    
    # blank channels you want to skip
    if skipChan:
        for chan in skipChan:
            mask[chan,:,:] = mask[chan,:,:]*0.0

    # fill in holes in individual channels
    for chan in np.arange(nchan):
        mask[chan,:,:] = ndimage.binary_fill_holes(mask[chan,:,:]) # fill holes in the mask

    maskhead = cube.header
    maskhead['BUNIT'] = ''
    cubemask = SpectralCube(data=mask.astype('short'),wcs=cube.wcs,header=maskhead)
    cubemask.write(os.path.join(outDir,outfile),format='fits',overwrite=True)


def deduplicate_keywords(hdr):
    for kw in ['CDELT4','CRPIX4','CRVAL4','CTYPE4']:
        if kw in hdr:
            hdr.remove(kw, ignore_missing=True)

    klist = [k for k in hdr.keys()]
    duplicates = [k for k in klist if klist.count(k) > 1]
    for kw in duplicates:
        val = hdr[kw]
        hdr.remove(kw, ignore_missing=True)
        hdr[kw] = val
    return(hdr)


def buildmasks(filename, nChan=2000, width=2e9, outdir=None,
               dilate=0, emissionfile=None):
    """Builds masks for use in DEGAS imaging pipeline. 

    Parameters
    ----------

    filename : str
        FITS filename of spectral cube mask. The file should be a
        binary mask with True / 1 indicating emission and False / 0
        otherwise.  This assumes the cube has a spectral axis in
        velocity and that the cube or has the metadata required to
        convert to velocity.  Note there is no checking of the
        spectral frame (LSRK, LSRD, BARY) and the conversion assumes
        radio Doppler convention.
  
    nChan : int
        Number of channels in output mask.  This should be larger than
        the number of channels in the DEGAS bandpass (1024)

    width : float
        Spectral width in Hz of the resulting mask.  This should be
        larger than the GBT bandwidth used (usually 1.5 GHz for DEGAS)

    outdir : str
        Directory for output masks to be stored in

    dilate : int
        Spectrally dilate the mask by this number of channels

    """


    if outdir is None:
       outdir = os.environ['DEGASDIR'] + 'masks/'
       
    if not os.access(outdir,os.W_OK):
        try:
            os.mkdir(outdir)
            print('Made directory {0}'.format(outdir))
        except OSError:
            try:
                os.mkdir('/'.join((outdir.split('/'))[0:-1])) # there may be a safer what to do this with os.path.split
                os.mkdir(outdir)
                print('Made directory {0}'.format(outdir))
            except:
                warnings.warn('Unable to make output directory '+outdir)
                raise
        except:
            warnings.warn('Unable to make output directory '+outdir)
            raise


    c = 299792.458
    # Read in original cube, ensure in velocity space
    s = SpectralCube.read(filename)
    s = s.with_spectral_unit(u.km / u.s, velocity_convention='radio')
    vmid = s.spectral_axis[len(s.spectral_axis)//2].value
    if emissionfile:
        cube = SpectralCube.read(emissionfile)
        cube = cube.with_spectral_unit(u.km / u.s,
                                       velocity_convention='radio')
        s = cube.with_mask(s > 0 * s.unit)
        s = s.with_fill_value(0 * s.unit)
        dtype = np.float
        outtype = np.float
    else:
        dtype = np.bool
        outtype = np.uint8
    # HCN_HCO+
    # Build a mask with a spectral width of 2 GHz and the same spatial 
    # dimensions as the original mask
  
    s_hcn = s.with_spectral_unit(u.Hz, rest_value=88.631847 * u.GHz)
    s_hcop = s.with_spectral_unit(u.Hz, rest_value=89.188518 * u.GHz)

    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=outtype)
    hdr = s_hcn.wcs.to_header()
    hdr['CRPIX3'] = 1000
    hdr['CDELT3'] = width / nChan
    hdr['CRVAL3'] = (89.188518 + 88.631847) / 2 * 1e9 * (1 - vmid / c)
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = mask.shape[0]
    hdr['NAXIS2'] = mask.shape[1]
    hdr['NAXIS3'] = mask.shape[2]
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 8
    hdr['EXTEND'] = 'T'

    hdr = deduplicate_keywords(hdr)
    w = wcs.WCS(hdr)
    maskcube = SpectralCube(mask, w, header=hdr)
    for zz in range(nChan):
        nu = maskcube.spectral_axis[zz]
        _, _, zz_hcn = s_hcn.wcs.wcs_world2pix(hdr['CRVAL1'],
                                               hdr['CRVAL2'],
                                               nu, 0)
        zz_hcn = int(zz_hcn)
        _, _, zz_hcop = s_hcop.wcs.wcs_world2pix(hdr['CRVAL1'],
                                                 hdr['CRVAL2'],
                                                 nu, 0)
        zz_hcop = int(zz_hcop)
        if 0 <= zz_hcn < s_hcn.shape[0]:
            mask[zz, :, :] = np.array(s_hcn.filled_data[zz_hcn, :, :],
                                      dtype=dtype)
        if 0 <= zz_hcop < s_hcop.shape[0]:
            mask[zz, :, :] = np.array(s_hcop.filled_data[zz_hcop, :, :],
                                      dtype=dtype)
    if dilate > 0 and (dtype == np.bool):
        mask = binary_dilation(mask,
                               np.ones((int(2 * dilate + 1), 1, 1),
                                       dtype=dtype))
    maskcube = SpectralCube(mask.astype(outtype), w, header=hdr)
    galname = os.path.split(filename)[1].split('_')[0]
    
    maskcube.write(outdir + galname+'.hcn_hcop.mask.fits',
                   overwrite=True)

    # C18O/13CO
    # Build a mask with a spectral width of 2 GHz and the same spatial
    # dimensions as the original mask

    s_13co = s.with_spectral_unit(u.Hz, rest_value=110.20135 * u.GHz)
    s_c18o = s.with_spectral_unit(u.Hz, rest_value=109.78217 * u.GHz)

    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=outtype)
    hdr = s_13co.wcs.to_header()
    hdr['CRPIX3'] = 1000
    hdr['CDELT3'] = width / nChan
    hdr['CRVAL3'] = (110.20135 + 109.78217) / 2 * 1e9 * (1 - vmid / c)
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = mask.shape[0]
    hdr['NAXIS2'] = mask.shape[1]
    hdr['NAXIS3'] = mask.shape[2]
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 8
    hdr['EXTEND'] = 'T'
    w = wcs.WCS(hdr)
    hdr = deduplicate_keywords(hdr)
    
    maskcube = SpectralCube(mask, w, header=hdr)
    for zz in range(nChan):
        nu = maskcube.spectral_axis[zz]
        _, _, zz_13co = s_13co.wcs.wcs_world2pix(hdr['CRVAL1'],
                                               hdr['CRVAL2'],
                                               nu, 0)
        zz_13co = int(zz_13co)
        _, _, zz_c18o = s_c18o.wcs.wcs_world2pix(hdr['CRVAL1'],
                                                 hdr['CRVAL2'],
                                                 nu, 0)
        zz_c18o = int(zz_c18o)
        if 0 <= zz_13co < s_13co.shape[0]:
            mask[zz, :, :] = np.array(s_13co.filled_data[zz_13co, :, :],
                                      dtype=dtype)
        if 0 <= zz_c18o < s_c18o.shape[0]:
            mask[zz, :, :] = np.array(s_c18o.filled_data[zz_c18o, :, :],
                                      dtype=dtype)
    if dilate > 0 and dtype == np.bool:
        mask = binary_dilation(mask,
                               np.ones((int(2 * dilate + 1), 1, 1),
                                       dtype=dtype))

    maskcube = SpectralCube(mask.astype(outtype), w, header=hdr)
    maskcube.write(outdir + galname + '.13co_c18o.mask.fits',
                   overwrite=True)
    
    # 12CO
    # Build a mask with a spectral width of 2 GHz and the same spatial
    # dimensions as the original mask

    s_12co = s.with_spectral_unit(u.Hz, rest_value=115.271204 * u.GHz)


    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=outtype)
    hdr = s_12co.wcs.to_header()
    hdr['CRPIX3'] = 1000
    hdr['CDELT3'] = width / nChan
    hdr['CRVAL3'] = (115.271204) * 1e9
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = mask.shape[0]
    hdr['NAXIS2'] = mask.shape[1]
    hdr['NAXIS3'] = mask.shape[2]
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 8
    hdr['EXTEND'] = 'T'
    w = wcs.WCS(hdr)
    hdr = deduplicate_keywords(hdr)
    
    maskcube = SpectralCube(mask, w, header=hdr)
    for zz in range(nChan):
        nu = maskcube.spectral_axis[zz]
        _, _, zz_12co = s_12co.wcs.wcs_world2pix(hdr['CRVAL1'],
                                               hdr['CRVAL2'],
                                               nu, 0)
        zz_12co = int(zz_12co)
        if 0 <= zz_12co < s_12co.shape[0]:
            mask[zz, :, :] = np.array(s_12co.filled_data[zz_12co, :, :],
                                      dtype=dtype)
    if dilate > 0 and dtype == np.bool:
        mask = binary_dilation(mask,
                               np.ones((int(2 * dilate + 1), 1, 1),
                                       dtype=dtype))

    maskcube = SpectralCube(mask.astype(outtype), w, header=hdr)
    maskcube.write(outdir + galname+'.12co.mask.fits', overwrite=True)

def build_setup_masks(**kwargs):
    OutputRoot = os.environ["DEGASDIR"]
    fl = glob.glob(OutputRoot + 'masks/NGC????_mask.fits')
    fl += glob.glob(OutputRoot + 'masks/IC????_mask.fits')
    for thisfile in fl:
        buildmasks(thisfile, **kwargs)
