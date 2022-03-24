##this script uses phangs-alma routines to derive moment 0 and 1 maps and error maps for DEGAS cubes##
## following EMPIRE, we will use CO as mask for generating moment products##
## what we have: cleaned degas cube smoothed to 15", mask made from CO smoothed to 15"
#what we want: mom0, mom1, and errormaps
#what we need: noise cube -- needed for estimating error maps, but perhaps can be removed after each run??

#modified from PHANGS code for DEGAS
import numpy as np
import scipy.ndimage as nd
import scipy.ndimage.morphology as morph
import scipy.stats as ss
from scipy.signal import savgol_coeffs
import astropy.wcs as wcs
import astropy.units as u
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.stats import mad_std
from spectral_cube import SpectralCube, Projection
import glob
from matplotlib import pyplot as plt
import inspect

import os
import sys
import ipdb

#based on PHANGS-example_make_products.py
#save the henerated noise cubes
def make_line_products(galaxy, line, inDir, maskDir, outDir, noise_kwargs):
    #logging.basicConfig(level=logging.DEBUG, 
    #                    handlers = [logging.FileHandler(os.path.join(outDir,'product.log'), mode='w'), logging.StreamHandler(sys.stdout)])
##########
#    make moment maps and associated error maps for one galaxy in DEGAS

#    galaxy: str, galaxy name, should be the same as the beginning of the file

#    inDir: input directory relative to the base directory, where cleaned cubes are stored

#    maskDir: mask directory, where regrided 3D CO masks are stored

#    outDir: output directory, where generated products are stored

#   do_map: whether to estimate noise at each pixel (False) or over a box (True)

#   do_spec: whether to estimate noise assuming same for all channels or over several channels (True, spec_box = 1 for now)

##########
    #check line input
    if line not  in ['HCN', 'HCOp', '13CO', 'C18O']:
        logger.error("Must be HCN, HCOp, 13CO or C18O!")

    # -------noise estimation from DEGAS cube----------#
    # grab DEGAS cube for the line and galaxy ---- might need to change the file names later to be less specific?
    this_base_infile = glob.glob(os.path.join(inDir,galaxy+'_'+line+'_rebase7_smooth1.3_hanning1*.fits'))[0]
    co_mask_file = glob.glob(os.path.join(maskDir,galaxy+'_12CO_mask_regrid.fits'))[0]
    rms_file = os.path.join(outDir, \
                            os.path.basename(this_base_infile).replace('.fits','_rms.fits'))
    print("Noise estimation for native res for ", this_base_infile)
    basecube= SpectralCube.read(this_base_infile)
    
    rmscube = recipe_degas_noise(incube=basecube, 
                                 mask=co_mask_file, 
                                 outfile=rms_file, 
                                 noise_kwargs=noise_kwargs,
                                 overwrite=True) 
   
    mom0_file = os.path.join(outDir, \
                             os.path.basename(this_base_infile).replace('.fits','_mom0.fits')) #need to change this later to NAME_LINE_RES_mom0.fits 
    emom0_file = os.path.join(outDir, os.path.basename(this_base_infile).replace('.fits','_emom0.fits')) #need to change this later to NAME_LINE_RES_emom0.fits
    mom0,emom0=write_moment0(
        basecube, rms=rmscube, channel_correlation=None, #not implemented yet in phangs-pipeline
        outfile=mom0_file, 
        errorfile=emom0_file,
        overwrite=True, unit=None,
        include_limits=True,
        line_width=30 * u.km / u.s,
        return_products=True)
    print ('mom0: ',mom0_file, '; error map:', emom0_file)

    mom1_file = os.path.join(outDir, os.path.basename(this_base_infile).replace('.fits','_mom1.fits')) #need to change this later to NAME_LINE_RES_mom1.fits 
    emom1_file = os.path.join(outDir, os.path.basename(this_base_infile).replace('.fits','_emom1.fits')) #need to change this later to NAME_LINE_RES_emom1.fits
    write_moment1(
        basecube, mom0, rms=rmscube, channel_correlation=None,
        outfile=mom1_file, 
        errorfile=emom1_file,
        overwrite=True, unit=None,
        return_products=False)
    print ('mom1: ',mom1_file, '; error map:', emom1_file)
    
    print ('DONE. Check ', outDir, 'for products.')
    os.chdir(outDir)
##########################################
########## Estimate noise cube #########
#####PHANGS: scNoiseRoutines.py###########
##########################################
mad_to_std_fac = 1.482602218505602

#grab noise estimation func from phangs
def mad_zero_centered(data, mask=None):
    """
    Estimates the noise in a data set using the median absolute
    deviation. Assumes a normal distribution and that the data are
    centered on zero. Excludes any masked regions.
    
    Parameters:
    -----------
    
    data : np.array

        Array of data (floats)
    
    Keywords:
    ---------
    
    mask : np.bool

        Boolean array with True indicating where data can be used in
        the noise estimate. (i.e., True is noise). If none is supplied
        a trivial mask is estimated.
        YSnote: This is opposite from the mask used in the main function to generate noise map!!!

    """

    # TBD: Add input error checking
    
    #DEGAS: input is spectral_cube, make it array only so it won't complain about too many indices
    #data = data.unmasked_data[:,:,:].value
    #Select finite data and exclude not-a-number
    where_data_valid = np.logical_and(~np.isnan(data), np.isfinite(data))

    # Make a trivial mask (to simplify the code) if none has been supplied
    if mask is None:
        mask = np.isfinite(data)
    else:
        #YS - changed from PHANGS
        ##flipping mask signs : True --> False for signals, False --> True for noise
        newmask = np.full(mask.shape,np.nan)
        newmask[mask==True]=False
        newmask[mask==False]=True
    
    # Check that we have enough data to estimate the noise
    #YS - changed from PHANGS
    nData = np.nansum(newmask)
    if nData == 0:
        logger.info('No data in mask. Returning NaN.')
        return(np.nan)

    # Estimate the significance threshold likely associated with
    # false positives based on the size of the data.

    sig_false = ss.norm.isf(0.5 / nData)

    # Make a first estimate of the noise based on the negatives.

    data_lt_zero = np.logical_and(
        np.less(data, 0, 
                where=where_data_valid,
                out=np.full(where_data_valid.shape,
                            False, dtype=bool)), newmask)
    
    mad1 = mad_to_std_fac * np.abs(np.median(data[data_lt_zero]))

    # Make a second estimate now including the positives less than
    # the false-positive threshold.

    data_lt_mad1 = np.logical_and(
        np.less(data, (sig_false * mad1),
                where=where_data_valid,
                out=np.full(where_data_valid.shape,
                            False, dtype=bool)), newmask)
    mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data_lt_mad1])))

    # Return this second estimate

    return(mad2)

def noise_cube(data, mask=None, 
               nThresh=10, iterations=1,
               do_map=True, do_spec=True, #DEGAS: TODO: try turning them off one at a time to see what happens (both off will give a single value)
               box=None, spec_box=None,  
               bandpass_smooth_window=0,
               bandpass_smooth_order=3,
               oversample_boundary=False):

    """

    Makes an empirical estimate of the noise in a cube assuming that
    it is normally distributed about zero. Treats the spatial and
    spectral dimensions as separable.
    
    Parameters:
    -----------
    
    data : np.array
        Array of data (floats)
    
    Keywords:
    ---------
    
    mask : np.bool

        Boolean array with False indicating where data can be 
        used in the noise estimate. (i.e., True is signal). 
        DEGAS: using masks generated from CO data from PHANGS.
    
    do_map : np.bool
    
        Estimate spatial variations in the noise. Default is True. If
        set to False, all locations in a plane have the same noise
        estimate.
        DEGAS: using Default = True.

    do_spec : np.bool
    
        Estimate spectral variations in the noise. Default is True. If
        set to False, all channels in a spectrum have the same noise
        estimate.
        DEGAS: using Default=True.

    box : int

        Spatial size of the box over which noise is calculated in
        pixels.  Default: no box, every pixel gets its own noise
        estimte.
        DEGAS: using Defaul now, may change later. 
    
    spec_box : int

        Spectral size of the box overwhich the noise is calculated.
        Default: no box, each channel gets its own noise estimate.
        DEGAS: using Defaul now, may change later.

    nThresh : int
        Minimum number of data to be used in an individual noise estimate.
        DEGAS: using 10 now , may change later.
    
    iterations : int
        Number of times to iterate the noise solution to force Gaussian 
        statistics.  Default: no iterations.
        DEGAS: using Defaul now, may change later.
    
    bandpass_smooth_window : int
        Number of channels used in bandpass smoothing kernel.  Defaults to 
        nChan / 4 where nChan number of channels.  Set to zero to suppress 
        smoothing. Uses Savitzky-Golay smoothing
        DEGAS: using Zero now to avoid smoothing.
        
    bandpass_smooth_order : int
        Polynomial order used in smoothing kernel.  Defaults to 3.
        DEGAS: currently not smoothing, so doesn't matter.
    
    """

    # TBD: add error checking
    
    # Create a mask that identifies voxels to be fitting the noise

    noisemask = np.isfinite(data)
    if mask is not None:
        #noisemask[mask] = False
        noisemask = mask.unmasked_data[:,:,:].value

    # Default the spatial step size to individual pixels

    step = 1
    halfbox = step // 2

    # If the user has supplied a spatial box size, recast this into a
    # step size that critically samples the box and a halfbox size
    # used for convenience.
    
    if box is not None:
        step = np.floor(box/2.5).astype(np.int)
        halfbox = int(box // 2)

    # Include all pixels adjacent to the spatial
    # boundary of the data as set by NaNs
    boundary = np.all(np.isnan(data), axis=0)

    if oversample_boundary:
        struct = nd.generate_binary_structure(2, 1)
        struct = nd.iterate_structure(struct, halfbox)
        rind = np.logical_xor(nd.binary_dilation(boundary, struct),
                              boundary)
        extray, extrax = np.where(rind)
    else:
        extray, extrax = None, None
        
    # If the user has supplied a spectral box size, use this to
    # calculate a spectral step size.

    if spec_box is not None:
        spec_step = np.floor(spec_box / 2).astype(np.int)
        boxv = int(spec_box // 2)
    else:
        boxv = 0

    # Default the bandpass smoothing window

    if bandpass_smooth_window is None:
        bandpass_smooth_window = 2 * (data.shape[0] // 8) + 1

    # Initialize output to be used in the case of iterative
    # estimation.

    noise_cube_out = np.ones_like(data)
      
    # Iterate
  
    for ii in np.arange(iterations):

        if not do_map:

            # If spatial variations are turned off then estimate a
            # single value and fill the noise map with this value.

            noise_value = mad_zero_centered(data, mask=noisemask)            
            noise_map = np.zeros(data.shape[1:]) + noise_value

        else:

            # Initialize map to be full of not-a-numbers
            noise_map = np.zeros(data.shape[1:]) + np.nan

            # Make a noise map

            xx = np.arange(data.shape[2])
            yy = np.arange(data.shape[1])

            # Sample starting at halfbox and stepping by step. In the
            # individual pixel limit this just samples every spectrum.

            xsamps = xx[halfbox::step]
            ysamps = yy[halfbox::step]
            xsampsf = (xsamps[np.newaxis,:]
                      * (np.ones_like(ysamps))[:,np.newaxis]).flatten()
            ysampsf = (ysamps[:,np.newaxis] * np.ones_like(xsamps)).flatten()

            for x, y in zip(xsampsf, ysampsf):
                # Extract a minicube and associated mask from the cube
                
                # extract each spectra from each box or pixel (if no box specified)
                minicube = data[:, (y-halfbox):(y+halfbox+1),
                               (x-halfbox):(x+halfbox+1)]
                minicube_mask = noisemask[:, (y-halfbox):(y+halfbox+1),
                                             (x-halfbox):(x+halfbox+1)]

                # If we have enough data, fit a noise value for this entry

                if np.sum(minicube_mask) > nThresh:
                    noise_map[y, x] = mad_zero_centered(minicube,
                                                        mask=minicube_mask)
            
            if extrax is not None and extray is not None:
                for x, y in zip(extrax, extray):

                    minicube = data[:, (y-halfbox):(y+halfbox+1),
                                    (x-halfbox):(x+halfbox+1)]
                    minicube_mask = noisemask[:, (y-halfbox):(y+halfbox+1),
                                              (x-halfbox):(x+halfbox+1)]

                    if np.sum(minicube_mask) > nThresh:
                        noise_map[y, x] = mad_zero_centered(minicube,
                                                            mask=minicube_mask)
            #ipdb.set_trace()
            noise_map[boundary] = np.nan

            # If we are using a box size greater than an individual pixel
            # interpolate to fill in the noise map.

            if halfbox > 0:

                # Note the location of data, this is the location
                # where we want to fill in noise values.
                data_footprint = np.any(np.isfinite(data), axis=0)

                # Generate a smoothing kernel based on the box size.
                kernel = Gaussian2DKernel(box / np.sqrt(8 * np.log(2)))

                # Make a weight map to be used in the convolution, in
                # this weight map locations with measured values have
                # unity weight. This without measured values have zero
                # weight.

                # wt_map = np.isfinite(noise_map).astype(np.float)
                # wt_map[boundary] = np.nan
                # Take an average weighted by the kernel at each
                # location.
                # noise_map[np.isnan(noise_map)] = 0.0
                # y, x = np.where(np.isfinite(noise_map))
                # import scipy.interpolate as interp
                # func = interp.interp2d(x, y, noise_map[y, x], kind='cubic')

                noise_map = convolve(noise_map, kernel, boundary='extend')

                # yy, xx = np.indices(noise_map.shape)
                # noise_map_beta = interp.griddata((y, x), noise_map[y,x],
                #                                  (yy, xx), method='cubic')
                # noise_map_beta = func(yy, xx)
                # noise_map_beta[boundary] = np.nan
                # noise_map = noise_map_beta
                # wt_map = convolve(wt_map, kernel, boundary='extend')

                # noise_map /= wt_map

                # Set the noise map to not-a-number outside the data
                # footprint.

                noise_map[~data_footprint] = np.nan
        # Initialize spectrum

        noise_spec = np.zeros(data.shape[0]) + np.nan

        if not do_spec:

            # If spectral variations are turned off then assume that
            # the noise_map describes all channels of the cube.
            #YSnote: basically cannot turn this off because everything will be multiplied by nan in the end
            #YS note: now the code it's update to spit out cube with the same noise_map across all channels
            
            #pass
            noise_spec = np.ones(data.shape[0]) #make everything = 1 so all channels will have the same nosie map
            
        else:

            # Loop over channels

            zz = np.arange(data.shape[0])
            for z in zz:

                # Idententify the range of channels to be considered
                # in this estimate.

                lowz = np.clip(z - boxv, 0, data.shape[0])
                hiz = np.clip(z + boxv + 1, 0, data.shape[0])
                
                # Extract a slab from the cube and normalize it by the
                # noise map. Now any measured noise variations are
                # relative to those in the noise map.

                slab = data[lowz:hiz, :, :] / noise_map[np.newaxis, :, :]
                slab_mask = noisemask[lowz:hiz, :, :]
                noise_spec[z] = mad_zero_centered(slab, mask=slab_mask)
                
            # Smooth the spectral variations in the noise.

            if bandpass_smooth_window > 0:

                # Initialize a Savitzky-Golay filter then run it over
                # the noise spectrum.

                kernel = savgol_coeffs(int(bandpass_smooth_window),
                                       int(bandpass_smooth_order))

                baddata = np.isnan(noise_spec)
                noise_spec = convolve(noise_spec, kernel,
                                      nan_treatment='interpolate',
                                      boundary='extend')
                noise_spec[baddata] = np.nan

                # Make sure that the noise spectrum is normalized by
                # setting the median to one.

                noise_spec /= np.nanmedian(noise_spec)

        # Combine the spatial and spectral variations into a
        # three-dimensional noise estimate.

        noise_cube = np.ones_like(data)
        #if noise_spec is nan, then everything is nan <----do_spec=False
        noise_cube *= (noise_map[np.newaxis, :]
                       * noise_spec[:, np.newaxis, np.newaxis])

        if iterations == 1:
            return(noise_cube)
        else:

            # If iterating, normalize the data by the current noise
            # estimate and scale the current noise cube by the new
            # estimate.

            data = data / noise_cube
            noise_cube_out *= noise_cube
        

    # If iterating return the iterated noise cube.

    return(noise_cube_out)

def recipe_degas_noise(
    incube=None,
    outfile=None,
    mask=None,
    noise_kwargs=None,
    return_spectral_cube=True, # DEGAS: return rms cube in SpectralCube format for moment generating
        overwrite=True): #DEGAS: overwrite on for testing
    """

    Wrap noise_cube with a set of preferred parameters for the
    PHANGS-ALMA CO work.
    
    Parameters:
    -----------
    
    cube : np.array

        Array of data (floats)
    
    Keywords:
    ---------
    
    mask : np.bool

        Boolean array with False indicating where data can be used in
        the noise estimate. (i.e., True is signal). 
        YSnote: this is opposite to mask definition in mad_zero_centered!!!!

    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking and work out inputs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


    if type(incube) is SpectralCube:
        cube = incube
    elif type(incube) == str:
        cube = SpectralCube.read(incube)
    else:
        print ("Input must be a SpectralCube object or a filename.")

    # Initialize an empty kwargs dictionary
    if noise_kwargs is None:
        noise_kwargs = {}

    # If no box is specified, default to one about two beams across
    if 'box' not in noise_kwargs:
        pixels_per_beam = cube.pixels_per_beam
        box = np.ceil(pixels_per_beam**0.5)
        noise_kwargs['box'] = box

    # Default to an odd bandpass smothing window 
    if 'bandpass_smooth_window' not in noise_kwargs:
        #spectral_smooth = np.ceil(cube.shape[0] / 5) // 2 * 2 + 1
        spectral_smooth = 0  #no smoothing for DEGAS for now
        noise_kwargs['bandpass_smooth_window'] = spectral_smooth

    if 'spec_box' not in noise_kwargs:
        noise_kwargs['spec_box'] = 1 #DEGAS: using default for now, maybe can use 3 as well
        #noise_kwargs['spec_box'] = None 

    if 'iterations' not in noise_kwargs:
        noise_kwargs['iterations'] = 1 #DEGAS: no iterating for now

    # Require a valid cube input as a mask
    if mask is not None:
        if type(mask) is SpectralCube:
            noise_kwargs['mask'] = mask
        elif type(mask) == type("hello"):
            noise_kwargs['mask'] = SpectralCube.read(mask)
        else:
            logger.error("Mask must be a SpectralCube object or a filename or None.")

    # Fill in the mask if it hasn't already been filled in.
    if 'mask' not in noise_kwargs:

        # Check if a non-trivial signal mask is attached to the cube
        if (np.sum(cube.mask.include())
            < np.sum(np.isfinite(cube.filled_data[:].value))):
            noise_kwargs['mask'] = cube.mask.include()
        else:
            noise_kwargs['mask'] = None

            
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Run the noise estimate
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    data = cube.filled_data[:].value
    badmask = np.isnan(data)
    badmask = nd.binary_dilation(badmask,
                                 structure=nd.generate_binary_structure(3, 2))
    data[badmask] = np.nan
    rms = noise_cube(data, 
                     **noise_kwargs)

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Write or return as requested
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # In this case can avoid a recast
    if not return_spectral_cube and (outfile is None):
        return(rms)

    # Recast from numpy array to spectral cube
    header = cube.header
    header['DATAMIN'] = np.nanmin(rms)
    header['DATAMAX'] = np.nanmax(rms)
    #header['COMMENT'] = 'Produced with PHANGS-ALMA pipeline version ' + version
    #if tableversion:
        #header['COMMENT'] = 'Galaxy properties from PHANGS sample table version ' + tableversion
    rms = SpectralCube(rms, wcs=cube.wcs, header=header,
                       meta={'BUNIT':cube.header['BUNIT']})
    
    # Write to disk, if desired
    if outfile is not None:
        rms.write(outfile, overwrite=overwrite)
        
    if return_spectral_cube:
        return(rms)
    else:
        return(rms.filled_data[:].value)

##########################################
########## Create Moment Maps #########
#####PHANGS: scMoments.py###########
##########################################

def _nicestr(quantity):
    if quantity.value == int(quantity.value):
        return(str(int(quantity.value))+' '+str(quantity.unit))
    else:
        return(str(quantity))
'''
def _func_and_kwargs_for_moment(moment_tag=None):
    """
    Return function name and defalt kwargs for a moment tag.
    """

    func = None
    kwargs = None
    if moment_tag is None:
        return(func,kwargs)

    if moment_tag == 'mom0':
        func = write_moment0
        kwargs ={'unit': u.K * u.km / u.s}
    elif moment_tag == 'mom1':
        func = write_moment1
        kwargs = {'unit': u.km / u.s}

    return(func, kwargs)

def moment_tag_known(moment_tag=None):
    """
    Test whether the programs know about a moment tag.
    """
    func, kwargs = _func_and_kwargs_for_moment(moment_tag)
    if func is None:
        return(False)
    return(True)


def moment_generator(
        cubein, mask=None, noise=None,
        moment=None, momkwargs=None,
        outfile=None, errorfile=None,
        channel_correlation=None,
        context=None, assignkunits=False):

    """
    Generate one moment map from input cube, noise, and masks.
    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Set up the call
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Get the relevant function and keyword arguments for this moment
    func, kwargs = _func_and_kwargs_for_moment(moment)
    if func is None:
        logging.error("Moment tag not recognized: "+str(moment))
        raise NotImplementedError
        return(None)

    # Add any user-supplied kwargs to the dictionary
    if momkwargs is not None:
        if type(momkwargs) != type({}):
            logging.error("Type of momkwargs should be dictionary.")
            raise NotImplementedError
        for this_kwarg in momkwargs:
            kwargs[this_kwarg] = momkwargs[this_kwarg]

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Read in the data
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Read in the cube (if needed)
    if type(cubein) is str:
        cube = SpectralCube.read(cubein)
    elif type(cubein) is SpectralCube:
        cube = cubein
    else:
        logging.error('Unrecognized input type for cubein')
        raise NotImplementedError

    cube.allow_huge_operations = True
        
    # Force Kelvin. We will be unit agnostic later.
    cube = cube.to(u.K)
    
    # Attach a mask if needed
    if mask is not None:
        if type(mask) is str:
            mask = SpectralCube.read(mask)
        elif type(mask) is SpectralCube:
            mask = mask
        else:
            logging.error('Unrecognized input type for mask')
            raise NotImplementedError

        # Ensure the mask is booleans and attach it to the cube. This
        # just assumes a match in astrometry. Could add reprojection
        # here or (better) build a masking routine to apply masks with
        # arbitrary astrometry.

        mask = np.array(mask.filled_data[:].value, dtype=np.bool)
        cube = cube.with_mask(mask, inherit_mask=False)

    # Read in the noise (if present)
    if noise is not None:        
        if type(noise) is str:
            noisecube = SpectralCube.read(noise)
        elif type(noise) is SpectralCube:
            noisecube = noise
        else:
            logging.error('Unrecognized input type for noise.')
            raise NotImplementedError

        noisecube.allow_huge_operations = True

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Call the moment generation
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Probably not needed anymore
    theseargs = (inspect.getfullargspec(func)).args

    if 'context' in theseargs:
        moment_map, error_map = func(
            cube, rms=noisecube,
            outfile=outfile, errorfile=errorfile,
            channel_correlation=channel_correlation,
            #context=context,
            **kwargs)
    else:
        moment_map, error_map = func(
            cube, rms=noisecube,
            outfile=outfile, errorfile=errorfile,
            channel_correlation=channel_correlation,
            **kwargs)
        
    return(moment_map, error_map)
'''
######################################
######### Helper functions  ##############
###### PHANGS: scDerivativeRoutines.py#######
# #####################################


# Force reduction in bit-depth to save space.
def writer(projection, filename, overwrite=True, dtype=np.float32):
    hdu = fits.PrimaryHDU(projection.hdu.data.astype(dtype),
                          header=projection._header)
    hdu.writeto(filename, overwrite=overwrite)

    
def update_metadata(projection, cube, error=False):

    keys = ['BMAJ', 'BMIN', 'BPA', 'JYTOK', 'VELREF',
            'TELESCOP', 'INSTRUME', 'ORIGIN', 'OBJECT',
            'TIMESYS','MJDREFI','MJDREFF','DATEREF']
    calling_name = inspect.getouterframes(inspect.currentframe())[1][3]
    btype_dict = {'write_moment0':'Moment0',
                  'write_moment1':'Moment1'}
    hdr = projection.header
    try:
        btype = btype_dict[calling_name]
        if error:
            btype = btype+' Error'
        hdr['BTYPE'] = btype
    except KeyError:
        btype = 'Product'

    for key in keys:
        try:
            hdr[key] = cube.header[key]
        except KeyError:
            pass

    # Check if the moment map is empty. If so, nanmax and nanmin
    # will not be finite and writing the header to disk will fail.
    if not np.isfinite(projection.filled_data[:].value).any():
        mx = 0.
        mn = 0.
    else:
        ind = np.isfinite(projection.filled_data[:])
        mx =  np.nanmax(projection.filled_data[ind])
        mn =  np.nanmin(projection.filled_data[ind])
        hdr['DATAMAX'] = mx.value
        hdr['DATAMIN'] = mn.value

    if 'moment_axis' in projection.meta.keys():
        idx = projection.meta['moment_axis']
    else:
        idx = 0 # Assume spectral moment

    collapse_name = cube.wcs.axis_type_names[::-1][idx]
    med_spaxis = np.abs(np.median(cube.spectral_axis[1:]
                                  - cube.spectral_axis[0:-1]))
    hdr['CHANWDTH'] = med_spaxis.value
    hdr.comments['CHANWDTH'] = 'Median channel width in {0}'.format(med_spaxis.unit)

    # Eliminate BEAM keyword in favor of BMIN, BMAJ, BPA
    try:
        del hdr['BEAM']
    except KeyError:
        pass

    #hdr['COMMENT'] = 'Produced with PHANGS-ALMA pipeline version ' + version
   # if tableversion:
       # hdr['COMMENT'] = 'Galaxy properties from PHANGS sample table version ' + tableversion
    hdr['COMMENT'] = (btype
                      + ' generated by collapsing cube over '
                      + collapse_name + ' axis.')
    
    comments = hdr['COMMENT']
    unique_comments = list(set(comments))
    try:
        del hdr['COMMENT']
    except KeyError:
        pass

    for comm in unique_comments:
        hdr.set('COMMENT', comm)

    projection._header = hdr
    return(projection)

def channel_width(cube):
    dv = np.median(np.abs(cube.spectral_axis[1:] 
                          - cube.spectral_axis[0:-1]))
    return(dv)

def build_covariance(spectrum=None,
                     rms=None,
                     channel_correlation=None,
                     index=None):
    """
    Build a covariance matrix from a channel_correlation vector
    
    Keywords:
    ---------
    
    spectrum : np.array
        One-dimensional array of spectrum values
        
    rms : np.array
        One-dimensional array containing the root-mean-squared error
        estimate for the values in the spectrum
    
    channel_correlation : np.array
        One-dimensional array containing the channel-to-channel 
        normalize correlation coefficients
    
    index : np.array
        Integer array indicating the spectral indices of the data 
        in the original cube
    """
    
    # Note that this assumes you are masking out values to make sure 
    # arrays stay the same shape as the input

    if index is None:
        index = np.arange(len(spectrum))
    if channel_correlation is None:
        return(np.diag(rms**2))
    if len(channel_correlation) == 1:
        return(np.diag(rms**2))
    distance = np.abs(index[:, np.newaxis] - index[np.newaxis, :])
    covar = rms[:, np.newaxis] * rms[np.newaxis, :]
    maxdist = len(channel_correlation)
    covar[distance >= maxdist] = 0
    covar[distance < maxdist] *= channel_correlation[distance[distance 
                                                              < maxdist]]
    return(covar)    

def calculate_channel_correlation(cube, length=1):
    """
    TBD - calculate the channel correlation.
    """
    raise NotImplementedError

# &%&%&&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Moment 0
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def write_moment0(
        cube, rms=None, channel_correlation=None,
        outfile=None, errorfile=None,
        overwrite=True, unit=None,
        include_limits=True,
        line_width=30 * u.km / u.s, #DEGAS default 30 km/s
        return_products=False):

    """Write out moment0 map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a moment0 map
    
    outfile : str
        File name of output file
        
    errorfile : str
        File name of map for the uncertainty
        
    rms : SpectralCube
        Root-mean-square estimate of the error.  This must have an estimate
        the noise level at all positions where there is signal, and only at 
        those positions.
        
    channel_correlation : np.array
        One-dimensional array containing the channel-to-channel 
        normalize correlation coefficients --> don't worry about this now
        
    overwrite : bool
        Set to True (the default) to overwrite existing maps if present. 
        
    unit : astropy.Unit
        Preferred unit for moment masks
    
    include_limits : bool
        If true: For masked lines of sight inside the data, set the
        moment0 value to 0 and the error to a value of 1sigma over
        line_width as specified below.

    line_width : astropy.Quantity
        Assumed line width for moment0 upper limit.  Default = 10 km/s  #set to 30km/s for DEGAS

    return_products : bool
        Return products calculated in the map

    """
    
    # Spectral cube collapse routine. Applies the masked, automatically
    mom0 = cube.moment0()
    valid = np.isfinite(mom0)
    if include_limits:
        observed = np.any(np.isfinite(cube._data), axis=0)
        mom0[np.logical_and(np.isnan(mom0), observed)] = 0.0

    # Handle the error.
    mom0err_proj = None

    if errorfile is not None and rms is None:
        logger.error("Moment 0 error requested but no RMS provided")

    if rms is not None:
        # Initialize the error map
        mom0err = np.empty(mom0.shape)
        mom0err.fill(np.nan)

        # Note the channel width
        dv = channel_width(cube)
        if include_limits:
            rmsmed = np.nanmedian(rms.filled_data[:].value, axis=0)
            mom0err[observed] = (rmsmed[observed]
                                 * (np.abs(line_width
                                           / dv).to(u.dimensionless_unscaled).value)**0.5)

        # Make a masked version of the noise cube
        rms = rms.with_mask(cube._mask.include(), inherit_mask=False)

        if channel_correlation is None:
            sumofsq = (rms * rms).sum(axis=0)
            mom0err[valid] = np.sqrt(sumofsq[valid])
        else:
            # Iterates over the cube one ray at a time
            yy, xx = np.where(valid)
            for y, x in zip(yy, xx):
                slc = (slice(None), slice(y,y+1,None), slice(x,x+1,None))
                mask = np.squeeze(cube.mask.include(view=slc))
                index = np.where(mask)[0]

                # One dimensional versions of the data and noise
                rms_spec = rms.flattened(slc).value
                spec = cube.flattened(slc).value

                # Build a covariance matrix given the channel correlation
                covar = build_covariance(
                    spectrum=spec, rms=rms_spec,
                    channel_correlation=channel_correlation,
                    index=index)

                # Collapse the covariance matrix into an integrated moment map
                mom0err[y, x] = (np.sum(covar))**0.5

        # Multiply by the channel width and assign correct units
        mom0err = u.Quantity(mom0err * dv.value,
                             cube.unit * dv.unit, copy=False)

        # Convert units if request
        if unit is not None:
            mom0err = mom0err.to(unit)
        
        # Convert from an array into a spectral-cube projection that
        # shares metadata with the moment map
        mom0err_proj = Projection(
            mom0err, wcs=mom0.wcs, header=mom0.header, meta=mom0.meta)

        #DEGAS - mask mom0 and emom0 based on SNR>=3 -- Work in Progress
        #right now the edge of the internal mask remains, need to figure out a way to remove it
        mom0[mom0.value <= 3*mom0err_proj.value]=np.nan
        mom0err_proj[np.isnan(mom0.value)]=np.nan
        mom0[np.isnan(mom0err_proj)]=np.nan #trim pb edges to match errormap

        # Write to disk if requested
        if errorfile is not None:
            mom0err_proj = update_metadata(mom0err_proj, cube, error=True)
            writer(mom0err_proj, errorfile, overwrite=overwrite)
            # mom0err_proj.write(errorfile, overwrite=overwrite)
    else:
        mom0err_proj = None

    # Convert units if requested
    if unit is not None:
        mom0 = mom0.to(unit)

    # If requested, write to disk
    if outfile is not None:
        mom0 = update_metadata(mom0, cube)
        writer(mom0, outfile, overwrite=overwrite)
        # mom0.write(outfile, overwrite=overwrite)

    # If requested, return
    if return_products:
        return(mom0, mom0err_proj)

# &%&%&&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Moment 1
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def write_moment1(
    cube, mom0, rms=None, channel_correlation=None,
    outfile=None, errorfile=None,
    overwrite=True, unit=None,
    return_products=True):
    """
    Write out moment1 map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a moment1 map
    
    outfile : str
        File name of output file
        
    errorfile : str
        File name of map for the uncertainty
        
    rms : SpectralCube
        Root-mean-square estimate of the error.  This must have an estimate
        the noise level at all positions where there is signal, and only at 
        those positions.
        
    channel_correlation : np.array
        One-dimensional array containing the channel-to-channel 
        normalize correlation coefficients
        
    overwrite : bool
        Set to True (the default) to overwrite existing maps if present. 
        
    unit : astropy.Unit
        Preferred unit for moment masks
        
    return_products : bool
        Return products calculated in the map
    """

    mom1 = cube.moment1()
    mom1err_proj = None
    spaxis = cube.spectral_axis.value

    if errorfile is not None and rms is None:
        logger.error("Moment 1 error requested but no RMS provided")

    if rms is not None:
        mom1err = np.empty(mom1.shape)
        mom1err.fill(np.nan)
        # Ensure the same mask applied to both.
        rms = rms.with_mask(cube._mask.include(), inherit_mask=False)
        valid = np.isfinite(mom1)

        if channel_correlation is None:
            vval = spaxis
            sum_T = cube.sum(axis=0).value
            numer = np.nansum(rms.filled_data[:].value**2
                              * (vval[:, np.newaxis, np.newaxis]
                                 - mom1.value[np.newaxis, :, :])**2,
                              axis=0)
            mom1err = (numer / sum_T**2)**0.5
            mom1err[np.isnan(mom1.value)] = np.nan
        else:
            yy, xx = np.where(valid)
            for y, x in zip(yy, xx):
                slc = (slice(None), slice(y,y+1,None), slice(x,x+1,None))
                mask = np.squeeze(cube.mask.include(view=slc))
                index = np.where(mask)[0]
                rms_spec = rms.flattened(slc).value
                spec = cube.flattened(slc).value
                covar = build_covariance(spectrum=spec,
                                         rms=rms_spec,
                                         channel_correlation=channel_correlation,
                                         index=index)

                vval = spaxis[index]
                sum_T = np.sum(spec)
                sum_vT = np.sum(spec * vval)

                jacobian = vval / sum_T - sum_vT / sum_T**2
                mom1err[y, x] = np.dot(
                    np.dot(jacobian[np.newaxis, :], covar),
                    jacobian[:, np.newaxis])**0.5
        mom1err = u.Quantity(mom1err, cube.spectral_axis.unit, copy=False)
        if unit is not None:
            mom1err = mom1err.to(unit)
        mom1err_proj = Projection(mom1err,
                                  wcs=mom1.wcs,
                                  header=mom1.header,
                                  meta=mom1.meta)

        #DEGAS - mask mom0 and emom0 based on SNR>=3 -- Work in Progress
        #right now the edge of the internal mask remains, need to figure out a way to remove it
        mom1[np.isnan(mom0.value)]=np.nan
        mom1err_proj[np.isnan(mom0.value)]=np.nan
        mom1[np.isnan(mom1err_proj)]=np.nan #trim pb edges

        if errorfile is not None:
            mom1err_proj = update_metadata(mom1err_proj, cube, error=True)
            writer(mom1err_proj, errorfile, overwrite=overwrite)
            # mom1err_proj.write(errorfile, overwrite=overwrite)
    
    if unit is not None:
        mom1 = mom1.to(unit)
    if outfile is not None:
        mom1 = update_metadata(mom1, cube)
        writer(mom1, outfile, overwrite=overwrite)
        # mom1.write(outfile, overwrite=True)

    if return_products and mom1err_proj is not None:
        return(mom1, mom1err_proj)
    elif return_products and mom1err_proj is None:
        return(mom1)

