import numpy as np
from astropy.io import fits
from astrodendro import Dendrogram
from astrodendro import pruning as p
from matplotlib import pyplot as plt
from astropy.stats import mad_std


def cubemask(filename,low_cut,peak_cut):
    #import datacube and info
    cube=fits.getdata(filename)
    cubehead=fits.getheader(filename)
    nchan=cubehead['NAXIS3'] #number of channels
    bmaj=cubehead['BMAJ'] # degree
    bmin=cubehead['BMIN'] # degree
    pxscale=np.abs(cubehead['CDELT1']) # degree per pixel
    beampx=(np.pi*bmaj*bmin/(4*np.log(2)))/pxscale**2 #gaussian beamsize in pixels
    madstd=mad_std(cube,ignore_nan=True) #rms estimated from MAD
    
    #get combined mask from all parent structures
    mask = np.zeros(cube.shape, dtype=bool)
    for i in range(nchan):
        #for each channel map
        slice=cube[i,:,:]
        #compute dendrogram with a min threshold (input), 1sigma contrast, 1 beamsize as lower limit and a peak value lower limit(input)
        #all distinct regions will need to have a peak value > 5sigma, or else it will be merged
        d=Dendrogram.compute(slice,min_value=low_cut*madstd,min_delta=1*madstd,min_npix=beampx,is_independent=p.min_peak(peak_cut*madstd))
        for t in d.trunk:
            #assemble channel masks
            mask[i,:,:] = mask[i,:,:] | t.get_mask()
            mask[i,:,:] = ndimage.binary_fill_holes(mask[i,:,:]) #fill holes in the mask
  
    
    #create and export mask cube and masked datacube
    cubemask = fits.PrimaryHDU(mask.astype('short'), cubehead)
    cubemask.writeto(filename.replace('.fits','_mask.fits'),overwrite=True)
    masked_cube=np.multiply(mask,cube)
    newcube = fits.PrimaryHDU(masked_cube, cubehead)
    newcube.writeto(filename.replace('.fits','_masked_data.fits'),overwrite=True)
