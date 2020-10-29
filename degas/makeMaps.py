#make S/N map using MAD from cube and mask
import numpy as np
import os
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube, Projection
def mapSN(galaxy, regridDir, outDir):
    cube=SpectralCube.read(os.path.join(regridDir, galaxy+'_12CO_regrid.fits'))
    cube=cube.with_spectral_unit(u.km / u.s)
    mask=SpectralCube.read(os.path.join(regridDir, galaxy+'_12CO_mask_regrid.fits'))
    mask=mask.with_spectral_unit(u.km / u.s)
    madstd=cube.mad_std(how='cube') #K #raw cube
    chanwidth=np.abs(cube.spectral_axis[0]-cube.spectral_axis[1]) #channel width is same for all channels, km/s
    masksum=mask.sum(axis=0) #map number of unmasked pixels 
    noise=np.sqrt(masksum)*(madstd*chanwidth) #moment0 error map, in K km/s
    #mask datacube
    masked_cube=cube.with_mask(mask==1.0*u.dimensionless_unscaled) 
    mom0 = masked_cube.moment(order=0)
    snmap=mom0/noise #should be unitless #mom0 is from masked cube
    snmap[snmap==np.inf]=np.nan #convert division by zero to nan
    snmap.write(os.path.join(outDir, galaxy.upper()+'_SNmap.fits'), overwrite=True)
    snmap.quicklook()
    plt.savefig(os.path.join(outDir, galaxy.upper()+'_SNmap.png'))
    plt.close()
    plt.clf()
    #get rid of parts of mom0 where S/N<3
    mom0cut=mom0.copy()
    sn=np.nan_to_num(snmap.value)
    mom0cut[sn<3.0]=np.nan #blank out low sn regions
    #use sigma-clipped mom0 as new 2D mask for the original cube to preserve noise in signal-free channel
    mom0mask=~np.isnan(mom0cut)
    masked_cube=cube.with_mask(mom0mask) #use for stacking later
    return mom0cut, masked_cube, mom0mask

#make stellarmass map
#TO DO: update once get new ancillary data from Sarah
def mapStellar(galaxy, mom0cut, regridDir, outDir):
    stellarhdu=fits.open(os.path.join(regridDir,galaxy+'_w1_stellarmass_regrid.fits'))[0]
    starmask=fits.getdata(os.path.join(regridDir,galaxy+'_w1_gauss15_stars_regrid.fits'))
    stellar=stellarhdu.data
    stellar[starmask==1.0]=np.nan #apply star mask
    stellar[np.isnan(mom0cut)]=np.nan #apply SN mask (SN >3)
    w=WCS(stellarhdu.header)
    stellarhdu.header['BUNIT']='Msun/pc^2' #not the correct unit!!
    stellarmap=Projection(stellar,header=stellarhdu.header,wcs=w) 
    stellarmap.quicklook()
    plt.savefig(os.path.join(outDir,galaxy+'_stellarmass.png'))
    plt.clf()
    plt.close()
    return stellarmap

#import sfr map from W4+FUV
def mapSFR(galaxy, mom0cut, regridDir, outDir):
    sfrhdu=fits.open(os.path.join(regridDir,galaxy+'_w4fuv_sfr_regrid.fits'))[0]
    sfr=sfrhdu.data
    sfr[np.isnan(mom0cut)]=np.nan
    w=WCS(sfrhdu.header)
    sfrhdu.header['BUNIT']='MSUN/YR/KPC^2'  #might change when getting new files from Sarah
    sfrmap=Projection(sfr,header=sfrhdu.header,wcs=w) 
    sfrmap.quicklook()
    plt.savefig(os.path.join(outDir,galaxy+'_sfr.png'))
    plt.clf()
    plt.close()
    return sfrmap
