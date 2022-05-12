import glob
import os
import numpy as np

from astropy.io import fits
from astropy import units as u

from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy.convolution import Box1DKernel

from reproject import reproject_interp
from astropy import wcs

import re

import ipdb

def fixBIMA(bimafile):
    '''

    Purpose: fix up the bima cubes

    The smoothed BIMA data doesn't have the outer regions
    appropriately masked. This routine copies the mask from the
    original data to the smoothed data. 

    They also don't have all the header information correct.

    '''

    #get smoothed bima data, bima cube has 4 axes
    cube = fits.getdata(bimafile)[0,:,:,:]
    cubehead = fits.getheader(bimafile)

    #get the non-smoothed data, bima cubes have 4 axes
    oldfile = bimafile.replace('_gauss15.fits','_k.fits')
    oldcube = fits.getdata(oldfile)[0,:,:,:]
    oldmask = oldcube==0.0 #extract the original '0'valued pb mask
    cube[oldmask] = np.nan #where the '0' mask was, mask as nan

    #now data is 3d delete info on 4th axis ('stokes')
    cubehead.remove('CDELT4')
    cubehead.remove('CRPIX4')
    cubehead.remove('CRVAL4')
    cubehead.remove('CTYPE4')

    # convert BIMA HEADERS from m/s to km/s
    if 'CUNIT3' not in cubehead:
        cubehead['CDELT3'] = (cubehead['CDELT3']/1000.0)
        cubehead['CRVAL3'] = (cubehead['CRVAL3']/1000.0)
        cubehead['CUNIT3'] = 'km/s'

    # add CUNIT to RA and DEC axis for BIMA headers.
    if 'CUNIT2' not in cubehead:
        cubehead['CUNIT2'] = 'deg'

    if 'CUNIT1' not in cubehead:
        cubehead['CUNIT1'] = 'deg'

    # Set BUNIT to appropriate value. Original in Jy/beam, but I think
    # it's mis-labeled.
    cubehead['BUNIT'] = 'K'

    fixedcube = fits.PrimaryHDU(cube, cubehead)#fixed cube will have 3axes only
    fixedcube.writeto(bimafile.replace('.fits','_fixhdr.fits'),overwrite=True)

    # Now trying to change VOPT to VRADIO via SpectralCube
    cube2 = SpectralCube.read(bimafile.replace('.fits','_fixhdr.fits'))
   
    ## change to frequency and interpolate to equal frequency intervals
    cube2freq = cube2.with_spectral_unit(u.GHz,rest_value=115.27120180*u.GHz)
    freq_axis = cube2freq.spectral_axis
    chan_width = np.mean(freq_axis[0:-2] - freq_axis[1:-1])
    new_freq_axis = freq_axis[0] - chan_width * np.arange(0,len(freq_axis))
    cubenewfreq = cube2freq.spectral_interpolate(new_freq_axis)
    
    # convert to radio
    cube2radio = cubenewfreq.with_spectral_unit(u.km/u.s,rest_value=115.27120180*u.GHz,velocity_convention='radio')
    cube2radio.write(bimafile.replace('.fits','_fixed.fits'),overwrite=True)

def fixHeracles(fitsimage):
    '''

    fixes up image headers for missing values etc so that we can
    regrid appropriately.
 
    '''

    # Using the astropy.io.fits to modify headers 
    # since SpectralCube disapproves of modifying headers as far 
    # as I can tell.
    f = fits.open(fitsimage)

    # Convert VELOCITY to VRAD for Heracles headers.
    if f[0].header['CTYPE3'] == 'VELOCITY':
        f[0].header['CTYPE3'] = 'VRAD' 
        
    newimage = fitsimage.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)

    cube = SpectralCube.read(newimage)
    cube_kms = cube.with_spectral_unit(u.km / u.s)

    cube_kms.write(newimage.replace('.fits','_kms.fits'),overwrite=True)

def fixOVRO(fitsimage,beam=15.0):
    '''
    fix up OVRO data.
    '''
    
    # open image
    cube = SpectralCube.read(fitsimage)    
    
    # convert to K
    kcube = cube.to(u.K)
    
    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = kcube.convolve_to(newBeam)
    
    # write out
    smoothCube.write(fitsimage.replace('.fits','_gauss15_fixed.fits'),
                     overwrite=True)
    

def fixNRO(fitsimage):
    '''
    fix up NRO data
    '''
    
    cube = SpectralCube.read(fitsimage)
    # add beam
    beam_cube = cube.with_beam(Beam(15.0*u.arcsec))
    beam_cube.write(fitsimage.replace('.FITS','_fixed.fits'),overwrite=True)
    

def fixExtraHERA(fitsimage,beam=15.0, eff=1.58621):
    '''

    Fix the extra heracles data.

    This data needs an efficiency correction added. From the IRAM 30m
    docs, the correction at 230GHz is:
    Beff = 0.59
    Feff= 0.92
    Ta*/TMB = Feff/Beff = 0.92/0.59 = 1.559322

    1/1.559322 = 0.6413

    From the header on some of Adam's data the value is 
    Feff / Beff = 0.92/0.58 = 1.58621

    1/1.58621 = 0.630433

    Ta*/1.58621 = Tmb

    '''
    
    f = fits.open(fitsimage)
    f[0].header['CUNIT3'] = 'm/s'
    newimage = fitsimage.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)
    f.close()
    
    # open image
    cube = SpectralCube.read(newimage)    

    cube_kms = cube.with_spectral_unit(u.km / u.s)

    if re.search('ngc3631',newimage):
        # extra channels with data
        subCube = cube_kms[72:379,:,:]
    elif re.search('ngc4030',newimage):
        subCube = cube_kms[78:384,:,:]
    else:
        subCube = cube_kms[:,:,:]

    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = subCube.convolve_to(newBeam)
    
    # smooth velocity
    smoothFactor = 4.0
    spSmoothCube = smoothCube.spectral_smooth(Box1DKernel(smoothFactor))

    # Interpolate onto a new axis
    spec_axis = spSmoothCube.spectral_axis
    chan_width = spec_axis[1]-spec_axis[0] # channels are equally spaced in velocity
    new_axis = np.arange(spec_axis[0].value,spec_axis[-1].value,smoothFactor*chan_width.value) * u.km/u.s

    interpCube = spSmoothCube.spectral_interpolate(new_axis,
                                                   suppress_smooth_warning=True)

    if eff:
        interpCube = interpCube / eff
        interpCube.meta['EFF'] = eff

    # write out
    interpCube.write(newimage.replace('.fits','_10kms_gauss15.fits'),
                     overwrite=True)

def fixExtraHERAfromAdam(fitsimage,beam=15.0):
    '''
    Purpose: fix up extra HERA data from Adam
    '''

    # switch header from M/S to m/s to fix up wcs read errors with SpectralCube
    f = fits.open(fitsimage)
    f[0].header['CUNIT3'] = 'm/s'
    newimage = fitsimage.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)
    f.close()
    
    # open image
    cube = SpectralCube.read(newimage)    
    
    # switch to km/s
    cube_kms = cube.with_spectral_unit(u.km / u.s)

    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = cube_kms.convolve_to(newBeam)

    # write out
    smoothCube.write(newimage.replace('.fits','_10kms_gauss15.fits'),overwrite=True)

def fixPhangs(image,beam=15.0):
    '''Purpose: fix up Phangs data so we can do analysis for DEGAS.

    4/1/2021: I needed to modify, spectral_cube/base_class.py,
    _get_filled_data function to explicitly use memory mapping. Hack
    courtesy of Eric Koch. Evidently the cube was memory mapped, but
    not the mask.

    return self._mask._filled(data=data, wcs=self._wcs, fill=fill,
                                  view=view, wcs_tolerance=self._wcs_tolerance,
                                  #use_memmap=use_memmap
                                  use_memmap=False # hack to get phangs data to work
                                 )

    '''
    #os.environ['TMPDIR'] = '/lustre/cv/users/akepley/tmp'

    #cube = SpectralCube.read(image)
    cube = SpectralCube.read(image, use_memmap=False, mode='readonly')

    if re.search('ngc2903',image):
        subCube = cube[59:332,:,:]
    else:
        subCube = cube[:,:,:]
    
    cube_kms = subCube.with_spectral_unit(u.km / u.s)

    smoothFactor = 4.0 
    cube_kms.allow_huge_operations=True
    spSmoothCube = cube_kms.spectral_smooth(Box1DKernel(smoothFactor),
                                            use_memmap=False,verbose=1)

    # Interpolate onto a new axis
    spec_axis = spSmoothCube.spectral_axis
    chan_width = spec_axis[1]-spec_axis[0] # channels are equally spaced in velocity
    new_axis = np.arange(spec_axis[0].value,spec_axis[-1].value,smoothFactor*chan_width.value) * u.km/u.s

    interpCube = spSmoothCube.spectral_interpolate(new_axis,
                                                   suppress_smooth_warning=True)

    newBeam = Beam(beam*u.arcsec)

    #interpCube.allow_huge_operations=True
    smoothCube = interpCube.convolve_to(newBeam)

    smoothCube.write(image.replace('_7p5as.fits','_10kms_gauss15.fits'),overwrite=True)

def fixJialu(image,beam=15.0):
    '''
    Purpose: fix up Jialu's IC0342 cube so we can use it for analysis for DEGAS
    '''

    #12CO rest frequency
    rest_freq_12co = 115.27120180*u.GHz
    
    f = fits.open(image)

    f[0].header['BUNIT'] = 'K'
    f[0].header['RESTFRQ'] = float(rest_freq_12co.to(u.Hz).value)    

    f.writeto(image.replace('.fits','_fixed.fits'), overwrite=True)

    image =  image.replace('.fits','_fixed.fits')

    # open image
    cube = SpectralCube.read(image)   

    cube_kms = cube.with_spectral_unit(u.km / u.s) 

    ## This should set header
    #cube_ms = cube.with_spectral_unit(u.m / u.s, rest_value=rest_freq_12co)

    # chop off bad edge
    subcube = cube_kms.subcube(xlo=0,xhi=134,ylo=0,yhi=150)

    # add beam
    beamcube = subcube.with_beam(Beam(8.0*u.arcsec))

    # smooth
    newBeam = Beam(beam*u.arcsec)
    smoothCube = beamcube.convolve_to(newBeam)

    # smooth spectrally.
    smoothFactor = 3.0
    spSmoothCube = smoothCube.spectral_smooth(Box1DKernel(smoothFactor))
    spec_axis = spSmoothCube.spectral_axis
    chan_width = spec_axis[1]-spec_axis[0] # channels are equally spaced in velocity
    new_axis = np.arange(spec_axis[0].value,spec_axis[-1].value,smoothFactor*chan_width.value) * u.km/u.s

    interpCube = spSmoothCube.spectral_interpolate(new_axis,
                                                   suppress_smooth_warning=False)

    # write out
    interpCube.write(image.replace('8arcsec_fixed.fits','10kms_gauss15.fits'),
                     overwrite=True)


def fixEMPIRE(fitsimage):
    '''
    Purpose: fix up EMPIRE data headers to work with SpectralCube
    '''

    # switch header from M/S to m/s to fix up wcs read errors with SpectralCube
    f = fits.open(fitsimage)
    f[0].header['CUNIT3'] = 'm/s'
    newimage = fitsimage.replace('.fits','_fixed.fits')
    f.writeto(newimage,overwrite=True)
    f.close()
    
    # open image
    cube = SpectralCube.read(newimage)    
    
    # switch to km/s
    cube_kms = cube.with_spectral_unit(u.km / u.s)

    # write out
    cube_kms.write(newimage.replace('.fits','_kms.fits'),overwrite=True)

#----------------------------------------------------------------------

def calc_stellar_mass(iracFile, MtoLFile, outFile):

    '''
    Calculate the stellar mass surface density for a given IRAC image
    and MtoL image.
    '''
    
    # open IRAC data file
    irac = fits.open(iracFile)
    iracData = irac[0].data
    iracHdr = irac[0].header
        
    # open MtoL data file
    MtoL = fits.open(MtoLFile)
    MtoL_regrid = reproject_interp(MtoL,output_projection=iracHdr,
                                   hdu_in=0,
                                   order='nearest-neighbor',
                                   return_footprint=False)
    irac.close()
    MtoL.close()

    # Equation 8 from Leroy+ 2019 (z0mgs paper)
    mstarirac1 = 3.5e2 * (MtoL_regrid / 0.5) * (iracData)
        
    # create output file
    outHdr = iracHdr
    outHdr['BUNIT']= 'MSUN/PC^2'
    outfits = fits.PrimaryHDU(data=mstarirac1,header=outHdr)
    outfits.writeto(outFile, overwrite=True)        

# ----------------------------------------------------------------------

def calc_LTIR(outFile, b24=None, b70=None, b100=None, b160=None, w3=None, w4=None):
    '''
    Calculate LTIR based on the calibration in Galametz+ 2013.
    '''

    import math

    # useful constants
    Lsun = 3.826e33 # erg/s
    c = 2.99792458e10 #cm/s

    # Units for input data are MJy/sr and units for Galametz+2013
    # constants are W/kpc^2, so I need to do a conversation. 

    # unit conversion factor to go from MJy/sr to W/kpc^2. Needs to be
    # multiplied by nu to get final answer.
    unitfactor = 1e6 * 1e-26 * 4.0 * math.pi * (3.0861e16*1000.0)**2
   
    # now put together the various combinations of data.
    if (b24 and b70 and b100 and b160):
        print('calculating LTIR using 24, 70, 100, and 160 micron data')

        # coefficients from Table 3 in Galametz+ 2013
        c24 = 2.051
        c70 = 0.521
        c100 = 0.294
        c160 = 0.934
        
        # read in each image
        f24 = fits.open(b24)
        f70 = fits.open(b70)
        f100 = fits.open(b100)
        f160 = fits.open(b160)
        
        # regrid data to same grid (use 24micron as base) and 
        # convert from MJy/sr to W/kpc^2
        basehdr = f24[0].header

        data24freq = c / (24.0*1e-4)
        data24 = f24[0].data * unitfactor * data24freq

        data70freq = c / (70.0*1e-4)
        data70_regrid = reproject_interp(f70,output_projection=basehdr,
                                         hdu_in=0,
                                         order='nearest-neighbor',
                                         return_footprint=False)
        data70_regrid = data70_regrid * unitfactor * data70freq

        data100freq = c /(100.0*1e-4)
        data100_regrid = reproject_interp(f100,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data100_regrid = data100_regrid * unitfactor * data100freq
        

        data160freq = c/(160.0*1e-4)
        data160_regrid = reproject_interp(f160,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data160_regrid = data160_regrid * unitfactor * data160freq
        
        
        f24.close()
        f70.close()
        f100.close()
        f160.close()

        # calculate LTIR using coefficients from Table 3 in Galametz+ 213
        STIR_Wkpc2 = c24*data24 + c70*data70_regrid+c100*data100_regrid + \
                     c160*data160_regrid

        # convert from W/kpc^2 to Lsun/pc^2 (which is what we typically plot)
        # W/kpc^2 * 1e7 = erg/s/kpc^2
        # erg/s/kpc^2 / Lsun / (1000.0)**2 = Lsun/pc^2
        STIR_Lsunpc2 = (STIR_Wkpc2 * 1e7) / Lsun / (1000.0**2) 

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    elif (b70 and b100 and b160):
        print('calculating LTIR using 70, 100, and 160 micron data')

        # coefficients from Table 3 in Galametz+ 2013
        c70 = 0.789
        c100 = 0.387
        c160 = 0.960
        
        # read in each image
        f70 = fits.open(b70)
        f100 = fits.open(b100)
        f160 = fits.open(b160)
        
        # regrid data to same grid (use 24micron as base) and 
        # convert from MJy/sr to W/kpc^2
        basehdr = f70[0].header

        data70freq = c / (70.0*1e-4)
        data70 = f70[0].data * unitfactor * data70freq
        
        data100freq = c /(100.0*1e-4)
        data100_regrid = reproject_interp(f100,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data100_regrid = data100_regrid * unitfactor * data100freq
        

        data160freq = c/(160.0*1e-4)
        data160_regrid = reproject_interp(f160,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data160_regrid = data160_regrid * unitfactor * data160freq
        
        
        f70.close()
        f100.close()
        f160.close()

        # calculate LTIR using coefficients from Table 3 in Galametz+ 213
        STIR_Wkpc2 = c70*data70 + c100*data100_regrid + c160*data160_regrid

        # convert from W/kpc^2 to Lsun/pc^2 (which is what we typically plot)
        # W/kpc^2 * 1e7 = erg/s/kpc^2
        # erg/s/kpc^2 / Lsun / (1000.0)**2 = Lsun/pc^2
        STIR_Lsunpc2 = (STIR_Wkpc2 * 1e7) / Lsun / (1000.0**2) 

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)
        
    elif (b24 and b100 and b160): ## COLOR TERM?
        print('calculating LTIR using 24, 100, and 160 micron data')

        # coefficients from Table 3 in Galametz+ 2013
        c24 = 2.708
        c100 = 0.734
        c160 = 0.739
        
        # read in each image
        f24 = fits.open(b24)
        f100 = fits.open(b100)
        f160 = fits.open(b160)
        
        # regrid data to same grid (use 24micron as base) and 
        # convert from MJy/sr to W/kpc^2
        basehdr = f24[0].header

        data24freq = c / (24.0*1e-4)
        data24 = f24[0].data * unitfactor * data24freq

        data100freq = c /(100.0*1e-4)
        data100_regrid = reproject_interp(f100,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data100_regrid = data100_regrid * unitfactor * data100freq
        

        data160freq = c/(160.0*1e-4)
        data160_regrid = reproject_interp(f160,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data160_regrid = data160_regrid * unitfactor * data160freq
        
        
        f24.close()        
        f100.close()
        f160.close()

        # calculate LTIR using coefficients from Table 3 in Galametz+ 213
        STIR_Wkpc2 = c24*data24 + c100*data100_regrid + \
                     c160*data160_regrid

        # convert from W/kpc^2 to Lsun/pc^2 (which is what we typically plot)
        # W/kpc^2 * 1e7 = erg/s/kpc^2
        # erg/s/kpc^2 / Lsun / (1000.0)**2 = Lsun/pc^2
        STIR_Lsunpc2 = (STIR_Wkpc2 * 1e7) / Lsun / (1000.0**2) 

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    elif (b100 and b160): ## COLOR TERM?
        print('calculating LTIR using 24, 100, and 160 micron data')

        # coefficients from Table 3 in Galametz+ 2013
        c100 = 1.239
        c160 = 0.620
        
        # read in each image
        f100 = fits.open(b100)
        f160 = fits.open(b160)
        
        # regrid data to same grid (use 24micron as base) and 
        # convert from MJy/sr to W/kpc^2
        basehdr = f100[0].header

        data100freq = c /(100.0*1e-4)
        data100 = f100[0].data * unitfactor * data100freq
        

        data160freq = c/(160.0*1e-4)
        data160_regrid = reproject_interp(f160,output_projection=basehdr,
                                          hdu_in=0,
                                          order='nearest-neighbor',
                                          return_footprint=False)
        data160_regrid = data160_regrid * unitfactor * data160freq
        

        f100.close()
        f160.close()

        # calculate LTIR using coefficients from Table 3 in Galametz+ 213
        STIR_Wkpc2 =  c100*data100 + c160*data160_regrid

        # convert from W/kpc^2 to Lsun/pc^2 (which is what we typically plot)
        # W/kpc^2 * 1e7 = erg/s/kpc^2
        # erg/s/kpc^2 / Lsun / (1000.0)**2 = Lsun/pc^2
        STIR_Lsunpc2 = (STIR_Wkpc2 * 1e7) / Lsun / (1000.0**2) 

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    elif (b24): ## monochromatic. :-(
        
        print('calculating LTIR using 24 micron data')

        # coefficients from Table 2 in Galametz+ 2013
        a = 0.919
        b = 3.786

        # read in each image
        f24 = fits.open(b24)
        
        # convert from MJy/sr to W/kpc^2
        basehdr = f24[0].header

        data24freq = c / (24.0*1e-4)
        data24 = f24[0].data * unitfactor * data24freq

        f24.close()   

         # calculate LTIR using coefficients from Table 3 in Galametz+ 213
        STIR_Wkpc2 = 10**(a * np.log10(data24) + b)

        # convert from W/kpc^2 to Lsun/pc^2 (which is what we typically plot)
        # W/kpc^2 * 1e7 = erg/s/kpc^2
        # erg/s/kpc^2 / Lsun / (1000.0)**2 = Lsun/pc^2
        STIR_Lsunpc2 = (STIR_Wkpc2 * 1e7) / Lsun / (1000.0**2) 

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    elif (w4): 

        print('calculating LTIR using WISE band 4 data')

        ## use wise band 4 to calculate LTIR.  

        ## Cluver+2017 bootstraps off the SINGS/KINGFISH data to show
        ## that the relationship between LTIR and L12micron(W3) is
        ## tighter and more linear than LTIR and L23micron (w4)
        
        ## Equation in Cluver+ 2017 
        ## (Figure 3, equation 2; See erratum for correct equation.).

        ## TODO: the equation is in luminosity, so  I think I need to 
        ## convert t surface brightness. Do this per pixel or per beam?
        a = 0.915
        b = 1.96

        # read in each image
        f23 = fits.open(w4)
        
        basehdr = f23[0].header

        ## input unit is MJy/sr for W3 data, which is a surface brightness.
        # convert from MJy/sr to W/kpc^2, then to Lsun/pc^2
        data23freq = c / (23.0*1e-4)
        data23 = f23[0].data * unitfactor * data23freq
        data23 = (data23 * 1e7) / Lsun / (1000.0**2) 

        f23.close()   

        # calculate LTIR using coefficients above
        STIR_Lsunpc2 = 10**(a * np.log10(data23) + b)

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    elif (w3): 

        print('calculating LTIR using WISE band 3 data')

        ## use wise band 3 to calculate LTIR.  

        ## Cluver+2017 bootstraps off the SINGS/KINGFISH data to show
        ## that the relationship between LTIR and L12micron(W3) is
        ## tighter and more linear than LTIR and L23micron (w4)
        
        ## Equation in Cluver+ 2017 
        ## (Figure 3, equation 1; See erratum for correct equation.).
        a = 0.889
        b = 2.21

        # read in each image
        f12 = fits.open(w3)
        
        basehdr = f12[0].header

        ## input unit is MJy/sr for W3 data, which is a surface brightness.
        # convert from MJy/sr to W/kpc^2, then to Lsun/pc^2
        data12freq = c / (12.0*1e-4)
        data12 = f12[0].data * unitfactor * data12freq
        data12 = (data12 * 1e7) / Lsun / (1000.0**2) 

        f12.close()   

         # calculate LTIR using coefficients above
        STIR_Lsunpc2 = 10**(a * np.log10(data12) + b)

        # write the output file
        outHdr = basehdr
        outHdr['BUNIT'] = 'Lsun/pc^2'
        outfits = fits.PrimaryHDU(data=STIR_Lsunpc2, header=outHdr)
        outfits.writeto(outFile,overwrite=True)

    else:
        print("Inappropriate combination of images to calculate LTIR.")
        
#----------------------------------------------------------------------

def colorCorrect_24micron(b24, b70, outFile, beam=30.0):
    '''

    color correct the 24micron only data with lower resolution
    70micron data.

    Date        Programmer      Description of Code
    ----------------------------------------------------------------------
    5/5/2022    A.A. Kepley     Original Code
    '''

    # get some names set up
    b70_dir = os.path.dirname(b70)
    b24_name = os.path.basename(b24)
    galname = b24_name.split('_')[0]

    # smooth 24micron data
    # TODO -- determine the beam to smooth to
    # from the 70micron data directly
    b24_smoothed = smoothCube(b24,b70_dir, beam=beam)
        
    # calculate LTIR of smoothed data
    LTIRFile = os.path.join(b70_dir, galname+'_LTIR_b24_b70.fits')
    calc_LTIR(LTIRFile, b24=b24_smoothed, b70=b70)

    # get necessary data
    LTIR = fits.open(outFile)[0].data
    b24_data = fits.open(b24)[0].data
    b24_hdr = fits.open(b24)[0].header
    b24_smoothed_data = fits.open(b24_smoothed)[0].data
    
    # calculate the 24micron to LTIR ratio for the smoothed data
    ratio = b24_smoothed_data / LTIR
    
    # correct the unsmoothed 24micron data
    LTIR_corrected = b24_data / ratio

    b24_hdr['BUNIT'] = 'Lsun/pc^2'

    # write out the results
    fits.writeto(outFile, LTIR_corrected, header=b24_hdr, overwrite=True)

#----------------------------------------------------------------------

def smoothCube(fitsimage,outDir, beam=15.0):
    '''
    Smoothes input cube to given resolution
    
    beam is the beam to smooth to in arcsec.

    returns name of new fits file

    '''

    from radio_beam import Beam

    # determine ndim and bunit for other data sets
    f = fits.open(fitsimage)
    ndim = f[0].header['NAXIS'] # has to have this, so not worth checking
    
    if 'BUNIT' in f[0].header:
        bunit = f[0].header['BUNIT']
    else:
        bunit = None

    f.close()

    newFits = os.path.join(outDir,os.path.basename(fitsimage).replace(".fits","_smooth.fits"))

    newBeam = Beam(beam*u.arcsec)

    if ndim == 3:
        cube = SpectralCube.read(fitsimage)
    elif ndim == 2:
        vhdu = fits.open(fitsimage)

        if bunit:
            bunit = vhdu[0].header['BUNIT']
            # fix up bunits
            bunit = bunit.lower().replace('msun','Msun') 
            bunit = bunit.replace('mjy','MJy')

        ## TODO -- switch from spectral cube to straight astropy?
        cube = Projection.from_hdu(vhdu)
        if (cube.unit == u.Unit("")) & bool(bunit):
            cube = cube * u.Unit(bunit)
    else:
        print("ndim not 2 or 3")
        return

    smoothCube = cube.convolve_to(newBeam)

    #ipdb.set_trace()

    smoothCube.write(newFits,overwrite=True)
    return newFits

#----------------------------------------------------------------------

def regridData(baseCubeFits, otherDataFits, outDir, mask=False):
    '''
    regrids one data set to match the wcs of the base data set, which
    is assumed to be a cube. The regridded data set can be either 2d
    or 3d.
    '''

    # open the base cube
    try:
        baseCube = SpectralCube.read(baseCubeFits)
    except:
        print("Can't read in " + baseCubeFits + ".")

    # determine ndim and bunit for other data sets
    f = fits.open(otherDataFits)
    ndim = f[0].header['NAXIS']
    bunit = f[0].header['BUNIT']
    # TODO: save beam information here for input into 2D header
    f.close()

    # output image name
    newFits = os.path.join(outDir,os.path.basename(otherDataFits).replace('.fits','_regrid.fits'))

    # now regrid images appropriately
    if ndim == 3:

        # read in cube
        otherCube = SpectralCube.read(otherDataFits)

        # interpolate velocity axis. This needs to be done first.
        regridCube = otherCube.spectral_interpolate(baseCube.spectral_axis)

        regridCube.allow_huge_operations=True
        newCube = regridCube.reproject(baseCube.header)

        if mask:
            # if greater than 0.0 set value to 1. otherwise 0.
            newdata = np.where(newCube.unitless_filled_data[:,:,:] > 0.0, 
                               1, 0)  
            finalCube = SpectralCube(newdata,newCube.wcs,
                                     mask=baseCube.mask)
        
        else:
            newdata = newCube.filled_data[:,:,:]
            finalCube = SpectralCube(newdata,newCube.wcs,mask=baseCube.mask)
        
        finalCube.with_spectral_unit(baseCube.spectral_axis.unit).with_beam(otherCube.beam).write(newFits,overwrite=True)
        
        return newFits

    elif ndim == 2:
        
        ## consider using Projection.from_hdu then Projection.reproject here?

        # regrid image
        newcube = reproject_interp(otherDataFits,
                                   output_projection=baseCube.wcs.dropaxis(2),
                                   shape_out=baseCube.wcs.dropaxis(2).array_shape, 
                                   order='nearest-neighbor',                     
                                   return_footprint=False)        

        if mask:
            # set to 1 where greater than 0.0.
            newdata = np.where(newcube > 0.0,1.0,0.0)
        else:
            newdata = newcube

        # get mask on original data
        baseMask = baseCube.get_mask_array()
        totalBaseMask = baseMask.max(axis=0)
        
        newdata[np.invert(totalBaseMask)] = np.nan
        
        # fix up header with input image units
        outHdr = baseCube.wcs.dropaxis(2).to_header()
        outHdr.append(('BUNIT',bunit))
        ## TODO -- add beam info from otherDataFits

        # write out regridded data
        fits.writeto(newFits,newdata,outHdr,overwrite=True)
        
        return newFits

    else:
        print("Number of dimensions of other data set is not 2 or 3.")

def simpleR21scale(infile, r21, r21_ref=None):
    '''

    do a simple scaling of 12CO(2-1) data based a single R21 values.

    Date        Programmer      Description of Code
    ----------------------------------------------------------------------
    10/20/2021  A.A. Kepley     Original Code

    '''
    
    hdu = fits.open(infile)
    hdr = hdu[0].header
    data_rescaled = hdu[0].data/r21

    hdr['R21'] = (r21, r21_ref)

    fits.writeto(infile.replace('12CO21','12CO10_r21_simple'),data_rescaled, header=hdr, overwrite=True)

    hdu.close()


def sfrR21scale(infile, sigmaSFRfile, galaxy=None,
                re=None, sfr=None, r21=None, r21_ref=None,
                r21_lim=0.4):
    '''

    scale the R21 value by the sigma_SFR value using relationship from
    Leroy+2021.

    Date        Programmer              Description of Code
    ----------------------------------------------------------------------
    10/20/21    A. Kepley               Original Code
    '''
    
    import math

    hdu = fits.open(infile)
    hdr_in = hdu[0].header
    data_in = hdu[0].data
    hdu.close()

    hdu = fits.open(sigmaSFRfile)
    hdr_sigmaSFR = hdu[0].header
    sigmaSFR = hdu[0].data
    hdu.close()

    hdr_in['R21'] = (r21, r21_ref)

    sigmaSFR_avg = 0.5 * sfr / ( math.pi * re**2)

    Q = np.log10(sigmaSFR / sigmaSFR_avg)

    # coefficients from Leroy+ 2021 Table 5
    a = 0.129
    b = 0.019

    logR21ratio = a * Q + b

    logR21 = logR21ratio + np.log10(r21)

    R21_scaled = 10**logR21

    R21_scaled[np.where(R21_scaled < r21_lim)] = r21_lim

    data_scaled = data_in / R21_scaled

    indir = os.path.dirname(infile)
    hdr_out = hdr_in
    hdr_out['BUNIT'] = ''

    fits.writeto(os.path.join(indir,galaxy+"_R21.fits"),R21_scaled,header=hdr_out,overwrite=True)

    fits.writeto(infile.replace("12CO21","12CO10_r21_sigmaSFR"),data_scaled, header=hdr_in, overwrite=True)

    
