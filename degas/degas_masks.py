from spectral_cube import SpectralCube
import astropy.units as u
from scipy.ndimage import map_coordinates
import numpy as np
import copy
import astropy.wcs as wcs
from astropy.io import fits

def buildmasks(filename, nChan=2000, width=2e9):
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

    """


    # Read in original cube, ensure in velocity space
    s = SpectralCube.read(filename)
    s = s.with_spectral_unit(u.km / u.s, velocity_convention='radio')
    
    # HCN_HCO+
    # Build a mask with a spectral width of 2 GHz and the same spatial 
    # dimensions as the original mask
  
    s_hcn = s.with_spectral_unit(u.Hz, rest_value=88.631847 * u.GHz)
    s_hcop = s.with_spectral_unit(u.Hz, rest_value=89.188518 * u.GHz)

    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=np.byte)
    hdr = s_hcn.wcs.to_header()
    hdr['CRPIX3'] = 1000
    hdr['CDELT3'] = width / nChan
    hdr['CRVAL3'] = (89.188518 + 88.631847) / 2 * 1e9
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = mask.shape[0]
    hdr['NAXIS2'] = mask.shape[1]
    hdr['NAXIS3'] = mask.shape[2]
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 8
    hdr['EXTEND'] = 'T'
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
                                      dtype=np.bool)
        if 0 <= zz_hcop < s_hcop.shape[0]:
            mask[zz, :, :] = np.array(s_hcop.filled_data[zz_hcop, :, :],
                                      dtype=np.bool)
    maskcube = SpectralCube(mask, w, header=hdr)
    maskcube.write(filename.replace('.fits', '.hcn_hcop_mask.fits'),
                   overwrite=True)

    # C18O/13CO
    # Build a mask with a spectral width of 2 GHz and the same spatial
    # dimensions as the original mask

    s_13co = s.with_spectral_unit(u.Hz, rest_value=110.20135 * u.GHz)
    s_c18o = s.with_spectral_unit(u.Hz, rest_value=109.78217 * u.GHz)

    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=np.byte)
    hdr = s_13co.wcs.to_header()
    hdr['CRPIX3'] = 1000
    hdr['CDELT3'] = width / nChan
    hdr['CRVAL3'] = (110.20135 + 109.78217) / 2 * 1e9
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = mask.shape[0]
    hdr['NAXIS2'] = mask.shape[1]
    hdr['NAXIS3'] = mask.shape[2]
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 8
    hdr['EXTEND'] = 'T'
    w = wcs.WCS(hdr)
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
                                      dtype=np.bool)
        if 0 <= zz_c18o < s_c18o.shape[0]:
            mask[zz, :, :] = np.array(s_c18o.filled_data[zz_c18o, :, :],
                                      dtype=np.bool)
    maskcube = SpectralCube(mask, w, header=hdr)
    maskcube.write(filename.replace(
        '.fits', '.13co_c18o_mask.fits'), overwrite=True)

    # 12CO
    # Build a mask with a spectral width of 2 GHz and the same spatial
    # dimensions as the original mask

    s_12co = s.with_spectral_unit(u.Hz, rest_value=115.271204 * u.GHz)


    mask = np.zeros((nChan, s.shape[1], s.shape[2]), dtype=np.byte)
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
    maskcube = SpectralCube(mask, w, header=hdr)
    for zz in range(nChan):
        nu = maskcube.spectral_axis[zz]
        _, _, zz_12co = s_12co.wcs.wcs_world2pix(hdr['CRVAL1'],
                                               hdr['CRVAL2'],
                                               nu, 0)
        zz_12co = int(zz_12co)
        if 0 <= zz_12co < s_12co.shape[0]:
            mask[zz, :, :] = np.array(s_12co.filled_data[zz_12co, :, :],
                                      dtype=np.bool)
    maskcube = SpectralCube(mask, w, header=hdr)
    maskcube.write(filename.replace(
        '.fits', '.12co_mask.fits'), overwrite=True)
