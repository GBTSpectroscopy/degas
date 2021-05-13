from spectral_stack import stacking #need to make sure Erik's code is installed
from spectral_cube import SpectralCube, Projection
import numpy as np
import numpy.ma as ma
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import colors
import aplpy
from astropy.table import Table, Column, vstack, QTable
import re
import glob 
from astropy.wcs import WCS
import os

def makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1', outname='test', release='DR1', sourceList=None):
    '''

    Calculate stacking results for each galaxy and make a fits table
    of the results.

    Inputs:

        regridDir: directory containing regridded data

        scriptDir: directory with DEGAS data base

        outDir: directory for output results.

        vtype: velocity type. Options are 'mom1' and 'peakVelocity'. Default: 'mom1'
    
        outname: basename for output fits table.


    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    10/30/2020  A.A. Kepley     Added comments and clarified inputs
    12/3/2020   A.A. Kepley     Modified to pull list of galaxies from 
                                degas base table.

    '''
    
    # get list of dr1 galaxies
    degas = Table.read(os.path.join(scriptDir,'degas_base.fits'))
    idx = degas[release] == 1

    if not sourceList:
        sourceList = degas[idx]['NAME']

    # go though each file, create a table, and attend it to the list of tables.
    tablelist=[]
    for galaxy in degas[idx]:
        if galaxy['NAME'] in sourceList:
            full_tab = makeGalaxyTable(galaxy, vtype, regridDir, outDir)
            tablelist.append(full_tab)

    # stack all the tables together
    table=vstack(tablelist)

    # Write out the table 
    table.write(os.path.join(outDir,outname+'_'+vtype+'.fits'),overwrite=True)

    return table

def makeGalaxyTable(galaxy, vtype, regridDir, outDir):
    '''

    make fitstable containing all lines from stacking results for each galaxy

    galaxy: data for galaxy we are processing from degas_base.fits

    vtype: velocity type that we are stacking on
    
    scriptDir: script directory -- AAK: if I read in degas_base.fits early,
    so I still need to pass this down.?
    
    regridDir: regrid directory
    x
    outDir: output directory

    
    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added comments
    4/29/2021   A.A. Kepley     Moved all calculations for galaxy here 
                                instead of repeating per line.
    5/6/2021    A.A. Kepley     Modified so that all lines are calculated at once.

    '''

    print("Processing " + galaxy['NAME'] + "\n")

    # Create associated maps needed for analysis.
    mom0cut, cubeCO, sn_mask = mapSN(galaxy, regridDir, outDir, sncut=3.0)  
    stellarmap = mapStellar(galaxy, mom0cut, regridDir, outDir)
    sfrmap = mapSFR(galaxy,mom0cut, regridDir, outDir)
    ltirmap = mapLTIR(galaxy,mom0cut,regridDir,outDir)
    R_arcsec, R_kpc, R_r25 = mapGCR(galaxy,mom0cut)

    velocity_file = galaxy['NAME']+'_12CO_'+vtype+'_regrid.fits' 
    vhdu = fits.open(os.path.join(regridDir,velocity_file))
    velocity = Projection.from_hdu(vhdu)

    # read in HCN
    linefile = os.path.join(regridDir,galaxy['NAME']+'_HCN_rebase3_smooth1.3_hanning1_maxnchan_smooth.fits')
    cubeHCN = SpectralCube.read(os.path.join(regridDir,linefile),mask=sn_mask)
    
    # read in HCO+
    linefile = os.path.join(regridDir,galaxy['NAME']+'_HCOp_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
    cubeHCOp = SpectralCube.read(os.path.join(regridDir,linefile),mask=sn_mask)

    # For NGC6946, skip 13CO and C18O since we don't have that data.
    if galaxy['NAME'] != 'NGC6946':
        # read in 13CO
        linefile = os.path.join(regridDir,galaxy['NAME']+'_13CO_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
        cube13CO = SpectralCube.read(os.path.join(regridDir,linefile),mask=sn_mask)

        # read in C18O
        linefile = os.path.join(regridDir,galaxy['NAME']+'_C18O_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
        cubeC18O = SpectralCube.read(os.path.join(regridDir,linefile),mask=sn_mask)
        
    else:
        cube13CO = None
        cubeC18O = None

    #get the full stack result for each line
    full_stack = makeStack(galaxy, regridDir, outDir,
                                       mom0cut = mom0cut, 
                                       cubeCO = cubeCO,
                                       cubeHCN = cubeHCN, cubeHCOp=cubeHCOp,
                                       cube13CO = cube13CO, cubeC18O = cubeC18O,
                                       velocity = velocity,
                                       sfrmap = sfrmap, ltirmap = ltirmap, 
                                       stellarmap=stellarmap, R_arcsec=R_arcsec) 
    # remove stacks that don't have CO spectra 
    nstack = len(full_stack)
    keepstack = np.full(nstack,True)
    
    for i in range(nstack):
        if np.all(full_stack['stack_profile_CO'][i] == 0):
            keepstack[i] = False

    full_stack = full_stack[keepstack]
        

    full_stack = addIntegratedIntensity(full_stack)

    # return the table and the stack.
    return full_stack


def makeStack(galaxy, regridDir, outDir,
              mom0cut = None,
              cubeCO = None,
              cubeHCN = None, cubeHCOp=None,
              cube13CO = None, cubeC18O = None,
              velocity = None, sfrmap = None, ltirmap = None,
              stellarmap = None, R_arcsec = None):
    '''
    
    make stacks for all lines and ancillary data for one galaxy
    output python dictionaries

    galaxy: degas_base.fits entry for galaxy we are processing
    
    line: line we are stacking on
    
    regridDir: directory with all the regridded data

    outDir: directory with all the output data.

    ASSUMPTIONS: 

    We are assuming that all linecubes have been smoothed 
    and regridded at this point, so the coordinate systems match.


     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added more comments plus minor changes to 
                                use the degas table 
    12/10/2020  A.A. Kepley     Added LTIR
    4/29/2021   A.A. Kepley     Moved up map calculations to makeStack to 
                                avoid doing multiple times.

    '''
 

    ## create intensity bins
    binmap, binedge, binlabels = makeBins(galaxy, mom0cut, 'intensity', outDir)  
    plotBins(galaxy, binmap, binedge, binlabels, 'intensity', outDir) 

    # do intensity stack
    # TODO CHECK UNITS
    #r_intensity=stackLines(galaxy, velocity, 
    #                  binmap, binedge, 'intensity', 'K km/s',
    #                  cubeCO = cubeCO,
    #                  cubeHCN = cubeHCN, cubeHCOp = cubeHCOp,
    #                  cube13CO = cube13CO, cubeC18O = cubeC18O,
    #                  sfrmap = sfrmap, ltirmap = ltirmap) 

    ## create stellar mass bin
    binmap, binedge, binlabels = makeBins(galaxy, stellarmap, 'stellarmass', outDir)  
    plotBins(galaxy, binmap, binedge, binlabels, 'stellarmass', outDir) 

    ## do stellar mass stack
    ## TODO CHECK UNITS
    #r_stellarmass=stackLines(galaxy, velocity,  
    #                    binmap, binedge,  'stellarmass', 'Msun/pc^2',
    #                    cubeCO = cubeCO,
    #                    cubeHCN = cubeHCN, cubeHCOp = cubeHCOp,
    #                    cube13CO = cube13CO, cubeC18O = cubeC18O,
    #                    sfrmap = sfrmap, ltirmap = ltirmap)

    ## create radius bins
    binmap, binedge, binlabels = makeRadiusBins(galaxy, R_arcsec, outDir) 
    plotBins(galaxy, binmap, binedge, binlabels, 'radius', outDir) 

    # stack on radius
    ## TODO CHECK UNITS
    r_radius=stackLines(galaxy, velocity,
                        binmap, binedge,  'radius',  'arcsec',
                        cubeCO = cubeCO,
                        cubeHCN = cubeHCN, cubeHCOp = cubeHCOp,
                        cube13CO = cube13CO, cubeC18O = cubeC18O,
                        sfrmap = sfrmap, ltirmap = ltirmap)
    r_radius.add_column(Column(np.tile('radius',len(r_radius))),name='bintype')

    ## TODO ADD OTHER STACKS ONCE I HAVE BINS
    # put the individual stack together. 
    #full_stack = vstack([r_intensity, r_stellarmass, r_radius])
    full_stack = r_radius
    
    full_stack.add_column(Column(np.tile(galaxy['NAME'],len(full_stack))),name='galaxy')

    return full_stack

def stackLines(galaxy, velocity, 
               binmap, binedge, bintype, binunit,
               cubeCO = None,
               cubeHCN = None,
               cubeHCOp = None,
               cube13CO = None,
               cubeC18O = None,
               sfrmap = None, 
               ltirmap = None,
               maxAbsVel = 250.0):

    '''
    Actually do the stacking.

    cube: data we are stacking (SpectralCube)

    galaxy: line from degas_base.fits table with galaxy information.
    
    velocity: velocity map to stack data using

    binmap: bin map

    binedge: bin edges
    
    binlabel: bin labels
    
    bintype: type of bin
    
    cubeCO: CO cube

    cubeHCN: HCN cube
    
    cubeHCOp: HCO+ cube
    
    cube13CO: 13CO cube
    
    cubeC18O: C18O cube

    sfrmap: star formation rate map

    ltirmap: LTIR map

    maxAbsVel: absolute maximum velocity to include in spectral in km/s. 
    Assumes line is centered at zero (which should be true for stacking).
    
    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/03/2020  A.A. Kepley     Added comments and moved bin creation up 
                                a level to simplify code.
    12/10/2020  A.A. Kepley     added LTIR calculation
    5/6/2021    A.A. Kepley     modified to fit all lines for one galaxy at 
                                once so I can use CO FWHM to calculate upper 
                                limits for other lines

    '''

    # get the relevant info on the galaxy and cube
    D = galaxy['DIST_MPC'] * u.Mpc
    pix_area = (np.radians(np.abs(cubeCO.header['CDELT1']))*D.to(u.kpc))**2  #kpc^2 ## TODO CHECK THIS CONVERSION

    # do the individual line stacks.
    # -----------------------------

    stack={}
    if cubeCO:
        ## fill in nans with 0 to avoid fft shift issue with nans
        cubeCO_nonan = cubeCO.with_fill_value(0) 
        
        # do the stack -- making sure the units work out right.

        stack['CO'], labelvals = stacking.BinByLabel(cubeCO_nonan,
                                                     binmap.value, velocity,
                                                     weight_map=None)
     
    if cubeHCN:
        ## fill in nans with 0 to avoid fft shift issue with nans
        cubeHCN_nonan = cubeHCN.with_fill_value(0) 
        
        # do the stack -- making sure the units work out right.
        stack['HCN'], labelvals = stacking.BinByLabel(cubeHCN_nonan,
                                                      binmap.value, velocity,
                                                      weight_map=None)

    if cubeHCOp:

        cubeHCOp_nonan = cubeHCOp.with_fill_value(0) 
        
        # do the stack -- making sure the units work out right.
        stack['HCOp'], labelvals = stacking.BinByLabel(cubeHCOp_nonan,
                                                       binmap.value, velocity,
                                                       weight_map=None)
        
    if cube13CO:
        cube13CO_nonan = cube13CO.with_fill_value(0) 
        
        # do the stack -- making sure the units work out right.
        stack['13CO'], labelvals = stacking.BinByLabel(cube13CO_nonan,
                                                       binmap.value, velocity,
                                                       weight_map=None)

    if cubeC18O:
        
        cubeC18O_nonan = cubeC18O.with_fill_value(0) 
        
        # do the stack -- making sure the units work out right.
        stack['C18O'], labelvals = stacking.BinByLabel(cubeC18O_nonan,
                                                       binmap.value, velocity,
                                                       weight_map=None)

    # putting the table together
    # -------------------------

    # first add bins
    t = {'bin_lower': binedge[0:-1], 
         'bin_upper': binedge[1:], 
         'bin_mean': (binedge[0:-1] + binedge[1:])/2.0}
    total_stack = QTable(t)
    
    # Then set up the structures for the bin-based profiles, etc.
    spectral_profile = {}
    for line in ['CO','HCN','HCOp','13CO','C18O']:
        spectral_profile[line] = np.zeros((len(labelvals),len(stack['CO'][0]['spectral_axis']))) 
 
    bin_label = np.zeros(len(labelvals))
    bin_area = np.zeros(len(labelvals)) * pix_area.unit
    sfr_mean = np.zeros(len(labelvals)) * sfrmap.unit
    ltir_mean = np.zeros(len(labelvals))* ltirmap.unit

    spectral_axis = np.zeros((len(labelvals),len(stack['CO'][0]['spectral_axis']))) * stack['CO'][0]['spectral_axis'].unit

    for i in range(len(labelvals)):
        spectral_axis[i,:] = stack['CO'][i]['spectral_axis']

        bin_label[i] = stack['CO'][i]['label']

        ## TODO CHECK THE CALCULATION BELOW FOR BINAREA
        bin_area[i]= float(sum(binmap[binmap==bin_label[i]].flatten()))*pix_area          
        # calculating the mean SFR as requested
        if sfrmap is not None:
            sfr_mean[i] = np.nanmean(sfrmap[binmap==bin_label[i]]) 
        # calculate the mean LTIR as request
        if ltirmap is not None:
            ltir_mean[i] = np.nanmean(ltirmap[binmap==bin_label[i]]) 
    
        # get the spectral profiles
        for line in ['CO','HCN','HCOp','13CO','C18O']:
            if line in stack.keys():
                spectral_profile[line][i,:] = stack[line][i]['spectrum'] 
            else:
                 spectral_profile[line][i,:] = np.full(len(stack['CO'][0]['spectral_axis']),np.nan)

    # add above items to the table
    total_stack.add_column(Column(spectral_axis),name='spectral_axis')
    for line in ['CO','HCN','HCOp','13CO','C18O']:
        total_stack.add_column(Column(spectral_profile[line], name='stack_profile_'+line,unit=cubeCO.unit))

    # TODO CHECK UNITS HERE
    total_stack.add_column(Column(bin_area,name='bin_area'))
    total_stack.add_column(Column(sfr_mean,name='sfr_mean')) # Msun/yr/kpc^2 -- units from input image
    total_stack.add_column(Column(sfr_mean * bin_area),name='sfr_total') # Msun/yr/kpc^2 * kpc^2 = Msun/Yr
    total_stack.add_column(Column(ltir_mean),name='ltir_mean') # Lsun/pc^2 -- units from input image
    total_stack.add_column(Column(ltir_mean * 1000.0**2 *bin_area),name='ltir_total') # Lsun/pc^2 * (1000 pc/kpc)^2 * kpc^2   = Lsun
        
    return total_stack

def addIntegratedIntensity(full_stack):
    '''
    calculate the integrated line flux and the associated noise

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    5/13/2021   A.A. Kepley     Original Code

    '''

    # get number of stacks
    nstack = len(full_stack)

    # initialize output arrays
    int_intensity_fit = {}
    int_intensity_fit_err = {}
    int_intensity_fit_uplim = {}

    int_intensity_sum = {}
    int_intensity_sum_err = {}
    int_intensity_sum_uplim = {}

    fwhm_fit = np.zeros(nstack)* full_stack['spectral_axis'].unit
    fwhm_sum = np.zeros(nstack)* full_stack['spectral_axis'].unit

    for line in ['CO','HCN','HCOp','13CO','C18O']:

        int_intensity_fit[line] = ma.zeros(nstack) * full_stack['stack_profile_CO'].unit * full_stack['spectral_axis'].unit
        int_intensity_fit_err[line] = ma.zeros(nstack)* full_stack['stack_profile_CO'].unit * full_stack['spectral_axis'].unit
        int_intensity_fit_uplim[line] = np.full(nstack,False)

        int_intensity_sum[line] = ma.zeros(nstack) * full_stack['stack_profile_CO'].unit * full_stack['spectral_axis'].unit
        int_intensity_sum_err[line] = ma.zeros(nstack) * full_stack['stack_profile_CO'].unit * full_stack['spectral_axis'].unit
        int_intensity_sum_uplim[line] = np.full(nstack,False)
        
    # calculate the integrated intensity from a fit and from a simple sum.
    for i in range(nstack):
        
        # gaussian fit
        (stack_int, stack_int_err, fwhm, uplim) = fitIntegratedIntensity(full_stack[i]['spectral_axis'], 
                                                                   full_stack[i]['stack_profile_CO'])

        int_intensity_fit['CO'][i] = stack_int
        int_intensity_fit_err['CO'][i] = stack_int_err
        int_intensity_fit_uplim[i] = uplim
        fwhm_fit[i] = fwhm

        for line in ['HCN', 'HCOp', '13CO','C18O']:

            if ( ( full_stack[i]['galaxy'] == 'NGC6946') & 
                 ( ( line == '13CO') | (line == 'C18O'))):
            
                int_intensity_fit[line][i] = np.nan
                int_intensity_fit_err[line][i] = np.nan
                int_intensity_fit_uplim[line][i] = True

            else:
            
                (stack_int, stack_int_err, fwhm, uplim) = fitIntegratedIntensity(full_stack[i]['spectral_axis'], 
                                                                                 full_stack[i]['stack_profile_'+line], fwhm = fwhm)

                int_intensity_fit[line][i] = stack_int
                int_intensity_fit_err[line][i] = stack_int_err
                int_intensity_fit_uplim[line][i] = uplim

            
        # straight sum
        (stack_sum, stack_sum_err, fwhm, uplim) = sumIntegratedIntensity(full_stack[i]['spectral_axis'], 
                                                                   full_stack[i]['stack_profile_CO'])

        int_intensity_sum['CO'][i] = stack_sum
        int_intensity_sum_err['CO'][i] = stack_sum_err
        int_intensity_sum_uplim[i] = uplim
        fwhm_sum[i] = fwhm

        for line in ['HCN', 'HCOp', '13CO','C18O']:
            
            if ( ( full_stack[i]['galaxy'] == 'NGC6946') & 
                 ( ( line == '13CO') | (line == 'C18O'))):
                
                int_intensity_sum[line][i] = np.nan
                int_intensity_sum_err[line][i] = np.nan
                int_intensity_sum_uplim[line][i] = True

            else:
            
            
                (stack_sum, stack_sum_err, fwhm, uplim) = sumIntegratedIntensity(full_stack[i]['spectral_axis'], 
                                                                                 full_stack[i]['stack_profile_'+line], 
                                                                             fwhm = fwhm)

                int_intensity_sum[line][i] = stack_sum
                int_intensity_sum_err[line][i] = stack_sum_err
                int_intensity_sum_uplim[line][i] = uplim


    # add results to data able.
    for line in ['CO','HCN','HCOp','13CO','C18O']:
        full_stack.add_column(Column(int_intensity_fit[line],name='int_intensity_fit_'+line))
        full_stack.add_column(Column(int_intensity_fit_err[line],name='int_intensity_fit_err_'+line))
        full_stack.add_column(Column(int_intensity_fit_uplim[line],name='int_intensity_fit_uplim_'+line))

        full_stack.add_column(Column(int_intensity_sum[line],name='int_intensity_sum_'+line))
        full_stack.add_column(Column(int_intensity_sum_err[line],name='int_intensity_sum_err_'+line))
        full_stack.add_column(Column(int_intensity_sum_uplim[line],name='int_intensity_sum_uplim_'+line))
        

    full_stack.add_column(Column(fwhm_fit),name='FWHM_fit')
    full_stack.add_column(Column(fwhm_sum),name='FWHM_sum')

    return full_stack


def sumIntegratedIntensity(spectral_axis, stack_profile, fwhm=None, maxAbsVel=250.0*u.km/u.s, snThreshold=3.0):
    '''
    calculate the straight sum of the integrated intensity.
    
    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    5/13/2021   A.A. Kepley     Original Code

    '''

    from astropy.modeling import models, fitting
    from scipy import integrate

    # default is to scale to normal distribution
    from scipy.stats import median_absolute_deviation as mad 

    chanwidth = spectral_axis[1] - spectral_axis[0]

    lineFreeChans = (spectral_axis > maxAbsVel ) | (spectral_axis < - maxAbsVel)

    # mad is already scaled to gaussian distribution
    noisePerChan = mad(stack_profile[lineFreeChans]) * stack_profile.unit

    lineChans = (spectral_axis < maxAbsVel ) & (spectral_axis > - maxAbsVel)
    
    if np.any(stack_profile[lineChans] > snThreshold * noisePerChan):
        # sum line

        # start off by fitting one Gassian
        amp_est = max(stack_profile)
        peak_cut = spectral_axis[stack_profile > amp_est * 0.5]
        fwhm = max(peak_cut) - min(peak_cut)
        sigma_est = fwhm/2.355
        init_g = models.Gaussian1D(amplitude = amp_est, stddev=sigma_est)
        fit_g = fitting.LevMarLSQFitter()
        result_g = fit_g(init_g, spectral_axis, stack_profile)

        newLineChans = ( (spectral_axis <  (result_g.mean + 3.0*result_g.stddev)) & 
                         (spectral_axis > (result_g.mean - 3.0*result_g.stddev)))
        
        stack_sum = np.sum(stack_profile[newLineChans]*chanwidth)
        stack_sum_err = np.sqrt(fwhm/chanwidth) * chanwidth * noisePerChan

        uplim=False

    elif fwhm:
        stack_sum_err = np.sqrt(fwhm/chanwidth) * chanwidth * noisePerChan
        stack_sum = snThreshold * stack_sum_err
        uplim = True

    else:
        stack_sum_err = np.nan * stack_profile.unit * spectral_axis.unit
        stack_sum = np.nan * stack_profile.unit * spectral_axis.unit
        fwhm = np.nan * spectral_axis.unit
        uplim = True
        
    return stack_sum, stack_sum_err, fwhm, uplim
        

def fitIntegratedIntensity(spectral_axis,stack_profile, fwhm=None, maxAbsVel=250 * u.km/u.s, snThreshold=3.0):
    '''
    
    calculate integrated intensity via a gaussian fit

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    5/13/2021   A.A. Kepley     Original Code

    '''
    
    # default is to scale to normal distribution
    from scipy.stats import median_absolute_deviation as mad 

    chanwidth = spectral_axis[1] - spectral_axis[0]

    lineFreeChans = (spectral_axis > maxAbsVel ) | (spectral_axis < - maxAbsVel)

    # mad is already scaled to gaussian distribution
    noisePerChan = mad(stack_profile[lineFreeChans]) * stack_profile.unit

    lineChans = (spectral_axis < maxAbsVel ) & (spectral_axis > - maxAbsVel)
    
    if np.any(stack_profile[lineChans] > snThreshold * noisePerChan):
        # fit line using Gaussian

        (stack_int, stack_int_err, fwhm) = fitGaussian(spectral_axis, stack_profile, noisePerChan)

        uplim = False

    elif fwhm:

        stack_int_err = np.sqrt(fwhm/chanwidth) * chanwidth * noisePerChan
        
        stack_int = snThreshold * stack_int_err

        uplim = True

    else:
        stack_int = np.nan * spectral_axis.unit * stack_profile.unit
        stack_int_err = np.nan* spectral_axis.unit * stack_profile.unit
        fwhm = np.nan* spectral_axis.unit 
        uplim = True

    return stack_int, stack_int_err, fwhm, uplim


def fitGaussian(spectral_axis, stack_profile, noisePerChan, maxAbsVel=250*u.km/u.s):
    '''

    Fit 1 or 2 peaked Gaussian to profile

    TODO: add some sort of diagnostic output for each fit?


    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    5/13/2021   A.A. Kepley     Original Code

    '''

    from astropy.modeling import models, fitting
    from scipy import integrate
    from scipy.stats import f

    # chanwidth
    chanwidth = spectral_axis[1] - spectral_axis[0]

    # construct weight vector
    weights = np.ones(len(stack_profile)) / noisePerChan.value

    # start off by fitting one Gassian
    amp_est = max(stack_profile)
    peak_cut = spectral_axis[stack_profile > amp_est * 0.5]
    fwhm_est = max(peak_cut) - min(peak_cut)
    sigma_est = fwhm_est/2.355
    init_g = models.Gaussian1D(amplitude = amp_est.value, stddev=sigma_est.value)

    fit_g = fitting.LevMarLSQFitter()

    result_g = fit_g(init_g, spectral_axis.value, stack_profile.value, weights=weights)

    # Now fit two Gaussians
    init_g1 = models.Gaussian1D(amplitude = result_g.amplitude/2.0, 
                            stddev = result_g.stddev/2.0, 
                            mean = result_g.mean+result_g.stddev)
    init_g2 = models.Gaussian1D(amplitude = result_g.amplitude/2.0, 
                            stddev = result_g.stddev/2.0, 
                            mean = result_g.mean-result_g.stddev)
    init_g1_g2 = init_g1 + init_g2

    fit_g1_g2 = fitting.LevMarLSQFitter()

    result_g1_g2 = fit_g1_g2(init_g1_g2, spectral_axis.value, stack_profile.value, weights=weights)

    # calculate reduced chi-square for each fit.
    
    chisquare_g, chisquare_r_g = chiSquare(stack_profile.value, result_g(spectral_axis.value),1.0/weights, nparams=3)
    chisquare_g1_g2, chisquare_r_g1_g2 = chiSquare(stack_profile.value, result_g1_g2(spectral_axis.value),1.0/weights, nparams=6)
    
    # calculate f-values by taking ratio of chisquare values
    fval = chisquare_r_g / chisquare_r_g1_g2


    # sf = survival function = 1 - cdf 
    # first parameter in f is the difference in the number of degrees
    # of freedom between the two fits (here 6-3). The second is the
    # number of degrees of freedom in the 2nd (2 gaussian
    # fit). Citation wikipedia article on F statistics and
    # regression. It's consistent with the description in Bevington
    # and Robinson.
    pvalue = f.sf(fval,6-3, (len(stack_profile)-6.0))
    
    # For diagnostics purposes.
    # print(chisquare_r_g, chisquare_r_g1_g2,fval,pvalue)

    # calculate integrated intensity from fits. If the pvalue is
    # small, we reject the null hypothesis that the double-gaussian
    # fits as well as a single gaussian.
    if pvalue < 0.05:
        stack_int = integrate.quad(result_g1_g2, -maxAbsVel.value, maxAbsVel.value)[0] * stack_profile.unit * spectral_axis.unit
        stack_int_err = np.sqrt(fwhm_est/chanwidth) * chanwidth * noisePerChan 
    else:
        stack_int =  integrate.quad(result_g, -maxAbsVel.value, maxAbsVel.value)[0] * stack_profile.unit * spectral_axis.unit
        stack_int_err = np.sqrt(fwhm_est/chanwidth) * chanwidth * noisePerChan

    
    return stack_int, stack_int_err, fwhm_est


def chiSquare(data, fit, error, nparams=1):
    '''
    Calculate reduced Chi squared

    Date        Programmer      Description of Changes
    ---------------------------------------------------------------------
    5/13/2021   A.A. Kepley     Original Code

    '''
    
    dof = len(fit) - nparams
    
    chisquare = np.sum((data-fit)**2/error**2)

    chisquare_reduced = chisquare / dof
    
    return chisquare, chisquare_reduced


def makeBins(galaxy, basemap,  bintype, outDir):
    '''
    Create bins for the data

    basemap: map used to create bins (SpectralCube Projection)

    bintype: type of bins ('intensity','stellarmass', 'radius')
    
    outDir: directory to write output bin image

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added more comments plus moved GCR map 
                                calculate up to other code.
    '''

    if bintype=='intensity' or bintype=='stellarmass':
        #Bin the basemap by brightness
        binnum=int(np.log(np.nanmax(basemap.value)/np.nanmin(basemap.value)))+1
        binedge = np.nanmin(basemap.value)*np.logspace(0, binnum, num=binnum+1, base=np.e) #create bins based on dynamic range of mom0
        bins = np.digitize(basemap.value, binedge) #this will automatically add an extra bin in the end for nan values
        binlabels=[""]+['{0:1.2f}'.format(i)+basemap.unit.to_string() for i in binedge] #need to add units to stellarmass map!!
    elif bintype=='radius':
        binnum=5
        binedge = np.zeros(binnum+1)
        binedge[1:]=np.logspace(-1, np.log10(0.5), num=binnum,base=10)#create bins based on r25
        bins = np.digitize(basemap, binedge) #this will automatically add an extra bin in the end for nan values
        binlabels=[""]+['{0:1.2f}'.format(i)+'R25' for i in binedge]
    else:
        raise Exception ("bintype should either be 'radius' or 'intensity' or 'stellarmass' ")
    # Blank NaN values
    bins[bins==len(binedge)] = 0 
    # make bins map 
    binmap=Projection(bins,wcs=basemap.wcs,header=basemap.header)
    binmap.write(os.path.join(outDir,galaxy['NAME'].upper()+'_binsby'+bintype+'.fits'), overwrite=True)
    return binmap, binedge, binlabels


def makeRadiusBins(galaxy, basemap,  outDir, beam=15.0):
    '''
    Create bins for the data

    basemap: map used to create bins (SpectralCube Projection)
    
    outDir: directory to write output bin image

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added more comments plus moved GCR map 
                                calculate up to other code.
    4/15/2021   A.A. Kepley     Changes bins to be beam width apart in radius.
    '''
    
    minrad = 0.0+beam/2.0
    maxrad = np.max(basemap).value + beam # want to add one bin beyond to capture max.

    binedge = np.arange(minrad, maxrad, beam)  
    binedge = np.insert(binedge,0,0)
    binedge = binedge * basemap.unit
    bins = np.digitize(basemap.value,binedge.value) 

    binlabels = ['{0:1.2f}'.format(i)+' arcsec' for i in binedge]

    # make bins map 
    binmap=Projection(bins,wcs=basemap.wcs,header=basemap.header)
    binmap.write(os.path.join(outDir,galaxy['NAME'].upper()+'_binsbyradius.fits'), overwrite=True)

    return binmap, binedge, binlabels

def plotBins(galaxy, binmap, binedge, binlabels, bintype, outDir):
    '''
    visualize binning

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added comments
    4/29/2021   A.A. Kepley     Simplified plotting code

    '''

    mapvals = np.arange(1.0,len(binedge)+1)    

    plt.subplot(projection=binmap.wcs)
    heatmap = plt.imshow(binmap.value,origin='lower')

    plt.colorbar(heatmap, boundaries=binedge, values=mapvals)

    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')

    plt.savefig(os.path.join(outDir,galaxy['NAME'].upper()+'_binsby'+bintype+'.png'))
    plt.clf()
    plt.close()


def mapSN(galaxy, regridDir, outDir, sncut=3.0):

    '''
    make S/N map using MAD from cube and mask    
    
    galaxy: line from degas_base.fits with galaxy properties

    regridDir: input directory with regridded data

    outDir: output directory

    sncut: keyword parameter specifying S/N cut
    
     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added comments and clarified inputs

    '''

    # read in cube
    cube=SpectralCube.read(os.path.join(regridDir, galaxy['NAME']+'_12CO_regrid.fits'))
    cube=cube.with_spectral_unit(u.km / u.s)

    # read in mask
    mask=SpectralCube.read(os.path.join(regridDir, galaxy['NAME']+'_12CO_mask_regrid.fits'))
    mask=mask.with_spectral_unit(u.km / u.s)

    # calculate noise
    madstd=cube.mad_std(how='cube') #K #raw cube
    chanwidth=np.abs(cube.spectral_axis[0]-cube.spectral_axis[1]) #channel width is same for all channels, km/s
    masksum=mask.sum(axis=0) #map number of unmasked pixels 
    noise=np.sqrt(masksum)*(madstd*chanwidth) #moment0 error map, in K km/s ## TODO: CHECK MATH HERE

    #mask datacube
    masked_cube=cube.with_mask(mask==1.0*u.dimensionless_unscaled) 
    mom0 = masked_cube.moment(order=0)
    snmap=mom0/noise #should be unitless #mom0 is from masked cube
    snmap[snmap==np.inf]=np.nan #convert division by zero to nan

    # write the resulting map to fits 
    snmap.write(os.path.join(outDir, galaxy['NAME'].upper()+'_SNmap.fits'), overwrite=True)
    snmap.quicklook()
    plt.savefig(os.path.join(outDir, galaxy['NAME'].upper()+'_SNmap.png'))
    plt.close()
    plt.clf()

    #get rid of parts of mom0 where S/N < S/N cut
    mom0cut=mom0.copy()
    sn=np.nan_to_num(snmap.value)
    mom0cut[sn<sncut]=np.nan #blank out low sn regions

    #use sigma-clipped mom0 as new 2D mask for the original cube to preserve noise in signal-free channel
    mom0mask=~np.isnan(mom0cut)
    masked_cube=cube.with_mask(mom0mask) #use for stacking later

    return mom0cut, masked_cube, mom0mask

def mapStellar(galaxy, mom0cut, regridDir, outDir):
    '''
    make stellarmass map

    
    galaxy: line from degas_base.fits table with galaxy information

    mom0cut: S/N cut on mom0

    regridDir: input data directory with all regridded data

    outDir: output data directory

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added comments and clarified inputs

    '''

    # open the stellar mass

    ## TODO : NEED TO UPDATE WITH NEW ANCILLAR DATA FROM SARAH
    stellarhdu=fits.open(os.path.join(regridDir,galaxy['NAME']+'_mstar_gauss15_regrid.fits'))[0]

    stellar=stellarhdu.data
    #stellar[starmask==1.0]=np.nan #apply star mask ## AAK: I think I can skip this.
    stellar[np.isnan(mom0cut)]=np.nan #apply SN mask (SN >3)
    w=WCS(stellarhdu.header)
    # stellarhdu.header['BUNIT']='Msun/pc^2' #not the correct unit!!  ## AAK: I think units are okay for the newdata

    stellarmap=Projection(stellar,header=stellarhdu.header,wcs=w) 
    stellarmap.quicklook()

    plt.savefig(os.path.join(outDir,galaxy['NAME']+'_stellarmass.png'))
    plt.clf()
    plt.close()

    return stellarmap

def mapLTIR(galaxy, mom0cut, regridDir, outDir):
    '''
    make LTIR map

    
    galaxy: line from degas_base.fits table with galaxy information

    mom0cut: S/N cut on mom0

    regridDir: input data directory with all regridded data

    outDir: output data directory

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    12/10/2020  A.A. Kepley      Original code based on mapStellar

    '''

    hdu=fits.open(os.path.join(regridDir,galaxy['NAME']+'_LTIR_gauss15_regrid.fits'))[0]

    data=hdu.data

    data[np.isnan(mom0cut)]=np.nan #apply SN mask (SN >3)

    LTIRmap=Projection(data,header=hdu.header,wcs=WCS(hdu.header)) 
    LTIRmap.quicklook()

    plt.savefig(os.path.join(outDir,galaxy['NAME']+'_LTIR.png'))
    plt.clf()
    plt.close()

    return LTIRmap

def mapSFR(galaxy, mom0cut, regridDir, outDir):
    '''
    import sfr map from W4+FUV

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020  A.A. Kepley     Added comments and clarified inputs
    '''
    
    sfrhdu=fits.open(os.path.join(regridDir,galaxy['NAME']+'_sfr_fuvw4_gauss15_regrid.fits'))[0]
    sfr=sfrhdu.data
    sfr[np.isnan(mom0cut)]=np.nan
    w=WCS(sfrhdu.header)
    
    ## TODO CHECK UNITS HERE
    sfrhdu.header['BUNIT']='MSUN/YR/KPC^2'  #might change when getting new files from Sarah
    sfrmap=Projection(sfr,header=sfrhdu.header,wcs=w) 
    sfrmap.quicklook()

    # save plot of map
    plt.savefig(os.path.join(outDir,galaxy['NAME']+'_sfr.png'))
    plt.clf()
    plt.close()
    
    return sfrmap

def mapGCR(galaxy,  basemap):
    '''
    Create map of galactic radii

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/03/2020  A.A. Kepley     Tweak to how input data is handled and 
                                added comments
    '''

    # get the data we need for the galaxy
    ra = galaxy['RA_DEG']
    dec = galaxy['DEC_DEG']
    inc = np.radians(galaxy['INCL_DEG'])
    pa = np.radians(galaxy['POSANG_DEG'])
    r25 = galaxy['R25_DEG'] * 3600.0 # arcsec
    Dmpc = galaxy['DIST_MPC']
    
    # get wcs
    w=basemap.wcs

    # get center
    x0,y0=w.all_world2pix(ra,dec,0,ra_dec_order=True)

    # get coordinates
    y=np.arange(0,np.shape(basemap)[0],1)
    x=np.arange(0,np.shape(basemap)[1],1)

    # create a 2d image of coordinates
    xx,yy=np.meshgrid(x,y)

    # calculate the radius in pixels
    xx_new=x0+(xx-x0)*np.cos(pa)+(yy-y0)*np.sin(pa) 
    yy_new=y0-(xx-x0)*np.sin(pa)+(yy-y0)*np.cos(pa) 
    R=np.sqrt((xx_new-x0)**2/np.cos(inc)**2+(yy_new-y0)**2) #pixels

    # now convert from pixels to actual units
    head=basemap.header
    pxscale=np.radians(np.abs(head['CDELT1'])) #radian/pixel

    R_arcsec=np.degrees(R*pxscale)*3600.0 * u.arcsec# map of GCR in arcsec
    R_kpc=R*pxscale*Dmpc*1000 * u.kpc# map of GCR in kpc
    R_r25=R_arcsec/r25 # map of GCR in units of R25

    Rarcsec_map=Projection(R_arcsec,header=head,wcs=w) 
    Rkpc_map=Projection(R_kpc,header=head,wcs=w) 
    Rr25_map = Projection(R_r25,header=head,wcs=w) 

    # map of galactocentric radius in unit of kpc, arcmin, in r25
    return Rarcsec_map, Rkpc_map, Rr25_map 






    


    
    

