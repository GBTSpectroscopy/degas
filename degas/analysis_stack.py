from spectral_stack import stacking #need to make sure Erik's code is installed
from spectral_cube import SpectralCube, Projection
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import colors
import aplpy
from astropy.table import Table, Column, vstack, hstack
import re
import glob 
from astropy.wcs import WCS
import os

def makeResultsFITSTable(regridDir, outDir, scriptDir, vtype='mom1', outname='test', release='DR1', sourceList=None):
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
    12/3/2020   A.A. Kepley     Modified to pull list of galaxies from degas base table.

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
            full_tab, meta=makeTable(galaxy, vtype, regridDir, outDir)
            tablelist.append(full_tab)
        

    # stack all the tables together
    table=vstack(tablelist)

    # AAK: Looks like the meta data comes from the last galaxy
    # processed. Is this going to be okay?
    table.meta=meta

    # Write out the table 
    table.write(os.path.join(outDir,outname+'_'+vtype+'.fits'),overwrite=True)

    return table

def makeTable(galaxy, vtype, regridDir, outDir):
    '''

    make fitstable containing all lines from stacking results for each galaxy

    galaxy: data for galaxy we are processing from degas_base.fits

    vtype: velocity type that we are stacking on
    
    scriptDir: script directory -- AAK: if I read in degas_base.fits early,
    so I still need to pass this down.?
    
    regridDir: regrid directory
    
    outDir: output directory

    
    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added comments

    '''

    # For NGC6946, skip 13CO/C18O since we don't have that data.
    if galaxy['NAME'] == 'NGC6946':
        linelist=['CO','HCN','HCOp']
    else:
        linelist=['CO','HCN','HCOp','13CO','C18O']

    # for each line stack the data.
    galtabs=[]
    for line in linelist: 

        #get the full stack result for each line
        full_stack, stack_meta=makeStack(galaxy, vtype, line, regridDir, outDir) 

        # combine the individual line stacks.
        linetabs=[]
        for i in range(len(full_stack)):
            stack=full_stack[i] #pick out stack for each bin type of this line
            tab=readStack(stack, line) #read stack of this bintype of this line into table format
            linetabs.append(tab)
        full_tab=vstack(linetabs) #vertically merge the stacks of all bintypes for each line
        galtabs.append(full_tab)

    # merge different lines together horizontally
    galtable=hstack(galtabs) 

    ## ADD empty 13CO and C18O columns for NGC6946
    ## '13CO_stack_profile','13CO_stack_sum','13CO_stack_noise','C18O_stack_profile','C18O_stack_sum','C18O_stack_noise'
    #if galaxy['NAME'] == 'NGC6946':
    #    for col in ['13CO_stack_sum','13CO_stack_noise','C18O_stack_sum','C18O_stack_noise']:
    #        galtable.add_column(np.zeros(len(galtable['HCN_stack_sum'])), name=col)
    #    for col in ['13CO_stack_profile','C18O_stack_profile']:
    #        galtable.add_column(np.zeros(np.shape(galtable['HCN_stack_profile'])), name=col)
        
            
        
    # add the name of the galaxy.
    galtable['galaxy']= galaxy['NAME']

    # return the table and the stack.
    return galtable, stack_meta


def readStack(stack,line):
    '''
    read stacking results (dictionary output) for each galaxy and each line

    stack: stack result
    line: line name

     Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    10/29/2020  Yiqing Song     Original Code
    12/3/2020   A.A. Kepley     Added more comments

    '''

    cols=[]
    nrow=len(stack[line+'_stack_sum'])
    rows=[]
    for key in stack.keys():
        data=stack[key]
        if key=='spectral_axis':
            col=Column(name=key,dtype=np.float, shape=data.shape)
        elif key == 'bin_type' or key=='bin_unit':
            col=Column(name=key,dtype='U10')
        elif key == line+'_stack_profile':
            col=Column(name=key, dtype=np.float, shape=data.shape[1])
        else:
            col=Column(name=key, dtype=np.float)
        cols.append(col)
    tab=Table()
    tab.add_columns(cols)
    for i in range(nrow):
        row=[]
        for key in stack.keys():
            data=stack[key]
            if key=='spectral_axis' or key=='bin_type' or key=='bin_unit':
                r=data
            elif key==line+'_stack_profile':
                r=data[i,:]
            else:
                r=data[i]
            row.append(r)
        tab.add_row(row)
    return tab


def makeStack(galaxy, vtype, line, regridDir, outDir):
    '''
    
    make stacks for all lines and ancillary data for one galaxy
    output python dictionaries

    galaxy: degas_base.fits entry for galaxy we are processing
    
    vtype: type of map we are stacking on
    
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

    '''

    #make SN map, sigma-clipped mom0, 2D masked cube, and galactocentric radius map.

    #TO DO : move to separate script??  only needs to be done once not every time we 
    # stack? But would need to add file i/o
    mom0cut, masked_co, sn_mask = mapSN(galaxy, regridDir, outDir, sncut=3.0)  
    stellarmap = mapStellar(galaxy,mom0cut, regridDir, outDir)
    sfrmap = mapSFR(galaxy,mom0cut, regridDir, outDir)
    ltirmap = mapLTIR(galaxy,mom0cut,regridDir,outDir)
    R_arcsec, R_kpc, R_r25 = mapGCR(galaxy,mom0cut)


    # Read in the line file for stacking
    if line=='CO':
        masked_line=masked_co
        sfr=sfrmap
        ltir=ltirmap
    else:
        if line=='HCN':
            linefile=os.path.join(regridDir,galaxy['NAME']+'_'+line+'_rebase3_smooth1.3_hanning1_maxnchan_smooth.fits')
        else:
            linefile=os.path.join(regridDir,galaxy['NAME']+'_'+line+'_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
        sfr=None
        ltir=None

        linecube=SpectralCube.read(os.path.join(regridDir,linefile))
        masked_line=linecube.with_mask(sn_mask)

    #import mom1 or  peakvel fits for stacking velocity
    ## TODO: this looks complicated. Can't i just open the velocity file?
    ## TODO: Also don't I only have to to this once for each galaxy?
    velocity_file=galaxy['NAME']+'_12CO_'+vtype+'_regrid.fits' #need to change
    vhdu=fits.open(os.path.join(regridDir,velocity_file))
    #vhdu[0].header['BUNIT']=masked_co.spectral_axis.unit.to_string()
    velocity=Projection.from_hdu(vhdu)
 
    ## filling in the meta information for the stack.
    stack_meta={}
    stack_meta['spectra_axis_unit']=masked_co.spectral_axis.unit.to_string()
    stack_meta['stacked_profile_unit']=masked_co.unit.to_string()
    stack_meta['stacksum_unit']=masked_co.unit.to_string()+masked_co.spectral_axis.unit.to_string()
    stack_meta['sfr_unit']='Msun/yr/kpc^2'
    stack_meta['bin_area']='kpc^2'

    ## create intensity bins
    binmap, binedge, binlabels = makeBins(galaxy, mom0cut, 'intensity', outDir)  
    cmap=plotBins(galaxy, binmap, binedge, binlabels, 'intensity', outDir) 

    # do intensity stack
    r_intensity=stack(line, masked_line, galaxy, velocity, 
                      binmap, binedge, 'intensity','K km/s',
                      sfrmap=sfr, ltirmap=ltir) 

    ## create stellar mass bin
    binmap, binedge, binlabels = makeBins(galaxy, stellarmap, 'stellarmass', outDir)  
    cmap=plotBins(galaxy, binmap, binedge, binlabels, 'stellarmass', outDir) 

    ## do stellar mass stack
    r_stellarmass=stack(line, masked_line, galaxy, velocity,  
                        binmap, binedge,  'stellarmass', 'Msun/pc^2',
                        sfrmap=sfr,ltirmap=ltir)

    ## create radius bins
    binmap, binedge, binlabels = makeRadiusBins(galaxy, R_arcsec, outDir) 
    cmap=plotBins(galaxy, binmap, binedge, binlabels, 'radius', outDir) 

    #breakpoint()

    # stack on radius
    r_radius=stack(line, masked_line, galaxy, velocity,
                   binmap, binedge,  'radius', 'arcsec',
                   sfrmap=sfr,ltirmap=ltir)
                 
    # put the individual stack together.
    full_stack=[r_intensity, r_stellarmass, r_radius]

    return full_stack, stack_meta


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

    R_arcsec=np.degrees(R*pxscale)*3600.0 # map of GCR in arcsec
    R_kpc=R*pxscale*Dmpc*1000 # map of GCR in kpc
    R_r25=R_arcsec/r25 # map of GCR in units of R25

    Rarcsec_map=Projection(R_arcsec,header=head,wcs=w) 
    Rkpc_map=Projection(R_kpc,header=head,wcs=w) 
    Rr25_map = Projection(R_r25,header=head,wcs=w) 

    # map of galactocentric radius in unit of kpc, arcmin, in r25
    return Rarcsec_map, Rkpc_map, Rr25_map 


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

    #breakpoint()
    
    # original code
    #binnum=5
    #binedge = np.zeros(binnum+1)
    #binedge[1:]=np.logspace(-1, np.log10(0.5), num=binnum,base=10)#create bins based on r25
    #bins = np.digitize(basemap, binedge) #this will automatically add an extra bin in the end for nan values
    #binlabels=[""]+['{0:1.2f}'.format(i)+'R25' for i in binedge]
    # Blank NaN values
    #bins[bins==len(binedge)] = 0 

    minrad = 0.0+beam/2.0
    maxrad = np.max(basemap).value

    binedge = np.arange(minrad, maxrad, beam)
    binedge = np.insert(binedge,0,0)
    bins = np.digitize(basemap,binedge) 

    binlabels = [""] + ['{0:1.2f}'.format(i)+' arcsec' for i in binedge]
    # Blank NaN values
    bins[bins==len(binedge)] = 0 

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

    '''
    

    #visualize binning 
    cmap= plt.cm.viridis
    cmaplist=[cmap(i) for i in range(cmap.N)]
    cmap=colors.LinearSegmentedColormap.from_list('Custom',cmaplist, cmap.N)
    norm=colors.BoundaryNorm(range(len(binedge)+1),cmap.N)

    heatmap=plt.imshow(binmap.value, cmap=cmap, norm=norm)
    plt.close()

    fig=aplpy.FITSFigure(binmap.hdu) #show coordinates
    fig.show_colorscale(cmap=cmap)
    cbar=plt.colorbar(heatmap)

    cbar.ax.set_yticklabels(binlabels) 
    cbar.ax.text(1.5,1.0/(2*len(binedge)),'NaN') 
    plt.savefig(os.path.join(outDir,galaxy['NAME'].upper()+'_binsby'+bintype+'.png'))
    plt.clf()
    plt.close()

    return cmap

def stack(line, cube, galaxy, velocity, 
          binmap, binedge, bintype, binunit,
          sfrmap=None, ltirmap=None,
          maxAbsVel = 250.0):

    '''
    Actually do the stacking.

    line: name of line we are stacking (string)

    cube: data we are stacking (SpectralCube)

    galaxy: line from degas_base.fits table with galaxy information.
    
    velocity: velocity map to stack data using

    binmap: bin map

    binedge: bin edges
    
    binlabel: bin labels
    
    bintype: type of bin
    
    binunit: units on bin [CAN I GET THIS FROM BINMAP?]

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

    '''

    # get the relevant info on the galaxy and cube
    Dmpc=galaxy['DIST_MPC']
    pix_area=(np.radians(np.abs(cube.header['CDELT1']))*Dmpc*1000)**2  #kpc^2
    nchan=cube.shape[0]

    cube_nonan = cube.with_fill_value(0)

    # do the stack -- making sure the units work out right.
    stack, labelvals = stacking.BinByLabel(cube_nonan,
                                           binmap.value, velocity,
                                           weight_map=None)
    
    # set up the results output
    bin_mean = np.zeros(len(stack))
    stacksum = np.zeros(len(stack))
    stacked_profile = np.zeros((len(stack),len(stack[0]['spectral_axis'])))
    binlabel = np.zeros(len(stack))
    bin_mean = np.zeros(len(stack))
    bin_lower = np.zeros(len(stack))
    bin_upper = np.zeros(len(stack))
    stacknoise = np.zeros(len(stack))
    sfr_mean = np.zeros(len(stack))
    sfr_total = np.zeros(len(stack))
    ltir_mean = np.zeros(len(stack))
    ltir_total = np.zeros(len(stack))
    bin_area = np.zeros(len(stack))

    # now put together the stacking results
    for i in range(len(stack)):

        d = stack[i]
        
        # profile
        stacked_profile[i,:] = d['spectrum']

        # integrated profile
        integral = np.abs(np.trapz(y=d['spectrum'], x=d['spectral_axis']) )#sum intensity under curve, K km/s
        stacksum[i] = integral.value

        # calculating the noise
        noise_region = np.logical_or(d['spectral_axis'].value < -abs(maxAbsVel),
                                     d['spectral_axis'].value > abs(maxAbsVel))
        channoise = np.nanstd(d['spectrum'][noise_region]) #noise from signal free channels, nned to multiple by sqrt(number of channels) and channel width 
        stacknoise[i] = channoise*np.sqrt(nchan)*np.abs(d['spectral_axis'].value[0]-d['spectral_axis'].value[1])

        # putting in the bins
        bin_lower[i] = binedge[i]
        bin_upper[i] = binedge[i+1]
        bin_mean[i] = (binedge[i]+binedge[i+1])/2
        #calculate the area of each bin in units of kpc^2    
        binlabel[i] = d['label']
        bin_area[i] = len(binmap[binmap==binlabel[i]].flatten())*pix_area    
        
        # calculating the mean SFR (as necessary)
        if sfrmap is not None:
            sfr_mean[i]=np.nanmean(sfrmap[binmap==binlabel[i]])

        # calculate the mean LTIR (as necessary)
        if ltirmap is not None:
            ltir_mean[i] = np.nanmean(ltirmap[binmap==binlabel[i]])
    
    ## AAK: Not quite sure I understand what's happening here.
    total_stack={}
    if line == 'CO': #CO is used to make stacking bins
        total_stack['spectral_axis'] = stack[0]['spectral_axis'].value
        total_stack['bin_lower'] = bin_lower
        total_stack['bin_upper'] = bin_upper
        total_stack['bin_mean'] = bin_mean
        total_stack['bin_type'] = bintype
        total_stack['bin_unit'] = binunit ### TODO: CAN I GET THIS FROM THE BINMAP DIRECTLY?

        if sfrmap is not None:
            total_stack['sfr_mean_w4fuv']=sfr_mean
            total_stack['sfr_total_w4fuv']=sfr_mean * bin_area

        if ltirmap is not None:
            total_stack['ltir_mean'] = ltir_mean
            total_stack['ltir_total'] = ltir_mean * bin_area

    total_stack[line+'_stack_profile']=stacked_profile
    total_stack[line+'_stack_sum']=stacksum
    total_stack[line+'_stack_noise']=stacknoise

    return total_stack


    


    
    

