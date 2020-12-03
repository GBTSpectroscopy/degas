from spectral_stack import stacking #need to make sure Erik's code is installed
from spectral_cube import SpectralCube, Projection
from degas.makeMaps import mapSN, mapSFR, mapStellar
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import colors
import aplpy
from astropy.table import Table, Column,vstack,hstack
import re
import glob 
from astropy.wcs import WCS
import os

def makeResultsFITSTable(regridDir, outDir, scriptDir, vtype='mom1', outname='test', release='DR1'):
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
    
    # TODO: CHANGE THIS FROM A GLOB TO A LIST OF DR1 GALAXIES.
    # get list of CO files     
    
    # get list of dr1 galaxies
    degas = Table.read(os.path.join(scriptDir,'degas_base.fits'))
    release = degas[release] == 1

    #    files=glob.glob(os.path.join(regridDir,'*_12CO_regrid.fits')) #get CO cube files for all galaxies

    # go though each file, create a table, and attend it to the list of tables.
    tablelist=[]
    for galaxy in degas[release]:
        full_tab, meta=makeTable(galaxy, vtype, scriptDir, regridDir, outDir)
        tablelist.append(full_tab)

    # stack all the tables together
    table=vstack(tablelist)

    # AAK: Looks like the meta data comes from the last galaxy
    # processed. Is this going to be okay?
    table.meta=meta

    # Write out the table 
    table.write(os.path.join(outDir,outname+'_'+vtype+'.fits'),overwrite=True)

def makeTable(galaxy, vtype, scriptDir, regridDir, outDir):
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

    # get name of galaxy we are processing.
    name = galaxy['NAME']
    
    # For NGC6946, skip 13CO/C18O since that data was never taken.
    if name =='NGC6946':
        linelist=['CO','HCN','HCOp']
    else:
        linelist=['CO','HCN','HCOp','13CO','C18O']

    # for each line stack the data.
    galtabs=[]
    for line in linelist: 
        #get the full stack result for each line
        full_stack, stack_meta=makeStack(galaxy, vtype, line, scriptDir, regridDir, outDir) 

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

    # add the name of the galaxy.
    galtable['galaxy']= name

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


#make stacks for all lines and ancillary data for one galaxy
#output python dictionaries
def makeStack(galaxy, vtype, line, scriptDir, regridDir, outDir):
    '''
    
    make stacks for all lines and ancillary data for one galaxy
    output python dictionaries

    galaxy: degas_base.fits entry for galaxy we are processing
    
    vtype: type of map we are stacking on
    
    line: line we are stacking on
    
    scriptDir: script directory -- if I already have read in
    degas_base.fits do I need to keep passing this along?

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

    '''

    cofile = os.path.join(regridDir,galaxy['NAME']+'_12CO_regrid.fits')
    name = galaxy['NAME']

    #make SN map & sigma-clipped mom0 & 2D masked cube --move to separate script. AAK: only needs to be done once not every time we stack?
    mom0cut, masked_co, sn_mask = mapSN(name, regridDir, outDir)  ## AAK: probably should add in here a keyword to control the mom0cut.
    stellarmap=mapStellar(name,mom0cut, regridDir, outDir)
    sfrmap=mapSFR(name,mom0cut, regridDir, outDir)
    # AAK: add radius map here too?

    # Read in the line file for stacking
    if line=='CO':
        masked_line=masked_co
        sfr=sfrmap
    else:
        if line=='HCN':
            linefile=os.path.join(regridDir,galaxy+'_'+line+'_rebase3_smooth1.3_hanning1_smooth.fits')
        else:
            linefile=os.path.join(regridDir,galaxy+'_'+line+'_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
        sfr=None
        linecube=SpectralCube.read(os.path.join(regridDir,linefile))
        masked_line=linecube.with_mask(sn_mask)

    #import mom1 or  peakvel fits for stacking velocity
    velocity_file=name+'_12CO_'+vtype+'_regrid.fits' #need to change
    vhdu=fits.open(os.path.join(regridDir,velocity_file))
    vhdu[0].header['BUNIT']=masked_co.spectral_axis.unit.to_string()
    velocity=Projection.from_hdu(vhdu)

    ## AAK: NOW I CAN JUST PASS THE RELEVANT TABLE DONE DOWN THE STACKING CODE BELOW???
    table=Table.read(os.path.join(scriptDir,'degas_base.fits'))

    ## filling in the meta information for the stack.
    stack_meta={}
    stack_meta['spectra_axis_unit']=masked_co.spectral_axis.unit.to_string()
    stack_meta['stacked_profile_unit']=masked_co.unit.to_string()
    stack_meta['stacksum_unit']=masked_co.unit.to_string()+masked_co.spectral_axis.unit.to_string()
    stack_meta['sfr_unit']='Msun/yr/kpc^2'
    stack_meta['bin_area']='kpc^2'

    # actually do the stack!

    # AAK: DO I ACTUALLY WANT TO CALCULATE THE STACKING BINS HERE WHEN I STACK?
    
    # Stack on line intensity
    r_intensity=stack(line,masked_line,galaxy, vtype, velocity, mom0cut,'mom0', 'intensity', table, regridDir, outDir, sfrmap=sfr) 

    # stack on stellar mass
    r_stellarmass=stack(line,masked_line,galaxy, vtype,  velocity,  stellarmap, 'stellarmass','stellarmass',table, regridDir, outDir, sfrmap=sfr)

    # stack on radius
    r_radius=stack(line,masked_line, galaxy, vtype,  velocity, mom0cut,'mom0', 'radius', table, regridDir, outDir, sfrmap=sfr)

    # put the individual stack together.
    full_stack=[r_intensity, r_stellarmass, r_radius]

    return full_stack, stack_meta

def makeBins(galaxy, basemap,  bintype, table, outDir):
    if bintype=='intensity' or bintype=='stellarmass':
        #Bin the basemap by brightness
        binnum=int(np.log(np.nanmax(basemap.value)/np.nanmin(basemap.value)))+1
        binedge = np.nanmin(basemap.value)*np.logspace(0, binnum, num=binnum+1, base=np.e) #create bins based on dynamic range of mom0
        bins = np.digitize(basemap.value, binedge) #this will automatically add an extra bin in the end for nan values
        binlabels=[""]+['{0:1.2f}'.format(i)+basemap.unit.to_string() for i in binedge] #need to add units to stellarmass map!!
    elif bintype=='radius':
        R_arcmin, R_kpc, R_r25=mapGCR(galaxy, table, basemap) #can choose from 3 different  Rmaps
        binnum=5
        binedge = np.zeros(binnum+1)
        binedge[1:]=np.logspace(-1, np.log10(0.5), num=binnum,base=10)#create bins based on r25
        bins = np.digitize(R_r25, binedge) #this will automatically add an extra bin in the end for nan values
        binlabels=[""]+['{0:1.2f}'.format(i)+'R25' for i in binedge]
    else:
        raise Exception ("bintype should either be 'radius' or 'intensity' or 'stellarmass' ")
    # Blank NaN values
    bins[bins==len(binedge)] = 0 
    # make bins map 
    binmap=Projection(bins,wcs=basemap.wcs,header=basemap.header)
    binmap.write(os.path.join(outDir,galaxy.upper()+'_binsby'+bintype+'.fits'), overwrite=True)
    return binmap, binedge, binlabels

def plotBins(galaxy, binmap, binedge, binlabels, bintype, outDir):
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
    plt.savefig(os.path.join(outDir,galaxy.upper()+'_binsby'+bintype+'.png'))
    plt.clf()
    plt.close()
    return cmap

def mapGCR(galaxy, table, basemap):
    w=basemap.wcs
    #import basic property, center, inclination, pa
    table=table[table['NAME']==galaxy.upper()]
    ra=table['RA_DEG'][0]
    dec=table['DEC_DEG'][0]
    inc=np.radians(table['INCL_DEG'][0])
    pa=np.radians(table['POSANG_DEG'][0])
    r25=table['R25_DEG'][0]*60 #arcmin
    Dmpc=table['DIST_MPC'][0]

    x0,y0=w.all_world2pix(ra,dec,0,ra_dec_order=True)
    y=np.arange(0,np.shape(basemap)[0],1)
    x=np.arange(0,np.shape(basemap)[1],1)
    xx,yy=np.meshgrid(x,y)
    xx_new=x0+(xx-x0)*np.cos(pa)+(yy-y0)*np.sin(pa) 
    yy_new=y0-(xx-x0)*np.sin(pa)+(yy-y0)*np.cos(pa) 
    R=np.sqrt((xx_new-x0)**2/np.cos(inc)**2+(yy_new-y0)**2) #pixels
    head=basemap.header
    pxscale=np.radians(np.abs(head['CDELT1'])) #radian/pixel
    R_kpc=R*pxscale*Dmpc*1000 # map of GCR in kpc
    R_arcmin=np.degrees(R*pxscale)*60 # map of GCR in arcmin
    R_r25=R_arcmin/r25

    return R_arcmin, R_kpc, R_r25 #map of galactocentric radius in unit of kpc, arcmin, in r25

def stack(line, cube, galaxy, vtype, velocity, basemap, basetype, bintype, table, regridDir, outDir, sfrmap):
    table=table[table['NAME']==galaxy.upper()]
    Dmpc=table['DIST_MPC'][0]
    pix_area=(np.radians(np.abs(cube.header['CDELT1']))*Dmpc*1000)**2  #kpc^2
    nchan=cube.shape[0]

     #base map is the map used to create bins, i.e. mom0 or stellarmass
    binmap, binedge,binlabels=makeBins(galaxy, basemap, bintype, table, outDir)   #create bins by intensity or radius
    cmap=plotBins(galaxy, binmap, binedge, binlabels, bintype, outDir) #plot binmap
    stack, labelvals = stacking.BinByLabel(cube, binmap.value, velocity,
                                           weight_map=None)
    xunit={'intensity':'K km/s','radius':'R25','stellarmass':'Msun/pc^2'}#unit for stellarmass needs to be updated!!
    colors=cmap(np.linspace(0,1,len(stack)+1))
    bin_mean=np.zeros(len(stack))
    stacksum=np.zeros(len(stack))
    stacked_profile=np.zeros((len(stack),len(stack[0]['spectral_axis'])))
    binlabel=np.zeros(len(stack))
    bin_mean=np.zeros(len(stack))
    bin_lower=np.zeros(len(stack))
    bin_upper=np.zeros(len(stack))
    stacknoise=np.zeros(len(stack))
    sfr_mean=np.zeros(len(stack))
    bin_area=np.zeros(len(stack))
    for i in range(len(stack)):
        d=stack[i]
        stacked_profile[i,:]=d['spectrum']
        integral=np.abs(np.trapz(y=d['spectrum'],x=d['spectral_axis']) )#sum intensity under curve, K km/s
        noise_region=np.logical_or(d['spectral_axis'].value < -250, d['spectral_axis'].value > 250)
        channoise=np.nanstd(d['spectrum'][noise_region]) #noise from signal free channels, nned to multiple by sqrt(number of channels) and channel width 
        stacknoise[i]=channoise*np.sqrt(nchan)*np.abs(d['spectral_axis'].value[0]-d['spectral_axis'].value[1])
        bin_lower[i]=binedge[i]
        bin_upper[i]=binedge[i+1]
        bin_mean[i]=(binedge[i]+binedge[i+1])/2
        stacksum[i]=integral.value
        binlabel[i]=d['label']
        if sfrmap is None:
            pass
        else:
            sfr_mean[i]=np.nanmean(sfrmap[binmap==binlabel[i]])
            #calculate the area of each bin in units of kpc^2
            bin_area[i]=len(binmap[binmap==binlabel[i]].flatten())*pix_area

    
    total_stack={}
    if line=='CO': #CO is used to make stacking bins
        total_stack['spectral_axis']=stack[0]['spectral_axis'].value
        total_stack['bin_lower']=bin_lower
        total_stack['bin_upper']=bin_upper
        total_stack['bin_mean']=bin_mean
        total_stack['bin_type']=bintype
        total_stack['bin_unit']=xunit[bintype]
        if sfrmap is None:
            pass
        else:
            total_stack['sfr_mean_w4fuv']=sfr_mean
            total_stack['bin_area']=bin_area
    total_stack[line+'_stack_profile']=stacked_profile
    total_stack[line+'_stack_sum']=stacksum
    total_stack[line+'_stack_noise']=stacknoise
    return total_stack


    


    
    

