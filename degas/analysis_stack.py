from spectral_stack import stacking #need to make sure Erik's code is installed
from spectral_cube import SpectralCube, Projection
from makeMaps import mapSN, mapSFR, mapStellar
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

scriptDir='/lustre/cv/users/akepley/degas/code/degas/scripts/'
regridDir='/lustre/cv/users/akepley/degas/IR5_regrid/'
outDir='/lustre/cv/users/ysong/degas/testdir_out/'

#main function: assemble a master fitstable from table output of each galaxy
#TO DO: change fits file name once code is officially working
def makeFits(vtype, scriptDir, regridDir, outDir):
    tablelist=[]
    files=glob.glob(os.path.join(regridDir,'*_12CO_regrid.fits')) #get CO cube files for all galaxies
    for f in files:
        full_tab, meta=makeTable(f, vtype, scriptDir, regridDir, outDir)
        tablelist.append(full_tab)
    table=vstack(tablelist)
    table.meta=meta
    table.write(outDir+'test_'+vtype+'_horizontal.fits',overwrite=True)
    return table
    
#make fitstable containing all lines from stacking results for each galaxy
def makeTable(cofile, vtype, scriptDir, regridDir, outDir):
    galaxy=os.path.basename(cofile).split('_')[0]
    if galaxy =='NGC6946':
        linelist=['CO','HCN','HCOp']
    else:
        linelist=['CO','HCN','HCOp','13CO','C18O']
    galtabs=[]
    for line in linelist: #for each line
        full_stack, stack_meta=makeStack(cofile, vtype, line, scriptDir, regridDir, outDir) #get the full stack result for each line
        linetabs=[]
        for i in range(len(full_stack)):
            stack=full_stack[i] #pick out stack for each bin type of this line
            tab=readStack(stack, line) #read stack of this bintype of this line into table format
            linetabs.append(tab)
        full_tab=vstack(linetabs) #vertically merge the stacks of all bintypes for each line
        galtabs.append(full_tab)
    galtable=hstack(galtabs) #merge different lines toegther horizontally
    galtable['galaxy']=galaxy.upper() #add galaxy name
    return galtable, stack_meta

#read stacking results (dictionary output) for each galaxy and each line
def readStack(stack,line):
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
def makeStack(cofile, vtype, line, scriptDir, regridDir, outDir):
    galaxy=os.path.basename(cofile).split('_')[0]
    #make SN map & sigma-clipped mom0 & 2D masked cube --move to separate script
    mom0cut, masked_co, sn_mask = mapSN(galaxy, regridDir, outDir) 
    stellarmap=mapStellar(galaxy,mom0cut, regridDir, outDir)
    sfrmap=mapSFR(galaxy,mom0cut, regridDir, outDir)
    #import mom1 or  peakvel fits for stacking velocity
    velocity_file=galaxy+'_12CO_'+vtype+'_regrid.fits' #need to change
    vhdu=fits.open(os.path.join(regridDir,velocity_file))
    vhdu[0].header['BUNIT']=masked_co.spectral_axis.unit.to_string()
    velocity=Projection.from_hdu(vhdu)

    table=Table.read(os.path.join(scriptDir,'degas_base.fits'))

    stack_meta={}
    stack_meta['spectra_axis_unit']=masked_co.spectral_axis.unit.to_string()
    stack_meta['stacked_profile_unit']=masked_co.unit.to_string()
    stack_meta['stacksum_unit']=masked_co.unit.to_string()+masked_co.spectral_axis.unit.to_string()
    stack_meta['sfr_unit']='Msun/yr/kpc^2'
    stack_meta['bin_area']='kpc^2'
    #all linecubes have been smoothed and regridded at this point
    if line=='CO':
        masked_line=masked_co
        sfr=sfrmap
    elif line=='HCN':
        linefile=os.path.join(regridDir,galaxy+'_'+line+'_rebase3_smooth1.3_hanning1_smooth.fits')
        linecube=SpectralCube.read(os.path.join(regridDir,linefile))
        masked_line=linecube.with_mask(sn_mask)
        sfr=None
    else:
        linefile=os.path.join(regridDir,galaxy+'_'+line+'_rebase3_smooth1.3_hanning1_smooth_regrid.fits')
        linecube=SpectralCube.read(os.path.join(regridDir,linefile))
        masked_line=linecube.with_mask(sn_mask)
        sfr=None


    r_intensity=stack(line,masked_line,galaxy, vtype, velocity, mom0cut,'mom0', 'intensity', table, regridDir, outDir, sfrmap=sfr) 

    r_stellarmass=stack(line,masked_line,galaxy, vtype,  velocity,  stellarmap, 'stellarmass','stellarmass',table, regridDir, outDir, sfrmap=sfr)

    r_radius=stack(line,masked_line, galaxy, vtype,  velocity, mom0cut,'mom0', 'radius', table, regridDir, outDir, sfrmap=sfr)

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


    


    
    

