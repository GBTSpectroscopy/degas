from spectral_stack import stacking #need to make sure Erik's code is installed
from spectral_cube import SpectralCube, Projection
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import colors
import aplpy
from astropy.table import Table, Column,vstack,hstack
import re
from astropy.wcs import WCS

#main function: assemble a master fitstable from table output of each galaxy
#TO DO: change fits file name once code is officially working
def makeFits(files, vtype, tabdir='./'):
    tablelist=[]
    for f in files:
        full_tab, meta=makeTable(f, vtype)
        tablelist.append(full_tab)
    table=vstack(tablelist)
    table.meta=meta
    table.write(tabdir+'test_'+vtype+'_horizontal.fits',overwrite=True)
    return table
    
#make fitstable containing all lines from stacking results for each galaxy
def makeTable(cubefile, vtype):
    galaxy=cubefile.split('_')[0]
    linelist=['CO','HCN','HCOp','13CO','C18O']
    galtabs=[]
    for line in linelist: #for each line
        full_stack, stack_meta=makeStack(cubefile, vtype, line) #get the full stack result for each line
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


#make stacks for all lines and ancillary data
#output python dictionaries
def makeStack(co_file, vtype, line, datadir='./'):
    galaxy=datadir+co_file.split('_')[0]
    co_cube = SpectralCube.read(co_file)
    co_cube = co_cube.with_spectral_unit(u.km / u.s)
    mask=SpectralCube.read(co_file.replace('.fits','_mask.fits'))
    #make SN map & sigma-clipped mom0 & 2D masked cube
    sn_map,mom0cut, masked_co, sn_mask = mapSN(galaxy, co_cube,mask) 
    stellarmap=mapStellar(galaxy,mom0cut)
    sfrmap=mapSFR(galaxy,mom0cut)
      
    #all linecubes have been smoothed and regridded at this point
    if line=='CO':
        linecube=co_cube
        masked_line=masked_co
        sfr=sfrmap
    else:
        linefile=galaxy+'_'+line+'_rebase3_smooth1.3_hanning1_smooth_regrid.fits'
        linecube=SpectralCube.read(linefile)
        linecube=linecube.with_spectral_unit(u.km/u.s)
        masked_line=linecube.with_mask(sn_mask)
        sfr=None

    r_intensity=stack(line,masked_line,galaxy, vtype, mom0cut,'mom0', 'intensity',sfr) 

    r_stellarmass=stack(line,masked_line,galaxy, vtype, stellarmap,'stellarmass','stellarmass',sfr)

    r_radius=stack(line,masked_line, galaxy, vtype, mom0cut,'mom0','radius',sfr)

    full_stack=[r_intensity, r_stellarmass, r_radius]

    stack_meta={}
    stack_meta['spectra_axis_unit']=linecube.spectral_axis.unit.to_string()
    stack_meta['stacked_profile_unit']=linecube.unit.to_string()
    stack_meta['stacksum_unit']=linecube.unit.to_string()+linecube.spectral_axis.unit.to_string()
    stack_meta['sfr_unit']='Msun/yr/kpc^2'
    stack_meta['bin_area']='kpc^2'
   
    #plotSpectra(galaxy.upper(),cube,mom0cut) #plot some central spectra
    return full_stack, stack_meta

#make stellarmass map
#TO DO: update once get new ancillary data from Sarah
def mapStellar(galaxy,mom0cut,plotdir='./',stellardatadir='./'):
    stellarhdu=fits.open(stellardatadir+galaxy+'_w1_stellarmass_regrid.fits')[0]
    starmask=fits.getdata(stellardatadir+galaxy+'_w1_gauss15_stars_regrid.fits')
    stellar=stellarhdu.data
    stellar[starmask==1.0]=np.nan
    stellar[np.isnan(mom0cut)]=np.nan
    w=WCS(stellarhdu.header)
    stellarhdu.header['BUNIT']='Msun/pc^2' #not the correct unit!!
    stellarmap=Projection(stellar,header=stellarhdu.header,wcs=w) 
    stellarmap.quicklook()
    plt.savefig(plotdir+galaxy+'_stellarmass.png')
    plt.clf()
    plt.close()
    return stellarmap

#import sfr map from W4+FUV
def mapSFR(galaxy,mom0cut,plotdir='./',sfrdatadir='./'):
    sfrhdu=fits.open(sfrdatadir+galaxy+'_w4fuv_sfr_regrid.fits')[0]
    starmask_fuv=fits.getdata(sfrdatadir+galaxy+'_fuv_gauss15_stars_regrid.fits')
    starmask_nuv=fits.getdata(sfrdatadir+galaxy+'_nuv_gauss15_stars_regrid.fits')
    sfr=sfrhdu.data
    sfr[starmask_fuv==1.0]=np.nan #mask out stars from fuv and nuv
    sfr[starmask_nuv==1.0]=np.nan
    sfr[np.isnan(mom0cut)]=np.nan
    w=WCS(sfrhdu.header)
    sfrhdu.header['BUNIT']='MSUN/YR/KPC^2'  
    sfrmap=Projection(sfr,header=sfrhdu.header,wcs=w) 
    sfrmap.quicklook()
    plt.savefig(plotdir+galaxy+'_sfr.png')
    plt.clf()
    plt.close()
    return sfrmap


#make S/N map using MAD from cube and mask
def mapSN(galaxy, cube, mask, datadir='./', plotdir='./'):
    madstd=cube.mad_std(how='cube') #K #raw cube
    chanwidth=np.abs(cube.spectral_axis[0]-cube.spectral_axis[1]) #channel width is same for all channels, km/s
    masksum=mask.sum(axis=0) #map number of unmasked pixels 
    noise=np.sqrt(masksum)*(madstd*chanwidth) #moment0 error map, in K km/s
    #mask datacube
    masked_cube=cube.with_mask(mask==1.0*u.dimensionless_unscaled) 
    mom0 = masked_cube.moment(order=0)
    snmap=mom0/noise #should be unitless #mom0 is from masked cube
    snmap[snmap==np.inf]=np.nan #convert division by zero to nan
    snmap.write(datadir+galaxy.upper()+'_SNmap.fits', overwrite=True)
    snmap.quicklook()
    plt.savefig(plotdir+galaxy.upper()+'_SNmap.png')
    plt.close()
    plt.clf()
    #get rid of parts of mom0 where S/N<3
    mom0cut=mom0.copy()
    sn=np.nan_to_num(snmap.value)
    mom0cut[sn<3.0]=np.nan #blank out low sn regions
    #use sigma-clipped mom0 as new 2D mask for the original cube to preserve noise in signal-free channel
    mom0mask=~np.isnan(mom0cut)
    masked_cube=cube.with_mask(mom0mask) #use for stacking later
    return snmap, mom0cut, masked_cube, mom0mask

def makeBins(galaxy, basemap,  bintype):
    if bintype=='intensity' or bintype=='stellarmass':
        #Bin the basemap by brightness
        binnum=int(np.log(np.nanmax(basemap.value)/np.nanmin(basemap.value)))+1
        binedge = np.nanmin(basemap.value)*np.logspace(0, binnum, num=binnum+1, base=np.e) #create bins based on dynamic range of mom0
        bins = np.digitize(basemap.value, binedge) #this will automatically add an extra bin in the end for nan values
        binlabels=[""]+['{0:1.2f}'.format(i)+basemap.unit.to_string() for i in binedge] #need to add units to stellarmass map!!
    elif bintype=='radius':
        R_arcmin, R_kpc, R_r25=mapGCR(galaxy,basemap) #can choose from 3 different  Rmaps
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
    binmap.write(datadir+galaxy.upper()+'_binsby'+bintype+'.fits', overwrite=True)
    return binmap, binedge, binlabels

def plotBins(galaxy, binmap, binedge, binlabels, bintype,plotdir='./'):
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
    plt.savefig(plotdir+galaxy.upper()+'_binsby'+bintype+'.png')
    plt.clf()
    plt.close()
    return cmap

def mapGCR(galaxy,basemap,degasdir='./'):
    w=basemap.wcs
    #import basic property, center, inclination, pa
    table=Table.read(degasdir+'degas_base.fits')
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

def stack(line,cube, galaxy, vtype, basemap, basetype, bintype, sfrmap=None, weightmap=None,degasdir='./',datadir='./'):
    #create bins by intensity or radius
    #base map is the map used to create bins, i.e. mom0 or stellarmass
    #import mom1 or  peakvel fits for stacking velocity
    table=Table.read(degasdir+'degas_base.fits')
    table=table[table['NAME']==galaxy.upper()]
    Dmpc=table['DIST_MPC'][0]
    pix_area=(np.radians(np.abs(cube.header['CDELT1']))*Dmpc*1000)**2  #kpc^2
    nchan=cube.shape[0]
    velocity_file=galaxy+'_'+vtype+'.fits'
    vhdu=fits.open(datadir+velocity_file)
    velocity=Projection.from_hdu(vhdu)
    vtype=velocity_file.split('_')[1].split('.')[0] #mom1 or peakvel
    binmap, binedge,binlabels=makeBins(galaxy, basemap, bintype)
    cmap=plotBins(galaxy, binmap, binedge, binlabels, bintype) #plot binmap
    stack, labelvals = stacking.BinByLabel(cube, binmap.value, velocity,
                                           weight_map=weightmap)
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


    


    
    

