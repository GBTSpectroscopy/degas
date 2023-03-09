#modified based on PHANGS MegaTable code
#generates hexgonal grid for sampling individual measurement on mom0 map for A SINGLE galaxy and A SINGLE transition 
import sys
sys.path.append('/lustre/cv/users/ysong/degas/phangs_megatable')
from core import VoronoiTessTable
from matplotlib import pyplot as plt
from matplotlib.patches import RegularPolygon
from astropy.io import fits
import os
import glob
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.table import Table, vstack
from astropy.io import ascii as at

#import info table 
degas=Table.read('/lustre/cv/users/ysong/degas/degas/scripts/degas_base.fits')
basedir='/lustre/cv/users/ysong/degas/'
productdir=os.path.join(basedir,'moments_IR6p1/')
IRdir=os.path.join(basedir,'IR6p1_regrid/')
outputdir=os.path.join(basedir,'test_hex_output/')
datadir=os.path.join(basedir,'GScompare/')

def GScompare(tablefile, dataDir=datadir, outDir=outputdir):
    f1, axs1=plt.subplots(1,2, figsize=(12,7), sharey=True)
    plt.subplots_adjust(top=0.7, wspace=0.4)
    f2, axs2=plt.subplots(1,2, figsize=(12,7), sharey=True)
    plt.subplots_adjust(top=0.7, wspace=0.4)
    filelist=glob.glob(dataDir+'*.ecsv')
    colors=plt.cm.tab20b(np.linspace(0,1,len(filelist))) 
    for (i, c) in zip(range(len(filelist)), colors): #plot comparison datasets, with unique color coding
        table=at.read(filelist[i])
        if 'LIR' in table.colnames and 'LHCN' in table.colnames: 
            axs1[0].loglog(table['LHCN'],table['LIR'], ls='', marker='o', c=c, alpha=0.6, label=table.meta['source'])
            axs1[1].loglog(table['LHCN'],table['LIR'], ls='', marker='o', c='gray', alpha=0.6)
            axs2[0].loglog(table['LHCN'],table['LIR']/table['LHCN'], ls='', marker='o', c=c, alpha=0.6, label=table.meta['source'])
            axs2[1].loglog(table['LHCN'],table['LIR']/table['LHCN'], ls='', marker='o', c='gray', alpha=0.6)

    mastertable=at.read(tablefile)
    galaxylist=list(set(mastertable['GALAXY'])) #actual list of galaxies meausred 
    axs1[0].loglog(mastertable['L_HCN'], mastertable['L_TIR'],marker='o',mfc='gray', mec='black', ls=' ', alpha=0.7, label='DEGAS')
    axs2[0].loglog(mastertable['L_HCN'], mastertable['L_TIR']/mastertable['L_HCN'],marker='o',mfc='gray', mec='black', alpha=0.7, ls=' ', label='DEGAS')
    # setup plot style as in analysis_plot for DR1 targets
    markers = ['o','v','^','s','>','D'] # 6 items
    colors = ['royalblue','forestgreen','darkorange','royalblue','crimson','rebeccapurple','darkcyan','darkmagenta']
    release = degas['DR1'] == 1
    nrelease = np.sum(release)
    markerlist = np.tile(markers,int(np.ceil(nrelease/len(markers))))
    markerlist = markerlist[0:nrelease]

    colorlist = np.tile(colors,int(np.ceil(nrelease/len(colors))))
    colorlist = colorlist[0:nrelease]
    # for each galaxy that is in DR1 that is supposed to be measured 
    for (galaxy, color, marker)  in zip(degas[release], colorlist, markerlist):
        if galaxy['NAME'] in galaxylist: # in case dr1 galaxy was not actually measured 
            hextab=mastertable[mastertable['GALAXY']==galaxy['NAME']]
            axs1[1].loglog(hextab['L_HCN'], hextab['L_TIR'],marker=marker,ls=' ',mfc=color, mec='black', alpha=0.7, label=galaxy['NAME'])
            axs2[1].loglog(hextab['L_HCN'], hextab['L_TIR']/hextab['L_HCN'],marker=marker,ls=' ',mfc=color, mec='black', alpha=0.7, label=galaxy['NAME'])
    xrange=np.linspace(0,1e10,1000)
    axs1[0].loglog(xrange, xrange*90, ls='--', c='black')
    axs1[0].loglog(xrange, xrange*900, ls='solid', c='black')
    axs1[0].loglog(xrange, xrange*9000, ls='-.', c='black')
    axs1[1].loglog(xrange, xrange*90, ls='--', c='black')
    axs1[1].loglog(xrange, xrange*900, ls='solid', c='black')
    axs1[1].loglog(xrange, xrange*9000, ls='-.', c='black')
    axs2[0].loglog(xrange, np.ones(len(xrange))*90, ls='--', c='black')
    axs2[0].loglog(xrange, np.ones(len(xrange))*900, ls='solid', c='black')
    axs2[0].loglog(xrange, np.ones(len(xrange))*9000, ls='-.', c='black')
    axs2[1].loglog(xrange, np.ones(len(xrange))*90, ls='--', c='black')
    axs2[1].loglog(xrange, np.ones(len(xrange))*900, ls='solid', c='black')
    axs2[1].loglog(xrange, np.ones(len(xrange))*9000, ls='-.', c='black')
    axs1[0].legend(ncol=3, columnspacing=0.6, loc='center', bbox_to_anchor=(0.5, 1.2))
    axs1[1].legend(ncol=4, columnspacing=0.6, loc='center',bbox_to_anchor=(0.5, 1.2))
    axs1[0].set_xlim(0, 1e10)
    axs1[1].set_xlim(0, 1e10)
    axs1[0].set_xlabel(r'log $L_{\rm HCN}$ (K km/s pc$^2) $') 
    axs1[0].set_ylabel(r'log $L_{\rm TIR}/L_{\rm HCN}$ ($L_\odot$/(K km/s pc$^2$))')
    axs1[1].set_xlabel(r'log $L_{\rm HCN}$ (K km/s pc$^2) $') 
    axs1[1].set_ylabel(r'log $L_{\rm TIR}$ ($L_\odot$)')
    axs2[0].legend(ncol=3, columnspacing=0.6, loc='center',bbox_to_anchor=(0.5,1.2))
    axs2[1].legend(ncol=4, columnspacing=0.6, loc='center',bbox_to_anchor=(0.5,1.2))
    axs2[0].set_xlabel(r'log $L_{\rm HCN}$ (K km/s pc$^2) $') 
    axs2[0].set_ylabel(r'log $L_{\rm TIR}/L_{\rm HCN}$ ($L_\odot$/(K km/s pc$^2$))')
    axs2[1].set_xlabel(r'log $L_{\rm HCN}$ (K km/s pc$^2) $') 
    axs2[1].set_ylabel(r'log $L_{\rm TIR}$ ($L_\odot$)')
    axs2[0].set_xlim(0, 1e10)
    axs2[1].set_xlim(0, 1e10)
    f1.savefig(outDir+'DEGAS_GS_LIRvLHCN.png')
    f2.savefig(outDir+'DEGAS_GS_LIRLHCNvLHCN.png')
    

def measure_hex_all(inDir=productdir, outDir=outputdir,  line=['HCN'], plot_grid=True):
    #get a list of galaxies that have data --> may consider changing to DR1 galaxy list  from degas base fits
    gallist=glob.glob(inDir+'*HCN*_mom0_masked.fits') 
    hextablist=[]
    for galfile in gallist:
        galaxy=os.path.basename(galfile).split('_')[0]
        hextab=measure_hex_single(galaxy, inDir=productdir, outDir=outputdir, line=line, plot_grid=plot_grid)
        
        newtab=at.read(outDir+galaxy+'_hex.ecsv')
        hextablist.append(newtab)
    hextab_all=vstack(hextablist)
    hextab_all.write(outDir+'DEGAS_hex_megatable.ecsv',overwrite=True) 
    return hextab_all

def measure_hex_single(galaxy, inDir=productdir, outDir=outputdir, line=['HCN'], plot_grid=True):
    '''
    use phangs megatable function to generate hex grid for a galaxy based on its nominal coordinates.
    galaxy: name of galaxy as appeared in degas_base.fits table
    inDir: where line data products live
    outDir: where to output final measurements, and plots

    line: str list, which line(s) to measure, default is HCN. If empty, then all lines are measured
    plot_grid: whether to plot the hex grid, default True
    '''
    #get basic info of galay
    idx=np.where(degas['NAME']==galaxy)[0][0]
    RA=degas['RA_DEG'][idx]
    DEC=degas['DEC_DEG'][idx]
    Dpc=degas['DIST_MPC'][idx]*1e6*u.pc
    
    #get data of galaxy
    TIRfile=glob.glob(IRdir+galaxy+'_LTIR_gauss15_regrid.fits')[0]
    TIRhdu=fits.open(TIRfile)[0]
    #beam=TIRhdu.header['BMAJ'] #deg
    beam=15.0/3600
    beam_area_pc2=((beam*3600/206265.0)*Dpc)**2*np.pi/4/np.log(2) #conversion factor

    #initialize hex grid 
    hextab=initialize_hex_grid(RA, DEC, beam, FOV=120/3600)
    hextab['GALAXY']=galaxy
    hextab['D_pc']=Dpc
    
    #measure IR emission
    hextab.resample_image(image=TIRhdu, unit='header', colname='sigL_TIR',fill_outside=np.nan)
    hextab['L_TIR']=hextab['sigL_TIR']*beam_area_pc2
    if line is None:
        line=['HCN', '13CO', 'C18O', 'HCOp']

    for l in line:
        mom0file=glob.glob(inDir+galaxy+'_'+l+'*_mom0_masked.fits')[0]
        mom0hdu=fits.open(mom0file)[0]
        hextab.resample_image(image=mom0hdu, unit='header', colname='S_'+l, fill_outside=np.nan)
        hextab['L_'+l]=hextab['S_'+l]*beam_area_pc2
        if plot_grid:
            fig=plot_hex_grid(hextab, mom0hdu, beam) 
            fig.savefig(outDir+galaxy+'_'+l+'_hexgrid.png')

    #trim rows with NaN values for HCN
    hextab=hextab[~np.isnan(hextab['S_HCN'])]
    hextab.write(outDir+galaxy+'_hex.ecsv', overwrite=True)

    return hextab

def initialize_hex_grid(RA_center, DEC_center, beam, FOV=120/3600):
    '''
    generate hexagon grid based on a center coordinate and cell size
    '''
    hextab=VoronoiTessTable(center_ra=RA_center, center_dec=DEC_center, fov_radius=FOV,seed_spacing=beam/2*np.sqrt(3),tile_shape='hexagon')
    return hextab

def plot_hex_grid(hextab, hdu, beam):
    '''
    plot grid over data 
    '''
    w=WCS(hdu)
    pxscale=hdu.header['CDELT2']
    f=plt.figure(figsize=(6,6)) #TO DO: add flag to turn on off the plotting
    #ax=f.add_subplot(111, projection=w) #<----try this 
    ax=f.add_subplot()
    ax.imshow(hdu.data, origin='lower')
    #ax.set_autoscale_on(True) #TO DO: add title 

    for x,y in zip(hextab['RA'].value,hextab['DEC'].value):
        #hexagon=RegularPolygon((x,y),numVertices=6, radius=beam/2, facecolor='none',edgecolor='b',transform=ax.get_transform(w))
        xp,yp=w.all_world2pix(x,y,0, ra_dec_order=True)
        hexagon=RegularPolygon((xp,yp),numVertices=6, radius=beam/2/pxscale, facecolor='none',edgecolor='b')
        ax.add_patch(hexagon)
        #ax.scatter(x,y,c='red',s=5,alpha=0.6,transform=ax.get_transform('fk5'))
    return f
    
