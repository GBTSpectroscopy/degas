import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
import glob
import os
import pandas as pd
from astropy import units as u
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import ascii as at

indir='/lustre/cv/users/ysong/degas/test_product/'
outdir='/lustre/cv/users/ysong/degas/test_product/figures/'
colors={'13CO':'rainbow','C18O':'rainbow','HCN':'rainbow','HCOp':'rainbow'}
sample=at.read('/lustre/cv/users/ysong/degas/degas/degas/data/dense_survey.cat')

def make_all_res(res, gal_list=None, line_list=None):
    if gal_list is None:
        gal_list=sample['NAME'].data #import list of DEGAS targets
    if line_list is None:
        line_list=['13CO', 'C18O', 'HCN', 'HCOp']

    for gal in gal_list:
        for line in line_list:
            print ('products for', gal, '...', line, 'generating...')
            make_single_line(gal, line, res=res, moments=[0,1])


def make_single_line(galaxy, line,  res, moments=[0,1]):
    # galaxy: name of galaxy, capitalized
    # line: '13CO', 'C18O', 'HCN', 'HCOp'
    # res: 'native' or 'smoothed'
    # moments: [0], [1], [0,1]

    for m in moments:
        momfits=glob.glob(indir+res+'/'+galaxy+'_'+line+'*_mom'+str(m)+'.fits')        
        emomfits=glob.glob(indir+res+'/'+galaxy+'_'+line+'*_emom'+str(m)+'.fits')

        if len(momfits) == 0:
            print ('No products. Moving on....')
        else:
            mom=fits.open(momfits[0])[0]
            emom=fits.open(emomfits[0])[0]

            fig=plt.figure(figsize=(8.5,4.5))
            f1=aplpy.FITSFigure(mom, figure=fig, subplot=[0.12, 0.1, 0.4, 0.8])
            f1.show_colorscale(cmap='rainbow',vmax=0.8*np.nanmax(mom.data), interpolation='bilinear')
            f1.add_beam()
            f1.beam.set_color('red')
            f1.add_colorbar('top')
            f1.colorbar.set_axis_label_text(mom.header['BUNIT'])
            f1.add_label(0.5, 0.8, line+'_mom'+str(m), relative=True)
            f1.add_label(0.5, 0.1, galaxy, relative=True,size='large',weight='medium')

            f2=aplpy.FITSFigure(emom, figure=fig, subplot=[0.57, 0.1, 0.4, 0.8])
            f2.show_colorscale(cmap='cubehelix',vmax=np.nanmax(emom.data), interpolation='bilinear')
            #f2.add_beam()
            #f2.beam.set_color('white')
            f2.add_colorbar('top')
            f2.colorbar.set_axis_label_text(emom.header['BUNIT'])
            f2.add_label(0.5, 0.8, line+'_emom'+str(m), relative=True)
            f2.add_label(0.5, 0.1, galaxy, relative=True, size='large', weight='medium')

            f2.axis_labels.hide_y()
            f2.tick_labels.hide_y()

            fig.savefig(outdir+res+'/'+galaxy+'_'+line+'_'+res+'_mom'+str(m)+'.png')
