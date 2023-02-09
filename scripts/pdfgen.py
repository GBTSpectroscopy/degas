import glob
import numpy as np
from pylatex import Document, Figure, NoEscape, Command, SubFigure, Package, HorizontalSpace, VerticalSpace, Center
from pylatex.base_classes.command import Options
import os
from PyPDF2 import PdfMerger
from astropy.table import Table

imdir = '/lustre/cv/users/akepley/degas/moments_IR6p1/figures/'

def make_multiple(gal_list):
    merger=PdfMerger()
    for gal in gal_list:
        newpdf=make_single(gal)
        merger.append(newpdf)
    merger.write(imdir+'/'+'IR6p1.pdf')
    merger.close()

def make_single(galaxy):

    files = glob.glob(imdir+'/'+galaxy+'')

    geometry_options = {"landscape":False, "lmargin":"1in", "tmargin":"1in"}
    doc = Document(geometry_options=geometry_options, page_numbers=False)
    doc.documentclass = Command('documentclass', options=['11pt'],arguments=['article'])
    with doc.create(Figure(position='b!')) as figure:
        for line in ['HCN', 'HCOp', '13CO', 'C18O']:
            line_mom0_fig = glob.glob(imdir+'/'+galaxy+'_'+line+'_mom0.png')
            line_mom1_fig = glob.glob(imdir+'/'+galaxy+'_'+line+'_mom1.png')
            if len(line_mom0_fig)==0:
                print (galaxy,': No products for ',line,'..skipping...')
            else:
                with figure.create(SubFigure(width=NoEscape(r'0.5\linewidth'),position='b')) as line_mom0:
                    line_mom0.add_image(line_mom0_fig[0], width=NoEscape(r'\linewidth'))
                figure.append(HorizontalSpace("1cm"))
                with figure.create(SubFigure(width=NoEscape(r'0.5\linewidth'),position='b')) as line_mom1:
                    line_mom1.add_image(line_mom1_fig[0], width=NoEscape(r'\linewidth'))
                figure.append(VerticalSpace("1cm"))
                figure.append("\n")
            
    doc.generate_pdf(imdir+'/'+galaxy, clean_tex=False)
    return imdir+'/'+galaxy+'.pdf'


##----make pdfs----##
degas = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))
degas_dr1 = degas[degas['DR1'] == 1]
#degas_dr1 = degas[degas['NAME'] == 'NGC2903']
gal_list = degas_dr1['NAME'] #import list of DEGAS targets

#make_multiple(gal_list, 'native')
make_multiple(gal_list)
