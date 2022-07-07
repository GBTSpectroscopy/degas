import glob
import numpy as np
from pylatex import Document, Figure, NoEscape, Command, SubFigure, Package, HorizontalSpace, VerticalSpace, Center
from pylatex.base_classes.command import Options
import os
from PyPDF2 import PdfFileMerger

imdir='/lustre/cv/users/ysong/degas/test_product/figures/'

def make_multiple(gal_list, res):
    merger=PdfFileMerger()
    for gal in gal_list:
        newpdf=make_single(gal, res)
        merger.append(newpdf)
    merger.write(imdir+res+'/'+'IR6p1_'+res+'.pdf')
    merger.close()

def make_single(galaxy, res):
    files=glob.glob(imdir+res+'/'+galaxy+'')
    geometry_options={"landscape":False, "lmargin":"1in", "tmargin":"1in"}
    doc=Document(geometry_options=geometry_options, page_numbers=False)
    doc.documentclass = Command('documentclass', options=['11pt'],arguments=['article'])
    with doc.create(Figure(position='b!')) as figure:
        for line in ['HCN', 'HCOp', '13CO', 'C18O']:
            line_mom0_fig=glob.glob(imdir+res+'/'+galaxy+'_'+line+'_'+res+'_mom0.png')
            line_mom1_fig=glob.glob(imdir+res+'/'+galaxy+'_'+line+'_'+res+'_mom1.png')
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
            
    doc.generate_pdf(imdir+res+'/'+galaxy+'_'+res, clean_tex=False)
    return imdir+res+'/'+galaxy+'_'+res+'.pdf'


##----make pdfs----##
allfiles=glob.glob(imdir+'native'+'/'+'*.png')
gal_list=[]
for file in allfiles:
    name=os.path.basename(file).split('_')[0]
    if name not in gal_list:
        gal_list.append(name)
make_multiple(gal_list, 'native')
make_multiple(gal_list, 'smoothed')
