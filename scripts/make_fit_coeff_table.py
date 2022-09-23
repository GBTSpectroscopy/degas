# Purpose: make latex table showing fit coefficients for paper
#
# Date          Programmer              Description of Changes
#----------------------------------------------------------------------
# 9/22/2022     A.A. Kepley             Original Code

import os
from astropy.table import Table
import numpy as np
from importlib import reload
from astropy import units as u
import math

## TODO: REMEMBER THAT RIGHT NOW I'M EXCLUDING NGC4414 FROM THE FITS I
## AM SHOWING.

#----------------------------------------------------------------------
#                       read in fit results
#----------------------------------------------------------------------

analysis_dir = os.environ['ANALYSISDIR']
stack_fits= Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1','stack_fits.fits'))
stack_fits_noNGC4414 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1','stack_fits_noNGC4414.fits'))
stack_fits_spatialR21 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21','stack_fits.fits'))
stack_fits_spatialR21 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21','stack_fits_noNGC4414.fits'))

#----------------------------------------------------------------------
#                       Put together table header
#----------------------------------------------------------------------

outfile = os.path.join(os.environ['ANALYSISDIR'],'tables','sample_table.tex')

header = '''\\begin{deluxetable}{}
\\tablewidth{0pt}
\\tabletypesize{\scriptsize}
\\tablecaption{Fit Coefficients \label{tab:fit_coeff}}
\\tablecolumns{}
\\tablehead{
       \colhead{} &
       \colhead{} &
       \multicolumn{2}{$R/R_{25}$} &
       \multicolumn{2}{$\log \left( \Sigma_* {\rm [M_\odot pc^{-2}]}$ ) &
       \multicolumn{2}{$\log \left( \Sigma_{mol} {\rm [M_\odot pc^{-2}]}$ ) \\
       \colhead{} &
       \colhead{\rtwoone} & 
       \colhead{m} &
       \colhead{b} &
       \colhead{m} &
       \colhead{b} &
       \colhead{m} &
       \colhead{b} }
\startdata
'''

footer = '''\endata
\\end{deluxetable}
'''

#----------------------------------------------------------------------
#                              Write Table
#----------------------------------------------------------------------

yvar_list = ['ratio_HCN_CO','ratio_ltir_mean_HCN']
xvar_list = ['r25','mstar','ICO']

for yvar in yvar_list:

    if yvar == 'ratio_HCN_CO':
        yvarstr = 'HCN-to-CO'
    elif yvar == 'ratio_ltir_mean_HCN':
        yvarstr = 'L$_TIR$/HCN'
    else:
        print('bin '+yvar+' not recognized.')
        

        
    # TODO: The below is really fugly. MAKE INTO A FUNCTION. TAKES STACK, XVAR AND YVAR AND GENERATES STRING
    datastr = ''
    data2str = ''
    for xvar in xvar_list:

        # get coefficients for constant R21 assumption
        r21str = 'constant'
        idx = (stack_fits['bin'] == xvar) & (stack_fits['column'] == yvar)        
        tmpstr = "${0:5.2f} \pm {1:5.2f}$ & ${2:5.2f} \pm {3:5.2f}$ & {4} \\\\ \n".format(stack_fits[idx]['slope'][0], stack_fits[idx]['slope_err'][0], stack_fits[idx]['intercept'][0], stack_fits[idx]['intercept_err'][0], stack_fits[idx]['chisq'][0])

        datastr = datastr + tmpstr
        
        
        # get coefficients for spatial R21 assumption
        r21str = 'spatial'
        idx = (stack_fits_spatialR21['bin'] == xvar) & (stack_fits_spatialR21['column'] == yvar)
        tmpstr = "${0:5.2f} \pm {1:5.2f}$ & ${2:5.2f} \pm {3:5.2f}$ & {4} \\\\ ".format(stack_fits_spatialR21[idx]['slope'][0], stack_fits_spatialR21[idx]['slope_err'][0], stack_fits_spatialR21[idx]['intercept'][0], stack_fits_spatialR21[idx]['intercept_err'][0], stack_fits_spatialR21[idx]['chisq'][0])

        data2str = data2str + tmpstr

    print(yvarstr + ' & ' + r21str + ' & ' + datastr)
    print(yvarstr + ' & ' + r21str + ' & ' + data2str)
