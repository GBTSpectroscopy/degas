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

def format_data_str(stack_fits, xvar, yvar):
    '''
    Purpose: output data string
    
    Input:
        stack: stack
        xvar: x variable
        yvar: y variable.

    Output:
        string formatted for appropriate output

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    9/29/2022   A.A. Kepley     Original Code
    '''

    idx = (stack_fits['bin'] == xvar) & (stack_fits['column'] == yvar)
    tmpstr = "${0:5.2f} \pm {1:5.2f}$ & ${2:5.2f} \pm {3:5.2f}$ & {4:5.2f}".format(stack_fits[idx]['slope'][0], stack_fits[idx]['slope_err'][0], stack_fits[idx]['intercept'][0], stack_fits[idx]['intercept_err'][0], stack_fits[idx]['chisq'][0])

    return tmpstr



#----------------------------------------------------------------------
#                       read in fit results
#----------------------------------------------------------------------

analysis_dir = os.environ['ANALYSISDIR']
stack_fits= Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1','stack_fits.fits'))
stack_fits_noNGC4414 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1','stack_fits_noNGC4414.fits'))
stack_fits_spatialR21 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21','stack_fits.fits'))
stack_fits_spatialR21_noNGC4414 = Table.read(os.path.join('/lustre/cv/users/akepley/degas/stack_IR6p1_spatialR21','stack_fits_noNGC4414.fits'))

r21_list = ['constant','spatial']

#----------------------------------------------------------------------
#                       Put together table header
#----------------------------------------------------------------------

outfile = os.path.join(os.environ['ANALYSISDIR'],'tables','degas_dr1_fits.tex')

header = '''\\begin{deluxetable}{llrrrrrrrrr}
\\tablewidth{0pt}
\\tabletypesize{\scriptsize}
\\tablecaption{Fit Coefficients \label{tab:fit_coeff}}
\\tablecolumns{11}
\\tablehead{
       \colhead{} &
       \colhead{} &
       \multicolumn{3}{c}{$R/R_{25}$\\tablenotemark{a}} &
       \multicolumn{3}{c}{$\log \left( \Sigma_* {\\rm [M_\odot \, pc^{-2}]} \\right)$ } &
       \multicolumn{3}{c}{$\log \left( I_{CO} {\\rm [K \, km \, s^{-1} pc^{2}]} \\right)$ } \\\\
       \colhead{} &
       \colhead{\\rtwoone} & 
       \colhead{m} &
       \colhead{b} &
       \colhead{$\chi^2$} &
       \colhead{m} &
       \colhead{b} &
       \colhead{$\chi^2$} &
       \colhead{m} &
       \colhead{b} &
       \colhead{$\chi^2$} }
\startdata
'''

footer = '''\enddata
\\tablenotetext{a}{These fits exclude NGC4414 since it it an extreme outlier on these plots.} 
\\end{deluxetable}
'''

#       \multicolumn{3}{c}{$\log \left( \Sigma_{mol} {\\rm [M_\odot pc^{-2}]} \\right)$ } \\\\

#----------------------------------------------------------------------
#                              Write Table
#----------------------------------------------------------------------

fout = open(outfile,'w')
fout.write(header)


yvar_list = ['ratio_HCN_CO','ratio_ltir_mean_HCN']
xvar_list = ['r25','mstar','ICO']

for yvar in yvar_list:

    if yvar == 'ratio_HCN_CO':
        yvarstr = 'HCN-to-CO'
    elif yvar == 'ratio_ltir_mean_HCN':
        yvarstr = 'L$_{TIR}$/HCN'
    else:
        print('bin '+yvar+' not recognized.')

    for r21 in r21_list:
        datastr = ''
        for xvar in xvar_list:
            if datastr != '':
                datastr = datastr + '  &  '

            if (r21 == 'constant') & (xvar == 'r25'):
                mystack = stack_fits_noNGC4414
            elif r21 == 'constant':
                mystack = stack_fits
            elif (r21 == 'spatial')& (xvar == 'r25'):
                mystack = stack_fits_spatialR21_noNGC4414
            elif (r21 == 'spatial' ):
                mystack = stack_fits_spatialR21
            else:
                print("r21 value not recognized: " +r21)

            tmpstr = format_data_str(mystack, xvar, yvar)
            datastr = datastr + tmpstr       
                
        fout.write(yvarstr + ' & ' + r21 + ' & ' + datastr + '  \\\\ \n')


fout.write(footer)
fout.close()
