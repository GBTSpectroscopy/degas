# Purpose: create latex table showing the source of the IR data for
# each galaxy in DR1.

import os
from astropy.table import Table
import numpy as np
from importlib import reload

#----------------------------------------------------------------------
#                           read in data base
#----------------------------------------------------------------------

scriptDir = os.environ['SCRIPTDIR']
degas = Table.read(os.path.join(scriptDir,'degas_base.fits'),format='fits')

degas.sort(['RA_DEG','DEC_DEG'])

idx = degas['DR1'] == 1

#----------------------------------------------------------------------
#                               prep table
#----------------------------------------------------------------------

outfile = os.path.join(os.environ['ANALYSISDIR'],'tables','IR_data_table.tex')

header = '''\\begin{deluxetable}{rrrrrrrr}
\\tablewidth{0pt}
\\tabletypesize{\scriptsize}
\\tablecaption{Source of Ancillary Data \label{tab:ir_data}}
\\tablecolumns{8}
\\tablehead{
        \colhead{} & 
        \colhead{} &
        \colhead{} &
        \colhead{} & 
        \multicolumn{4}{c}{\LTIR} \\\\
        \colhead{Name} & 
        \colhead{\\twelveCO} &
        \colhead{Transition} &
        \colhead{\LTIR($24\micron$)} & 
        \colhead{$24\micron$} &
        \colhead{$70\micron$} & 
        \colhead{$100\micron$} &
        \colhead{$160\micron$} }

\startdata
'''

footer = '''\enddata
\\tablenotetext{a}{Data at lower resolution (30\\arcsec) used to color correct 24\micron\ data.}
\\tablecomments{Heracles: \citet{Leroy2009Heracles:Survey}, GBT: Li et al. in prep, PHANGS: \citet{Leroy2021}, everyHERACLES: Schruba et al. in prep, Wilson: Wilson et al. in prep, BENDO: \citet{Bendo2012}, LVL: \citet{Dale2009}, SINGS: \citet{Kennicutt2003}, KINGFISH: \citet{Kennicutt2011}, HRS: \citet{Boselli2010}, VNGS: \citet{Bendo2012a}, z0mgs: \citet{Leroy2019}}
\\end{deluxetable}
'''
#----------------------------------------------------------------------
#                              Write table
#----------------------------------------------------------------------

f = open(outfile,'w')

f.write(header)

irbands = ['IR_24micron_only','IR_24micron','IR_70micron','IR_100micron','IR_160micron']

for galaxy in degas[idx]:

    # indicate which galaxies have no data for a particular band
    for band in irbands:
        if not galaxy[band]:
            galaxy[band] = '\\nodata'

    if galaxy['MASK'] == 'everyHERACLES_Andreas':
        codata = 'everyHERACLES'
    else:
        codata = galaxy['MASK']
    
    if ((codata == 'everyHERACLES') or (codata == 'PHANGS') or (codata == 'HERACLES')):
        coline = '\co{12}{2}{1}'
    else:
        coline = '\co{12}{1}{0}'
  
    if galaxy['IR_70micron'] == 'LVL_gauss30':
        IR_70micron = 'LVL \\tablenotemark{a}'
    else:
        IR_70micron = galaxy['IR_70micron']

    line = "{0} & {1} & {2} & {3}& {4} & {5}& {6} & {7}\\\\ \n".format(galaxy['NAME'],
                                                       codata,
                                                       coline,
                                                       galaxy['IR_24micron_only'],
                                                       galaxy['IR_24micron'],
                                                       IR_70micron,
                                                       galaxy['IR_100micron'],
                                                       galaxy['IR_160micron'])

    f.write(line)

f.write(footer)
f.close()
