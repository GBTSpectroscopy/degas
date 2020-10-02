# Purpose: Create latex table showing the sample properties for the paper
#
# Date          Programmer              Description of Changes
#----------------------------------------------------------------------
# 3/24/2020     A.A. Kepley             Original Code

import os
from astropy.table import Table
import numpy as np
from importlib import reload
from astropy import units as u
from astropy.coordinates import SkyCoord
import math

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

outfile = os.path.join(os.environ['ANALYSISDIR'],'tables','sample_table.tex')

header = '''\\begin{deluxetable}{rrrrrrrcrrr}
%\\rotate
\\tablewidth{0pt}
\\tabletypesize{\scriptsize}
\\tablecaption{ DEGAS DR1 Sample Properties \label{tab:sample}}
\\tablecolumns{11}
\\tablehead{
        \colhead{} &
        \colhead{R.A. (J2000)} &
        \colhead{Dec. (J2000)} &
        \colhead{Distance} & 
        \colhead{Inclination} &
        \colhead{Position Angle} & 
        \colhead{} &
        \colhead{R$_{25}$ } & 
        \colhead{$\log$ M$_*$} &
        \colhead{$\log$ SFR} \\\\ 
        \colhead{Name} &
        \colhead{h:m:s} &
        \colhead{$\degr$:$\\arcmin$:$\\arcsec$} & 
        \colhead{Mpc} & 
        \colhead{$\degr$} &
        \colhead{$\degr$} &
        \colhead{Ref.} &
        \colhead{$\\arcmin$} &
        \colhead{M$_\odot$} &
        \colhead{M$_\odot$ yr$^{-1}$}
}
\startdata
'''

footer = '''\enddata
\\tablecomments{The values for right ascension, declination, and R$_{25}$ are taken from \citet{Makarov2014HyperLEDA.Distances}. Distances are taken from \citet{Tully2009TheDatabase}  with the exception of NGC6946, whose distance was taken from \citet{Anand2018A6946}. he quantities $\log$~M$_*$ and $\log$~SFR  are taken from \citet{Leroy2019AGalaxies}. The references for the inclination and position angles of the galaxies are indicated in the table as follows: (1) \citet{Leroy2009Heracles:Survey}  (2)\citet{Makarov2014HyperLEDA.Distances}  (3) P. Lang et al., (submitted to ApJ)  (4) \citet{DeBlok2008High-resolutionThings}. }


\\end{deluxetable}
''' 

#----------------------------------------------------------------------
#                              Write table
#----------------------------------------------------------------------

f = open(outfile,'w')

f.write(header)

for galaxy in degas[idx]:
    

    c = SkyCoord(ra=galaxy['RA_DEG']*u.degree,dec=galaxy['DEC_DEG']*u.degree,
                 frame='fk5')

    if galaxy['REF_INCL'] == 'HERACLESVALUE':
        refno = '1'
    elif galaxy['REF_INCL'] == 'LEDA':
        refno = '2'
    elif galaxy['REF_INCL'] == 'LANGMEIDT19':
        refno = '3'
    elif galaxy['REF_INCL'] == 'DEBLOK08':
        refno = '4'
    else:
        print("I don't recognize the literature source for the inclination")

    ## should this code be in degas_base??
    dist_err = (10**galaxy['E_DIST_DEX'] - 1.0) * galaxy['DIST_MPC']

    datastr = "{0} & {1:02d}:{2:02d}:{3:06.3f}  &  {4:= 03d}:{5:02d}:{6:05.2f} &   ${8:4.1f}\pm{9:4.1f}$ & ${11:5.1f}\pm{12:5.1f}$ &   ${14:5.1f}\pm{15:5.1f}$ & {16:s} &   ${17:5.2f}\pm{18:5.2f}$  &  ${20:5.1f}\pm{21:5.1f}$ &   ${20:5.1f}\pm{21:5.1f}$ \\\\ \n".format(galaxy['NAME'], #0
                                                       int(c.ra.hms[0]), # 1
                                                       int(c.ra.hms[1]), # 2
                                                       c.ra.hms[2], #3
                                                       int(c.dec.dms[0]), #4
                                                       abs(int(c.dec.dms[1])), #5
                                                       abs(c.dec.dms[2]), # 6
                                                       galaxy['REF_POS'], #7
                                                       galaxy['DIST_MPC'], #8
                                                       dist_err, ## 9; 
                                                       galaxy['REF_DIST'], #10
                                                       round(galaxy['INCL_DEG'],1), #11
                                                       round(galaxy['E_INCL'],1), #12
                                                       galaxy['REF_INCL'], #13
                                                       round(galaxy['POSANG_DEG'],1), #14
                                                       round(galaxy['E_POSANG'],1), #15
                                                       #galaxy['REF_POSANG'][0].decode('ascii'), #16
                                                       refno, #16
                                                       galaxy['R25_DEG']*60.0, #17 deg to arcmin
                                                       galaxy['E_R25']*60.0, #18; deg to arcmin
                                                       galaxy['REF_R25'], #19
                                                       galaxy['LOGMSTAR'], #20
                                                       galaxy['E_LOGMSTAR'], #21
                                                       galaxy['LOGSFR'],  #22
                                                       galaxy['E_LOGSFR'])  #23  
    
    f.write(datastr)


f.write(footer)

f.close()
