# Purpose: Generate figure set of stacks for paper
#
# Date          Programmer              Original Code
#----------------------------------------------------------------------
# 4/29/2022     A.A. Kepley             Original Code

import os
import glob
from astropy.table import Table

fignum = '2'

plotDir = os.path.join(os.environ['ANALYSISDIR'],'plots')

degas_db = Table.read(os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits'))

idx = degas_db['DR1'] == 1


## do I want to use the degas data base here? Or go off copied plots. probably degas data base
# stacklist = glob.glob(os.path.join(plotDir,'*_

f = open(os.path.join(plotDir,'stack_fig_set.tex'),'w')

f.write("\\figsetstart \n")
f.write("\\figsetnum{"+fignum+"}\n")
f.write("DEGAS DR1 Stacks\n")

i = 1
for galaxy in degas_db[idx]['NAME']:
    
    if galaxy == 'IC0342':
        continue

    filename = galaxy+'_radius_stacks.pdf'
    filecaption = 'Stacks in 15\arcsec\ wide radial bins for ' + galaxy +'.'

    f.write("\\figsetgrpstart\n")
    f.write("\\figsetgrpnum{"+fignum+"."+str(i)+"}\n")
    f.write("\\figsetgrptitle{"+galaxy+"}\n")
    f.write("\\figsetplot{"+filename+"}\n")
    f.write("\\figsetgrpnote{"+filecaption+"}\n")
    f.write("\\fitsetgrpend\n")

f.write("\\fitsetend\n")

f.close()

