# Purpose: create a copy of the degas data base and insert the EMPIRE
# parameters for the overlap galaxies.

import os
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

scriptDir = os.environ['SCRIPTDIR']

inFile = os.path.join(scriptDir, 'degas_base.fits')
galpropFile = os.path.join(os.environ['ANALYSISDIR'],'ancillary_data','empire','galaxy_prop.txt')


outFile = os.path.join(scriptDir, 'degas_base_empireparams.fits')


degas = Table.read(inFile)
galprop = Table.read(galpropFile,format='ascii.commented_header')

degas_new = degas

for gal in galprop:
    if gal['Galaxy'] in degas_new['NAME']:
        degas_new['DIST_MPC'][degas_new['NAME'] == gal['Galaxy']] = gal['D']
        degas_new['INCL_DEG'][degas_new['NAME'] == gal['Galaxy']] = gal['i']
        degas_new['POSANG_DEG'][degas_new['NAME'] == gal['Galaxy']] = gal['P.A.']
        degas_new['R25_DEG'][degas_new['NAME'] == gal['Galaxy']] = gal['r_25'] / 60.0 # convert from arcmin to deg

        c = SkyCoord(gal['R.A.'],gal['Decl.'],frame='fk5',unit=(u.hourangle,u.deg))
        degas_new['RA_DEG'][degas_new['NAME'] == gal['Galaxy']] = c.ra.value
        degas_new['DEC_DEG'][degas_new['NAME'] == gal['Galaxy']] = c.dec.value


degas_new.write(outFile)
