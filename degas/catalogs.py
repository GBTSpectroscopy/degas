import subprocess
import glob
import os
import warnings
from astropy.table import Table, join, Column
import numpy as np
from astropy.utils.data import get_pkg_data_filename
import sys

def updateLogs(output='ObservationLog.csv',release=None):
    if release is None:
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=1wFgxwpDHuYsb7ISWUyBg1xGO6SraRwCIlm1u_fHlrPs&output=csv'"
        # The returns from subprocess are the error codes from the OS
        # If 0 then it worked so we should return True
        return not subprocess.call(command,shell=True)

    if 'DR1' in release:
        filename = get_pkg_data_filename('data/ObservationLog_DR1.csv',
                                         package='degas')
        command = 'cp '+filename+' ./ObservationLog.csv'
        return not subprocess.call(command,shell=True)

def loadCatalog(release=None, CatalogFile=None):
    if CatalogFile is None:
        CatalogFile = get_pkg_data_filename('./data/dense_survey.cat',
                                            package='degas')
    Catalog = Table.read(CatalogFile, format='ascii')
    return(Catalog)

def validateScanNumbers(catalog, release='QA0'):

    for observation in catalog:
        if ('Map' in observation['Scan Type']) and (observation[release]):
            Nscans = observation['End Scan'] - observation['Start Scan'] + 1
            if Nscans - observation['Nrows'] == 4:
                warnings.warn("Assuming Vane Cal at start and end")
                refscans = [observation['Start Scan'],
                            observation['End Scan'] -1]
                startscan = observation['Start Scan'] + 2
                endscan = observation['End Scan'] - 2
            elif Nscans - observation['Nrows'] == 2:
                warnings.warn("Assuming Vane Cal at start only")
                refscans = [observation['Start Scan']]
                startscan = observation['Start Scan'] + 2
                endscan = observation['End Scan']
            elif Nscans - observation['Nrows'] == 0:
                warnings.warn("Number of scans = Number of Mapped Rows: no VaneCal")
            else:
                warnings.warn("Inconsistent number of scan rows")
                print(observation)
                raise
            
def parseLog(logfile='ObservationLog.csv'):
    """
    Ingests a CSV log file into an astropy table
    """
    try:
        from astropy.table import Table
    except:
        warnings.warn('DEGAS Pipeline requires astropy.  Try Anaconda!')
        return
    t = Table.read(logfile)

    # Cull to full rows
    # idx = ~t['Project'].mask
    # t = t[idx]

    # Convert from Google Booleans to Python Booleans

    qa0 = np.zeros(len(t), dtype=bool)
    dr1 = np.zeros(len(t), dtype=bool)
    for idx, row in enumerate(t):
        qa0[idx] = ('TRUE' in row['QA0'])
        dr1[idx] = ('TRUE' in row['DR1'])
    qa0col = Column(qa0, dtype=bool, name='QA0')
    dr1col = Column(dr1, dtype=bool, name='DR1')

    t.replace_column('QA0', qa0col)
    t.replace_column('DR1', dr1col)

    for row in t:
        if (('Map' in row['Scan Type'])
            and ('science' in row['Source Type'])
            and (row['QA0'])):
            try:
                first = int(row['First Row Mapped'])
                last = int(row['Last Row Mapped'])
            except:
                print(row)
    validateScanNumbers(t)
    return(t)
