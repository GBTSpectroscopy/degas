import subprocess
import glob
import os
import warnings
from astropy.table import Table, join
import numpy as np
from astropy.utils.data import get_pkg_data_filename

def updateLogs(output='ObservationLog.csv',release=None):
    if 'DR1' in release:
        filename = get_pkg_data_filename('data/ObservationLog_DR1.csv',
                                         package='degas')
        command = 'cp '+filename+' ./ObservationLog.csv'
        return not subprocess.call(command,shell=True)
    if release is None:
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=1wFgxwpDHuYsb7ISWUyBg1xGO6SraRwCIlm1u_fHlrPs&output=csv'"
        # The returns from subprocess are the error codes from the OS
        # If 0 then it worked so we should return True
        return not subprocess.call(command,shell=True)

def loadCatalog(release=None):
    CatalogFile = get_pkg_data_filename('./data/dense_survey.cat',
                                        package='degas')
    Catalog = Table.read(CatalogFile, format='ascii')
    return(Catalog)

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
    idx = ~t['Project'].mask
    t = t[idx]

    # Convert from Google Booleans to Python Booleans
    for drkey in ['QA2','DR1','QA1','QA0']:
        t[drkey] = t[drkey].data.data=='TRUE'

    return(t)
