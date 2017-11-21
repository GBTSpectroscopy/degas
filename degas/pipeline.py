import os
import subprocess
import glob
import warnings
from astropy.time import Time
from . import catalogs
from . import arguspipe
import numpy as np

def reduceSession(session=1, overwrite=False, release = 'all', 
                  project='17B-151',
                  outputDir='/lustre/pipeline/scratch/DEGAS/'):
    """
    Function to reduce all data using the GBT-pipeline.
    
    reduceAll(overwrite=False, release='all')

    release : string
        Variable that selects which set of data is to be reduced. 
        Default value is 'all', while 'DR1' generates the Data Release 1, and 
        hopefully 'DR2' will be available in the near future.
    overwrite : bool
        If True it will overwrite files.
    """

    catalogs.updateLogs(release=release)
    RegionCatalog = catalogs.loadCatalog()
    Log = catalogs.parseLog()
    SessionRows = (Log['Session'] == session) * (Log['Project'] == project)
    Log = Log[SessionRows]
    uniqSrc = RegionCatalog['NAME']
    if outputDir is None:
        cwd = os.getcwd()
    else:
        cwd = outputDir

    for region in uniqSrc:
        if region != 'none':
            try:
                os.chdir(cwd+'/'+region)
            except OSError:
                os.mkdir(cwd+'/'+region)
                os.chdir(cwd+'/'+region)
            LogRows = Log['Source'] == region
            if np.any(LogRows):
                wrapper(region=region, overwrite = overwrite,
                        release=release, obslog = LogRows,
                        startdate=Log[LogRows][0]['Date (UT)'],
                        enddate=Log[LogRows][-1]['Date (UT)'])
                os.chdir(cwd)

def wrapper(logfile='ObservationLog.csv',galaxy='NGC2903',
            overwrite=False, startdate = '2015-01-1',
            enddate='2020-12-31',release='all',obslog = None):
    """
    This is the DEGAS pipeline which chomps the observation logs and
    then batch calibrates the data.  It requires AstroPy because
    their tables are pretty nifty.

    wrapper(logfile='../ObservationLog.csv', galaxy='NGC2903')

    galaxy : string
        Galaxy name as given in logs
    logfile : string
        Full path to CSV version of the logfile (optional)
    obslog : astropy.Table
        Table representing an already parsed observation log
    overwrite : bool
        If True, carries out calibration for files already present on disk.
    startdate : string
        representation of date in format YYYY-MM-DD for beginning calibration
    enddate : string
        date in format YYYY-MM-DD for ending calibration
    release : string
        name of column in the log file that is filled with boolean
        values indicating whether a given set of scans belongs to the data
        release.
    If a logfile or obslog isn't specified, logs will be retrieved from Google.
    """

    StartDate = Time(startdate)
    EndDate = Time(enddate)
    if not os.access(logfile,os.R_OK):
        catalogs.updateLogs(release=release)

    if obslog is None:
        t = catalogs.parseLog(logfile=logfile)
    else:
        t = obslog

    for observation in t:
        ObsDate = Time(observation['Date (UT)'])
        if (galaxy in observation['Source']) & \
                (ObsDate >= StartDate) & (ObsDate <= EndDate) & \
                (observation[release]):
            print(observation['Date (UT)'])
            Nscans = observation['End Scan'] - observation['Start Scan'] + 1
            if Nscans - observation['Nrows'] == 4:
                warnings.warn("Assuming Vane Cal at start and end")
                refscans = [observation['Start Scan'],
                            observation['End Scan'] -1]
                startscan = observation['Start Scan'] + 2
                endscan = observation['End Scan'] - 2
            if Nscans - observation['Nrows'] == 2:
                warnings.warn("Assuming Vane Cal at start only")
                refscans = [observation['Start Scan']]
                startscan = observation['Start Scan'] + 2
                endscan = observation['End Scan']
            if Nscans - observation['Nrows'] < 2:
                warnings.warn("Number of scans = Number of Mapped Rows: no Vane Cal")
                raise

            # TODO: Beam gains?
            if str(observation['Alternate Data Directory']) == '--':
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=startscan,
                           EndScan=endscan,
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           overwrite=overwrite)
            else :
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           RawDataDir=observation['Special RawDir'],
                           overwrite=overwrite)

def doPipeline(SessionNumber=7,StartScan = 27, EndScan=44,
               RefScans=[25, 45],
               Galaxy='NGC2903', Window='0',
               Project='17B-151',
               Setup='HCN/HCO+',
               OptionDict = None,
               RawDataDir = None,
               Gains=None,
               OutputRoot = None, overwrite=False):
    """
    This is the basic DEGAS pipeline which in turn uses the gbt pipeline.
    """
    if RawDataDir is None:
        RawDataDir = '/home/sdfits/'
        # RawDataDir = '/lustre/pipeline/scratch/DEGAS/rawdata/'
    if Gains is None:
        Gains = '1,'*16
        Gains = Gains[0:-1]
    SessionDir = 'AGBT' + Project.replace('-','_') + '_' + str(SessionNumber).zfill(2) 
    SessionSubDir = SessionDir + '.raw.vegas/'
    print('Reducing '+SessionDir)
                  
    # Set default pipeline options as a dictionary
    if OutputRoot is None:
        OutputRoot = '/lustre/pipeline/scratch/DEGAS/'
    # Try to make the output directory
    OutputDirectory = OutputRoot + Galaxy + '/' + Setup.replace('/','_')

    if not os.access(OutputDirectory,os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            print('Made directory {0}'.format(OutputDirectory))
        except OSError:
            try:
                os.mkdir('/'.join((OutputDirectory.split('/'))[0:-1]))
                os.mkdir(OutputDirectory)
                print('Made directory {0}'.format(OutputDirectory))
            except:
                warnings.warn('Unable to make output directory '+OutputDirectory)
                raise
        except:
            warnings.warn('Unable to make output directory '+OutputDirectory)
            raise

    arguspipe.calscans(RawDataDir + SessionDir + '/' + SessionSubDir, 
                       start=StartScan,
                       stop=EndScan,
                       refscans=RefScans,
                       outdir=OutputDirectory)

    # Clean up permissions
    fl = glob.glob(OutputDirectory + '/*fits')
    for thisfile in fl:
        os.chown(fl,0774)


    
    # for bank in BankNames:
    #     # Loop over each feed and polarization
    #     # we check if a pipeline call is necessary.
    #     for feed in ['0','1','2','3','4','5','6']:
    #         for pol in ['0','1']:
    #             FilesIntact = True
    #             if not overwrite:
    #                 outputfile = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.fits'.\
    #                     format(StartScan,EndScan,Window,feed,pol,SessionNumber)
    #                 FilesIntact = FilesIntact and os.path.exists(OutputDirectory+'/'+outputfile)
    #                 if FilesIntact:
    #                     print('Data for Polarization {0} of Feed {1} appear on disk... skipping'.format(pol,feed))
    #             #
    #             if (not FilesIntact) or (overwrite):
    #                 InputFile = RawDataDir+SessionDir+'AGBT16B_278_'+\
    #                     str(SessionNumber).zfill(2)+\
    #                     '.raw.vegas.{0}.fits'.format(bank)
    #                 command = 'gbtpipeline -i '+InputFile
    #                 for key in OptionDict:
    #                     command = command+' '+key+' '+OptionDict[key]
    #                 command = command+' --feed '+feed+' --pol '+pol
    #                 print(command)
    #                 subprocess.call(command,shell=True)

    #                 indexname    = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.index'.\
    #                     format(StartScan,EndScan,Window,feed,pol)
    #                 outindexname = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.index'.\
    #                     format(StartScan,EndScan,Window,feed,pol,SessionNumber)
    #                 try:
    #                     os.rename(indexname,OutputDirectory+'/'+outindexname)
    #                 except:
    #                     pass

    #                 filename   = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.fits'.\
    #                     format(StartScan,EndScan,Window,feed,pol)
    #                 outputfile = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.fits'.\
    #                     format(StartScan,EndScan,Window,feed,pol,SessionNumber)
    #                 try:
    #                     os.rename(filename,OutputDirectory+'/'+outputfile)
    #                     os.chown(OutputDirectory+'/'+outputfile,0774)
    #                 except:
    #                     pass
                    
