import os
import subprocess
import glob
import warnings
from astropy.time import Time
import pkg_resources
import numpy as np
from distutils.version import LooseVersion

from . import catalogs
from . import arguspipe

def reduceAll(release='QA0',
              update=False,
              overwrite=False,
              outputDir='/lustre/pipeline/scratch/DEGAS/',
              **kwargs):
    """
    This pulls logs and tries to reduce everything that's not already
    on disk

    Keywords:

    release : str
        Selects which data files to consider for batch reduction.
    
    update : bool 
        Overwrite only if the current pipeline version is larger than
        the version with which the file was generated.

    overwrite: bool
        Overwrite no matter what
        
    """
    currentVersion = pkg_resources.get_distribution("degas").version
    catalogs.updateLogs(release=release)
    ObjectCatalog = catalogs.loadCatalog()
    Log = catalogs.parseLog()
    SessionRows = np.where(['Map' in row['Scan Type'] for row in Log])
    Log = Log[SessionRows]

    uniqSrc = ObjectCatalog['NAME']
    if outputDir is None:
        cwd = os.getcwd()
    else:
        cwd = outputDir
    # TODO: Take out the [::-1]; just there for debug
    for galaxy in uniqSrc[::-1]:
        if galaxy != 'none':
            try:
                os.chdir(cwd+'/'+galaxy)
            except OSError:
                os.mkdir(cwd+'/'+galaxy)
                os.chdir(cwd+'/'+galaxy)
            LogRows = Log['Source'] == galaxy
            # TODO Fill in overwrite functionality and version checking
            for row in Log[LogRows]:
                filesIntact = False # Assume files aren't there
                setup = row['Setup'].replace('/','_')
                querystring = (setup + '/'
                               + '*' + row['Project']
                               + '*' 
                               + 'sess{0}'.format(row['Astrid Session'])
                               + '*fits')
                # Then look for the files
                ExtantFiles = glob.glob(querystring)
                if ExtantFiles:
                    matchdata = np.zeros_like(ExtantFiles, dtype='bool')
                    for idx, thisfile in enumerate(ExtantFiles):
                        pieces = thisfile.split('_')
                        for place, piece in enumerate(pieces):
                            if piece == 'scan':
                                startplace = place
                        dataStart = int(pieces[startplace + 1])
                        dataEnd = int(pieces[startplace + 2])
                        feed = int(pieces[startplace + 4][-1])
                        version = pieces[-1]
                        version = LooseVersion(version[0:-5]) # Cut off .fits
                        if ((dataStart >= row['Start Scan']) and
                            (dataEnd <= row['End Scan'])):
                            if not update:
                                matchdata[idx] = True
                            if update and version >= currentVersion:
                                matchdata[idx] = True
                    # Assume we need all 16 files to be there for this
                    # to work
                    filesIntact = (matchdata.sum() >= 16)
                # Kill off files that are going to be overwritten
                if overwrite and filesIntact:
                    gottaGo = ExtantFiles[matchdata]
                    for thisfile in gottaGo:
                        os.remove(thisfile)

                # Actually do the calibration
                if (overwrite) or (not filesIntact):
                    wrapper(galaxy=galaxy, overwrite = overwrite,
                            release=release, obslog = [row],
                            startdate=row['Date (UT)'],
                            enddate=row['Date (UT)'],
                            project=row['Project'],
                            **kwargs)
                
            # if np.any(LogRows):
            #     wrapper(galaxy=galaxy, overwrite = overwrite,
            #             release=release, obslog = Log[LogRows],
            #             startdate=Log[LogRows][0]['Date (UT)'],
            #             enddate=Log[LogRows][-1]['Date (UT)'],
            #             project=Log[LogRows][-1]['Project'],
            #             **kwargs)
            os.chdir(cwd)


def reduceSession(session=1, overwrite=False, release = 'QA2', 
                  project='17B-151',
                  outputDir='/lustre/pipeline/scratch/DEGAS/',
                  **kwargs):
    """
    Function to reduce single-session data using the GBT-pipeline.
    
    reduceSession(overwrite=False, release='QA2')

    release : string
        Variable that selects which set of data is to be reduced. 
        Default value is 'QA2', while 'DR1' generates the Data Release 1, and 
        hopefully 'DR2' will be available in the near future.
    overwrite : bool
        If True it will overwrite files.
    """

    catalogs.updateLogs(release=release)
    ObjectCatalog = catalogs.loadCatalog()
    Log = catalogs.parseLog()
    SessionRows = ((Log['Astrid Session'] == session) *
                   (Log['Project'] == project))
    Log = Log[SessionRows]
    uniqSrc = ObjectCatalog['NAME']
    if outputDir is None:
        cwd = os.getcwd()
    else:
        cwd = outputDir

    for galaxy in uniqSrc:
        if galaxy != 'none':
            try:
                os.chdir(cwd+'/'+galaxy)
            except OSError:
                os.mkdir(cwd+'/'+galaxy)
                os.chdir(cwd+'/'+galaxy)
            LogRows = Log['Source'] == galaxy
            if np.any(LogRows):
                wrapper(galaxy=galaxy, overwrite = overwrite,
                        release=release, obslog = Log[LogRows],
                        startdate=Log[LogRows][0]['Date (UT)'],
                        enddate=Log[LogRows][-1]['Date (UT)'],
                        project=project,
                        **kwargs)
                os.chdir(cwd)

def wrapper(logfile='ObservationLog.csv',galaxy='NGC2903',
            overwrite=False, startdate = '2015-01-1',
            project='17B-151',
            enddate='2020-12-31',release='QA2',obslog=None,
            **kwargs):
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
                (observation[release]) & \
                ('Map' in observation['Scan Type']):

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
                warnings.warn("Number of scans = Number of Mapped Rows: no VaneCal")
                raise
            

            # TODO: Beam gains?
            if str(observation['Alternate Data Directory']) == '--':
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=startscan,
                           EndScan=endscan,
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           Project=project,
                           overwrite=overwrite, **kwargs)
            else :
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           RawDataDir=observation['Special RawDir'],
                           Project=project,
                           overwrite=overwrite, **kwargs)

def doPipeline(SessionNumber=7,StartScan = 27, EndScan=44,
               RefScans=[25, 45],
               Galaxy='NGC2903', Window='0',
               Project='17B-151',
               Setup='HCN/HCO+',
               OptionDict = None,
               RawDataDir = None,
               Gains=None,
               OffType='median',
               OutputRoot = None, overwrite=False, **kwargs):
    """
    This is the basic DEGAS pipeline which in turn uses the gbt pipeline.
    """
    if RawDataDir is None:
        # RawDataDir = '/home/sdfits/'
        RawDataDir = '/lustre/pipeline/scratch/DEGAS/rawdata/'
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
    suffixname='_{2}_sess{0}_v{1}'.format(SessionNumber,
                                          pkg_resources.get_distribution("degas").version,
                                          Project)

    arguspipe.calscans(RawDataDir + SessionDir + '/' + SessionSubDir, 
                       start=StartScan,
                       stop=EndScan,
                       refscans=RefScans,
                       outdir=OutputDirectory,
                       OffType=OffType,
                       suffix=suffixname)
                                       
    # Clean up permissions
    fl = glob.glob(OutputDirectory + '/*fits')
    for thisfile in fl:
        os.chmod(thisfile,0664)


    
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
                    
