import os
import subprocess
import glob
import warnings
from astropy.time import Time
import pkg_resources
import numpy as np
from packaging import version

from . import catalogs
from gbtpipe import ArgusCal
from astropy.coordinates import SkyCoord
import astropy.units as u

# Needed for spatial mask based offs
import astropy.io.fits as fits
import astropy.wcs as wcs

import functools

def reduceAll(release='QA0',
              galaxyList=None,
              update=False,
              overwrite=False,
              outputDir=None,
              getmask=True,
              OffType='linefit',
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
    currentVersion = version.parse(currentVersion)
    catalogs.updateLogs(release=release)
    ObjectCatalog = catalogs.loadCatalog()
    Log = catalogs.parseLog()
    SessionRows = np.where(['Map' in row['Scan Type'] for row in Log])
    Log = Log[SessionRows]

    uniqSrc = list(ObjectCatalog['NAME'])
    if outputDir is None:
        cwd = os.getcwd()
    else:
        cwd = outputDir
    
    if not galaxyList:
        galaxyList = uniqSrc

    for galaxy in uniqSrc:
        if ((galaxy != 'none') and (galaxy in galaxyList)):
            print('Reducing galaxy: '+galaxy+'\n')
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
                        versionstr = pieces[-1]
                        dataVersion = version.parse(versionstr[0:-5])
                        # Cut off .fits
                        if ((dataStart >= row['Start Scan']) and
                            (dataEnd <= row['End Scan'])):
                            if not update:
                                matchdata[idx] = True
                            if update and dataVersion >= currentVersion:
                                matchdata[idx] = True
                    # Assume we need all 16 files to be there for this
                    # to work
                    filesIntact = (matchdata.sum() >= 16)

                # Kill off files that are going to be overwritten
                if overwrite and filesIntact:
                    gottaGo = np.array(ExtantFiles)[matchdata]
                    for thisfile in gottaGo:
                        os.remove(thisfile)

                # Actually do the calibration
                if (overwrite) or (not filesIntact):
                    wrapper(galaxy=galaxy, overwrite = overwrite,
                            release=release, obslog = [row],
                            startdate=row['Date (UT)'],
                            enddate=row['Date (UT)'],
                            project=row['Project'],
                            OffType=OffType,
                            **kwargs)
                
            # if np.any(LogRows):
            #     wrapper(galaxy=galaxy, overwrite = overwrite,
            #             release=release, obslog = Log[LogRows],
            #             startdate=Log[LogRows][0]['Date (UT)'],
            #             enddate=Log[LogRows][-1]['Date (UT)'],
            #             project=Log[LogRows][-1]['Project'],
            #             **kwargs)
            os.chdir(cwd)





def reduceSession(session=1, overwrite=False, release = 'QA0', 
                  project='17B-151',
                  outputDir='/lustre/pipeline/scratch/DEGAS/',
                  OffType='linefit',
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
            except FileNotFoundError:
                os.mkdir(cwd+'/'+galaxy)
                os.chdir(cwd+'/'+galaxy)
            LogRows = Log['Source'] == galaxy
            if np.any(LogRows):
                wrapper(galaxy=galaxy, overwrite = overwrite,
                        release=release, obslog = Log[LogRows],
                        startdate=Log[LogRows][0]['Date (UT)'],
                        enddate=Log[LogRows][-1]['Date (UT)'],
                        project=project,
                        OffType=OffType,
                        **kwargs)
                os.chdir(cwd)

def wrapper(logfile='ObservationLog.csv',galaxy='NGC2903',
            overwrite=False, startdate = '2015-01-1',
            project='17B-151',
            enddate='2020-12-31',release='QA2',obslog=None,
            OffType='linefit',
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
                raise Exception
            

            # TODO: Beam gains?
            if str(observation['Alternate Data Directory']) == '--':
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=startscan,
                           EndScan=endscan,
                           BadScans=observation['Bad Scans'],
                           BadFeeds=observation['Bad Feeds'],
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           Project=project,
                           OffType=OffType,
                           overwrite=overwrite, **kwargs)
            else :
                doPipeline(SessionNumber=observation['Astrid Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           BadScans=observation['Bad Scans'],
                           BadFeeds=observation['Bad Feeds'],
                           RefScans=refscans,
                           Galaxy=observation['Source'],
                           Setup=observation['Setup'],
                           RawDataDir=observation['Special RawDir'],
                           Project=project,
                           OffType=OffType,
                           overwrite=overwrite, **kwargs)

def doPipeline(SessionNumber=7,StartScan = 27, EndScan=44,
               RefScans=None,
               Galaxy='NGC2903', Window='0',
               Project='17B-151',
               Setup='HCN/HCO+',
               BadScans=None,
               BadFeeds=None,
               OptionDict = None,
               RawDataDir = None,
               Gains=None,
               OffType='linefit',
               OutputRoot = None, overwrite=False, 
               MaskName=None,
               CatalogFile=None,
               **kwargs):
    """
    This is the basic DEGAS pipeline which in turn uses the gbt pipeline.
    """
    if RefScans is None:
        RefScans = [25, 45]
    if RawDataDir is None:
        try:
            RawDataDir = os.environ["DEGASDIR"]+'/rawdata/'
        except KeyError:
            RawDataDir = '/lustre/pipeline/scratch/DEGAS/rawdata/'

    if Gains is None:
        Gains = '1,'*16
        Gains = Gains[0:-1]
    SessionDir = 'AGBT' + Project.replace('-','_') + '_' + str(SessionNumber).zfill(2) 
    SessionSubDir = SessionDir + '.raw.vegas/'
    print('Reducing '+SessionDir)
    
    # Set default pipeline options as a dictionary
    if OutputRoot is None:
        try:
            OutputRoot = os.environ["DEGASDIR"]
        except KeyError:
            OutputRoot = '/lustre/pipeline/scratch/DEGAS/'

    # Try to make the output directory
    OutputDirectory = os.path.join(OutputRoot,Galaxy,Setup.replace('/','_'))

    if not os.access(OutputDirectory,os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            print('Made directory {0}'.format(OutputDirectory))
        except OSError:
            try:
                os.mkdir('/'.join((OutputDirectory.split('/'))[0:-1])) # there may be a safer what to do this with os.path.split
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

    ObjectCatalog = catalogs.loadCatalog(CatalogFile=CatalogFile)
    ThisGalaxy = ObjectCatalog[ObjectCatalog['NAME'] == Galaxy]
    galaxy_center = SkyCoord(ThisGalaxy['RA'], ThisGalaxy['DEC'],
                             unit=(u.hour, u.deg), frame='fk5')

    
    if BadScans:
        # issue here is that single value is read as a int, but
        # multiple values as string. 
        if BadScans.isascii():
            # the below may or may not work. We'll find out!!!
            BadScanArray = np.array(BadScans.split(','),dtype=np.int)
        else:
            BadScanArray = np.array([BadScans])

    else:
        BadScanArray = []

    if BadFeeds:
        # issue here is that single value is read as a int, but
        # multiple values as string. 
        if BadFeeds.isascii():
            BadFeedArray = np.array(BadFeeds.split(','),dtype=np.int)
        else:           
            BadFeedArray = np.array([BadFeeds])
    else:
        BadFeedArray = []

    # if MaskName is None:
    #     # AAK: Modified to add 12CO to the mask names since that's
    #     # what my code generates.
    #     MaskName = os.environ["DEGASDIR"] + '/masks/{0}_12CO_mask.fits'.format(Galaxy.upper())

    suffixdict = {'HCN/HCO+':'.hcn_hcop.mask.fits',
                  '13CO/C18O':'.13co_c18o.mask.fits',
                  '12CO': '.12co.mask.fits'}

    if MaskName is None:
        TestMaskName = os.environ["DEGASDIR"] + '/masks/EMISSION_MASKS/' + Galaxy + suffixdict[Setup]
        if os.access(TestMaskName, os.R_OK):
            print('Found Mask:', TestMaskName)
            MaskName = TestMaskName

    if os.access(MaskName, os.R_OK):
        print('Using spatial masking to set OFF positions: '+ MaskName)
        warnings.warn('Using spatial masking to set OFF positions')
        maskhdu = fits.open(MaskName.format(Galaxy.upper()))
        # mask = ~np.any(maskhdu[0].data, axis=0)
        mask = (maskhdu[0].data).astype(np.bool)
        mask = maskhdu[0].data
        w = wcs.WCS(maskhdu[0].header)
        offselect = functools.partial(ArgusCal.SpatialSpectralMask,
                                      mask=mask, wcs=w, floatvalues=True,
                                      offpct=50)

    else:
        warnings.warn('No mask found. Using zone of avoidance masking')
        offselect = functools.partial(ArgusCal.ZoneOfAvoidance,
                                      center=galaxy_center)
                                      
    # offselect = ArgusCal.NoMask
    ArgusCal.calscans(RawDataDir + SessionDir + '/' + SessionSubDir, 
                      start=StartScan,
                      stop=EndScan,
                      refscans=RefScans,
                      badscans=BadScanArray,
                      badfeeds=BadFeedArray,
                      outdir=OutputDirectory,
                      OffSelector=offselect,
                      OffType=OffType,
                      suffix=suffixname,
                      center=galaxy_center,
                      radius=1 * u.arcmin)
    
    # Clean up permissions
    fl = glob.glob(OutputDirectory + '/*fits')
    for thisfile in fl:
        os.chmod(thisfile,0o664)


    
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
                    
