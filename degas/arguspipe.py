import gbtpipe
import numpy as np
import glob
import os
import fitsio
import copy
import warnings

# CRUFT
#/users/rmaddale/bin/getForecastValues -freqList 89 -typeList Opacity -elev 90 -timeList 53441.4
# tau = w.retrieve_zenith_opacity(53441.4,88e9,log=log)
# print tau
# tau = w.retrieve_zenith_opacity(53441.4,89e9,log=log,forcecalc=True)
# print tau


def makelogdir():
    if not os.path.exists('log'):
        os.mkdir('log')


def gettsys(cl_params, row_list, thisfeed, thispol, thiswin, pipe,
            weather=None, log=None):
    if not weather:
        weather = gbtpipe.Weather()
    if not log:
        log = gbtpipe.Logging(cl_params, 'gbtpipeline')
    cal = gbtpipe.Calibration()

    #Assume refscan is this scan and the next scan
    thisscan = cl_params.refscans[0]
    row1 = row_list.get(thisscan, thisfeed, thispol, thiswin)
    row2 = row_list.get(thisscan+1, thisfeed, thispol, thiswin)

    ext = row1['EXTENSION']
    rows = row1['ROW']
    columns = tuple(pipe.infile[ext].get_colnames())
    integ1 = gbtpipe.ConvenientIntegration(pipe.infile[ext][columns][rows], log=log)

    ext = row2['EXTENSION']
    rows = row2['ROW']
    columns = tuple(pipe.infile[ext].get_colnames())
    integ2 = gbtpipe.ConvenientIntegration(pipe.infile[ext][columns][rows], log=log)

    vec1 = integ1.data['DATA']
    vec2 = integ2.data['DATA']
    if 'Vane' in integ1.data['CALPOSITION'][0] and\
            'Observing' in integ2.data['CALPOSITION'][0]:
        onoff = np.median((vec1-vec2)/vec2)
        # Mean over time to get a vector of vaneCounts at each freq.
        vaneCounts = np.mean(vec1, axis=0)
    elif 'Vane' in integ2.data['CALPOSITION'][0] and\
            'Observing' in integ1.data['CALPOSITION'][0]:
        onoff = np.median((vec2-vec1)/vec1)
        vaneCounts = np.mean(vec2, axis=0)
    else:
        import pdb; pdb.set_trace()

    timestamps = integ1.data['DATE-OBS']
    mjd = np.mean(np.array([pipe.pu.dateToMjd(stamp)
                            for stamp in timestamps]))

    elevation = np.mean(integ1.data['ELEVATIO'])

    twarm = np.mean(integ1.data['TWARM']+273.15)
    tbg = 2.725
    tambient = np.mean(integ1.data['TAMBIENT'])
    avgfreq = np.mean(integ1.data['OBSFREQ'])
    tatm = weather.retrieve_Tatm(mjd, avgfreq, log=log,
                                 forcecalc=True)
    zenithtau = weather.retrieve_zenith_opacity(mjd, avgfreq, log=log,
                                                forcecalc=True,
                                                request='Opacity')

    tau = cal.elevation_adjusted_opacity(zenithtau, elevation)
    tcal = (tatm - tbg) + (twarm - tatm) * np.exp(tau)
    tsys = tcal / onoff
    return tcal, vaneCounts, tsys


#cl_params.mapscans = list(np.linspace(113,136,136-113+1).astype('int'))
#cl_params.refscans = [111,]


def calscans(inputdir, start=82, stop=105, refscans = [80],
             outdir=None, log=None, loglevel='warning',
             OffFrac=0.25, OffType='linefit'):

    if os.path.isdir(inputdir):
        fitsfiles = glob.glob(inputdir + '/*fits')
        if len(fitsfiles) == 0:
            warnings.warn("No FITS files found in input directory")
            return False
    elif os.path.isfile(inputdir):
        warnings.warn("Input name is a file and not a directory.")
        warnings.warn("Blindly reducing everything in this directory.")
        infilename = inputdir
        inputdir = os.getcwd()
    else:
        warnings.warn("No file or directory found for inputs")
        return False

    if not outdir:
        outdir = os.getcwd()

    makelogdir()
    if not log:
        log = gbtpipe.Logging('gbtpipeline')
        if loglevel=='warning':
            log.logger.setLevel(30)

    w = gbtpipe.Weather()
    cl_params = gbtpipe.initParameters(inputdir)
    cl_params.mapscans = list(np.linspace(start,stop,
                                          stop-start+1).astype('int'))
    cl_params.refscans = refscans
    if os.path.isdir(cl_params.infilename):
        log.doMessage('INFO', 'Infile name is a directory')
        input_directory = cl_params.infilename

        # Instantiate a SdFits object for I/O and interpreting the
        # contents of the index file
        sdf = gbtpipe.SdFits()

        # generate a name for the index file based on the name of the
        # raw SDFITS file.  The index file simply has a different extension
        directory_name = os.path.basename(cl_params.infilename.rstrip('/'))
        indexfile = cl_params.infilename + '/' + directory_name + '.index'
        try:
            # create a structure that lists the raw SDFITS rows for
            #  each scan/window/feed/polarization
            row_list, summary = sdf.parseSdfitsIndex(indexfile,
                                                     cl_params.mapscans)
        except IOError:
            log.doMessage('ERR', 'Could not open index file', indexfile)
            log.close()
            return False
            # sys.exit()

    for infilename in glob.glob(input_directory + '/' +
                                os.path.basename(input_directory) +
                                '*.fits'):
            log.doMessage('DBG', 'Attempting to calibrate',
                          os.path.basename(infilename).rstrip('.fits'))

            # change the infilename in the params structure to the
            # current infile in the directory for each iteration
            cl_params.infilename = infilename

            # copy the cl_params structure so we can modify it during
            # calibration for each seperate file.
            command_options = copy.deepcopy(cl_params)
            sdf = gbtpipe.SdFits()
            cal = gbtpipe.Calibration()
            indexfile = sdf.nameIndexFile(command_options.infilename)
            row_list, summary = sdf.parseSdfitsIndex(indexfile,
                                                     mapscans=command_options.mapscans)
            feedlist = (row_list.feeds())
            for thisfeed in feedlist:
                thispol = 0 # TODO: generalize to different POL/WIN
                thiswin = 0
                command_options = copy.deepcopy(cl_params)
                try:
                    pipe = gbtpipe.MappingPipeline(command_options,
                                                   row_list,
                                                   thisfeed,
                                                   thispol,
                                                   thiswin,
                                                   None, outdir=outdir)
                except KeyError:
                    pipe = None
                    pass
                if pipe:
                    tcal, vaneCounts, tsys = gettsys(cl_params, row_list,
                                                     thisfeed, thispol,
                                                     thiswin, pipe,
                                                     weather=w, log=log)
                    log.doMessage('INFO', 'Feed: {0}, Tsys (K): {1}'.format(
                            thisfeed,
                            tsys))

                    for thisscan in cl_params.mapscans:
                        print("Now Processing Scan {0} for Feed {1}".format(\
                                thisscan, thisfeed))
                        rows = row_list.get(thisscan, thisfeed,
                                            thispol, thiswin)
                        ext = rows['EXTENSION']
                        rows = rows['ROW']
                        columns = tuple(pipe.infile[ext].get_colnames())
                        integs = gbtpipe.ConvenientIntegration(\
                            pipe.infile[ext][columns][rows], log=log)
                        timestamps = integs.data['DATE-OBS']
                        elevation = np.median(integs.data['ELEVATIO'])
                        mjds = np.array([pipe.pu.dateToMjd(stamp)
                                         for stamp in timestamps])
                        avgfreq = np.median(integs.data['OBSFREQ'])
                        zenithtau = w.retrieve_zenith_opacity(np.median(mjds),
                                                              avgfreq,
                                                              log=log,
                                                              forcecalc=True)
                        tau = cal.elevation_adjusted_opacity(zenithtau,
                                                             elevation)

                        # Do the calibration here
                        ON = integs.data['DATA']

                        if OffType == 'median':
                            OFF = np.median(ON, axis=0) # Empirical bandpass
                            OFF.shape += (1,)
                            OFF = OFF * np.ones((1, ON.shape[0]))
                            OFF = OFF.T
                        if OffType == 'linefit':
                            xaxis = np.linspace(-0.5,0.5,ON.shape[0])
                            xaxis.shape += (1,)
                            xaxis = xaxis * np.ones((1,ON.shape[1]))
                            cutidx = np.int(OffFrac * ON.shape[0])
                            xsub = np.r_[xaxis[0:cutidx,:], xaxis[-cutidx:,:]]
                            ONsub = np.r_[ON[0:cutidx,:], ON[-cutidx:,:]]

                            MeanON = np.nanmean(ONsub, axis=0)
                            MeanON.shape += (1,)
                            MeanON = MeanON * np.ones((1, ONsub.shape[0]))
                            MeanON = MeanON.T

                            MeanX = np.nanmean(xsub, axis=0)
                            MeanX.shape += (1,)
                            MeanX = MeanX * np.ones((1, ONsub.shape[0]))
                            MeanX = MeanX.T

                            slope = (np.nansum(((xsub - MeanX)
                                                * (ONsub - MeanON)),
                                               axis=0)
                                     / np.nansum((xsub - MeanX)**2,
                                                 axis=0))
                            slope.shape += (1,)
                            slope = slope * np.ones((1, ON.shape[0]))
                            MeanON = MeanON[0,:]
                            MeanON.shape += (1,)
                            MeanON = MeanON * np.ones((1,ON.shape[0]))
                            MeanON = MeanON.T
                            OFF = slope.T * xaxis + MeanON
                        medianOFF = np.nanmedian(OFF, axis=0)

                        # Now construct a scalar factor by taking
                        # median OFF power (over time) and compare to
                        # the mean vaneCounts over time.  Then take
                        # the median of this ratio and apply.

                        scalarOFFfactor = np.median(medianOFF /
                                                    (vaneCounts - medianOFF))
                        TA = (tcal * scalarOFFfactor
                              * (ON - OFF) / (OFF))

                        medianTA = np.median(TA, axis=1)
                        medianTA.shape = (1,) + medianTA.shape
                        medianTA = np.ones((ON.shape[1], 1)) * medianTA
                        TAstar = TA - medianTA.T
                        for ctr, row in enumerate(rows):
                            row = gbtpipe.Integration(
                                pipe.infile[ext][columns][row])
                            row.data['DATA'] = TAstar[ctr,:]
                            row.data['TSYS'] = tsys
                            row.data['TUNIT7'] = 'Ta*'
                            pipe.outfile[-1].append(row.data)
                    pipe.outfile.close()
    return True
