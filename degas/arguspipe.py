import gbtpipe
import numpy as np
import glob
import os
import fitsio
import copy
from gbtpipe import SdFits

# CRUFT
#/users/rmaddale/bin/getForecastValues -freqList 89 -typeList Opacity -elev 90 -timeList 53441.4
# tau = w.retrieve_zenith_opacity(53441.4,88e9,log=log)
# print tau
# tau = w.retrieve_zenith_opacity(53441.4,89e9,log=log,forcecalc=True)
# print tau


def gettsys(cl_params, row_list, thisfeed, thispol, thiswin, pipe,
            weather=None, log=None):
    if not weather:
        weather = gbtpipe.Weather()
    if not log:
        log = gbtpipe.Logging()
    cal = gbtpipe.Calibration()

    #Assume refscan is this scan and the next scan
    thisscan = cl_params.refscans[0]
    row1 = row_list.get(thisscan, thisfeed, thispol, thiswin)
    row2 = row_list.get(thisscan+1, thisfeed, thispol, thiswin)

    ext = row1['EXTENSION']
    rows = row1['ROW']
    columns = tuple(pipe.infile[ext].get_colnames())
    integ1 = gbtpipe.ConvenientIntegration(pipe.infile[ext][columns][rows])
    
    ext = row2['EXTENSION']
    rows = row2['ROW']
    columns = tuple(pipe.infile[ext].get_colnames())
    integ2 = gbtpipe.ConvenientIntegration(pipe.infile[ext][columns][rows])

    vec1 = integ1.data['DATA']
    vec2 = integ2.data['DATA']
    if 'Vane' in integ1.data['CALPOSITION'][0] and\
            'Observing' in integ2.data['CALPOSITION'][0]:
        onoff = np.median((vec1-vec2)/vec2)
    elif 'Vane' in integ2.data['CALPOSITION'][0] and\
            'Observing' in integ1.data['CALPOSITION'][0]:
        onoff = np.median((vec2-vec1)/vec1)
    else:
        import pdb; pdb.set_trace

    timestamps = integ1.data['DATE-OBS']

    mjd = np.mean(np.array([pipe.pu.dateToMjd(stamp)
                            for stamp in timestamps]))

    elevation = np.mean(integ1.data['ELEVATIO'])

    twarm = np.mean(integ1.data['TWARM']+273.15)
    tbg = 2.725
    tambient = np.mean(integ1.data['TAMBIENT'])
    avgfreq = np.mean(integ1.data['OBSFREQ'])
    tatm = weather.retrieve_Tatm(mjd, avgfreq, log=log, forcecalc=True)
    zenithtau = weather.retrieve_zenith_opacity(mjd,avgfreq,log=log, forcecalc=True)
    
    tau = cal.elevation_adjusted_opacity(zenithtau, elevation)
    tcal = (tatm - tbg) + (twarm - tatm) * np.exp(tau) 
    tsys = tcal / onoff
    return tsys


#cl_params.mapscans = list(np.linspace(113,136,136-113+1).astype('int'))
#cl_params.refscans = [111,]


def calscans(inputdir, start=82, stop=105, refscans = [80], outdir=None):
    if not outdir:
        outdir = os.getcwd()

    log = gbtpipe.Logging()
    w = gbtpipe.Weather()
    cl_params = gbtpipe.initParameters(inputdir)
    cl_params.mapscans = list(np.linspace(start,stop,stop-start+1).astype('int'))
    cl_params.refscans = refscans
    if os.path.isdir(cl_params.infilename):
        log.doMessage('INFO', 'Infile name is a directory')
        input_directory = cl_params.infilename

        # Instantiate a SdFits object for I/O and interpreting the
        #  contents of the index file
        sdf = SdFits()

        # generate a name for the index file based on the name of the
        #  raw SDFITS file.  The index file simply has a different extension
        directory_name = os.path.basename(cl_params.infilename.rstrip('/'))
        indexfile = cl_params.infilename + '/' + directory_name + '.index'
        try:
            # create a structure that lists the raw SDFITS rows for
            #  each scan/window/feed/polarization
            row_list, summary = sdf.parseSdfitsIndex(indexfile, cl_params.mapscans)
        except IOError:
            log.doMessage('ERR', 'Could not open index file', indexfile)
            sys.exit()

    for infilename in glob.glob(input_directory + '/' +
                                os.path.basename(input_directory) +
                                '*.fits'):
            log.doMessage('DBG', 'Attempting to calibrate', 
                          os.path.basename(infilename).rstrip('.fits'))

            # change the infilename in the params structure to the
            #  current infile in the directory for each iteration
            cl_params.infilename = infilename

            # copy the cl_params structure so we can modify it during calibration
            # for each seperate file.
            command_options = copy.deepcopy(cl_params)
            sdf = SdFits()
            cal = gbtpipe.Calibration()
            indexfile = sdf.nameIndexFile(command_options.infilename)
            row_list, summary = sdf.parseSdfitsIndex(indexfile, 
                                                     mapscans=command_options.mapscans)
            feedlist = (row_list.feeds())
            print(feedlist)
            for thisfeed in feedlist:
                thispol = 0
                thiswin = 0
                command_options = copy.deepcopy(cl_params)
                try:
                    pipe = gbtpipe.MappingPipeline(command_options, 
                                                   row_list,
                                                   thisfeed,
                                                   thispol,
                                                   thiswin,
                                                   None)
                except KeyError:
                    pipe = None
                    pass
                if pipe:
                    tsys = gettsys(cl_params, row_list, thisfeed, 
                                   thispol, thiswin, pipe, weather=w, log=log)
                    log.doMessage('INFO', 'Feed: {0}, Tsys (K): {1}'.format(thisfeed,
                                                                            tsys))

                    for thisscan in cl_params.mapscans:
                        print(thisscan)
                        rows = row_list.get(thisscan, thisfeed, thispol, thiswin)

                        ext = rows['EXTENSION']
                        rows = rows['ROW']
                        columns = tuple(pipe.infile[ext].get_colnames())
                        integs = gbtpipe.ConvenientIntegration(\
                            pipe.infile[ext][columns][rows])
                        timestamps = integs.data['DATE-OBS']
                        elevation = np.median(integs.data['ELEVATIO'])
                        mjds = np.array([pipe.pu.dateToMjd(stamp) 
                                         for stamp in timestamps])
                        avgfreq = np.median(integs.data['OBSFREQ'])
                        zenithtau = w.retrieve_zenith_opacity(np.median(mjds), avgfreq,
                                                              log=log, forcecalc=True)
                        tau = cal.elevation_adjusted_opacity(zenithtau, elevation)

                        # Do the calibration here
                        ON = integs.data['DATA']
                        OFF = np.median(ON, axis=0) # Empirical bandpass
                        OFF.shape += (1,)
                        OFF = OFF * np.ones((1, ON.shape[0]))
                        OFF = OFF.T
                        TA = tsys * (ON-OFF)/(OFF)
                        medianTA = np.median(TA, axis=1)
                        medianTA.shape = (1,) + medianTA.shape
                        medianTA = np.ones((ON.shape[1], 1)) * medianTA
                        TA -= medianTA.T * np.exp(tau)
                        
                        for ctr, row in enumerate(rows):
                            row = gbtpipe.Integration(pipe.infile[ext][columns][row])
                            row.data['DATA'] = TA[ctr,:]
                            row.data['TSYS'] = tsys
                            row.data['TUNIT7'] = 'Ta*'
                            pipe.outfile[-1].append(row.data)
                    pipe.outfile.close()
        # gbtpipe.preview_zenith_tau(log, row_list, command_options, 
        #                            row_list.feeds(), row_list.windows(),
        #                            row_list.pols())
        

# sdf = gbtpipe.SdFits()
# indexfile = sdf.nameIndexFile(params.infilename)
# row_list, summary = sdf.parseSdfitsIndex(indexfile, params.mapscans)
# gbtpipe.runPipeline(params)

