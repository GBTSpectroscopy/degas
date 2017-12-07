import gbtpipe
import glob
import os
from . import catalogs

def gridGalaxy(galaxy='IC0342', setup='13CO_C18O',
               datadir='/lustre/pipeline/scratch/DEGAS/',
               overwrite=True, release='QA2', edgetrim = 100,
               basebuff = 64):

    # Note that we also use a few channels in the middle.

    filelist = glob.glob(datadir + galaxy + '/' +
                         setup + '/*fits')
    import pdb; pdb.set_trace()
    OutputDirectory = datadir + galaxy + '/images/'
    if not os.access(OutputDirectory, os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            os.chdir(OutputDirectory)
        except OSError:
            raise
    if '12CO' in setup:    
        gbtpipe.griddata(filelist,
                         startChannel=edgetrim,
                         endChannel=1024-edgetrim,
                         baselineRegion = [slice(edgetrim,
                                                 edgetrim+basebuff,1),
                                           slice(1024-edgetrim-basebuff,
                                                 1024-basebuff,1)],
                         flagRMS=True, plotTimeSeries=True,
                         flagRipple=True, pixPerBeam=4.0)

    else:
        gbtpipe.griddata(filelist,
                         startChannel=edgetrim,
                         endChannel=1024-edgetrim,
                         baselineRegion = [slice(edgetrim,
                                                 edgetrim+basebuff,1),
                                           slice(448,576,1),
                                           slice(1024-edgetrim-basebuff,
                                                 1024-basebuff,1)],
                         flagRMS=True, plotTimeSeries=True,
                         flagRipple=True, pixPerBeam=4.0)


    
