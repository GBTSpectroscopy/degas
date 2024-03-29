import gbtpipe
import glob
import os
import pkg_resources
import warnings
from . import catalogs
from . import postprocess


def gridBlocks(release='QA0', edgetrim = 100,
               basebuff = 64,
               datadir='/lustre/pipeline/scratch/DEGAS/',
               setup='13CO_C18O'):

    """
    This builds tiny maps suitable for QA1 of blocks
    """

    pipeversion = pkg_resources.get_distribution("degas").version
    filelist = glob.glob(datadir + galaxy + '/' +
                         setup + '/*fits')

def gridGalaxy(galaxy='IC0342', setup='13CO_C18O',
               datadir='/lustre/pipeline/scratch/DEGAS/',
               overwrite=True, release='QA1', edgetrim = 100,
               basebuff = 64, plotTimeSeries=True, PostprocOnly=False,
               scanblorder=11, posblorder=7, pipeversion=None,
               **kwargs):

    # Override automatic version; needed for development sometimes.
    if pipeversion is None:
        pipeversion = pkg_resources.get_distribution("degas").version

    setup_dict = {'13CO_C18O':'13co_c18o',
                  'HCN_HCO+':'hcn_hcop',
                  '12CO':'12co'}

    filelist = glob.glob(os.path.join(datadir,galaxy,setup,"*.fits"))
    OutputDirectory = os.path.join(datadir,galaxy,"images")

    if not os.access(OutputDirectory, os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            os.chdir(OutputDirectory)
        except OSError:
            raise
    else:
        os.chdir(OutputDirectory)
    if '12CO' in setup:
        filename = galaxy + '_' + setup + '_v{0}'.format(pipeversion)
        maskfile= (datadir + '/masks/BOOLEAN_MASKS/' + galaxy
                   + '.' + setup_dict[setup]
                   + '.mask.fits')

        if not PostprocOnly:
            gbtpipe.griddata(filelist,
                             startChannel=edgetrim,
                             endChannel=1024-edgetrim,
                             outdir=OutputDirectory,
                             flagSpike=True, spikeThresh=3,
                             flagRMS=True,  plotTimeSeries=plotTimeSeries,
                             flagRipple=True, rippleThresh=1.2,
                             pixPerBeam=4.0,
                             rmsThresh=1.1,
                             robust=False,
                             blorder=scanblorder,
                             plotsubdir='timeseries/',
                             windowStrategy='cubemask',
                             maskfile=maskfile,
                             outname=filename, **kwargs)

        postprocess.cleansplit(filename + '.fits',
                               spectralSetup=setup,
                               HanningLoops=1,
                               spatialSmooth=1.3,
                               maskfile=maskfile,
                               blorder=posblorder,
                               **kwargs)
    else:
        filename = galaxy + '_' + setup + '_v{0}'.format(pipeversion)
        maskfile= (datadir + '/masks/BOOLEAN_MASKS/' + galaxy
                   + '.' + setup_dict[setup]
                   + '.mask.fits')
        
        if not PostprocOnly:
            if os.path.exists(maskfile):
                warnings.warn("Processing with mask {0}".format(maskfile))

                gbtpipe.griddata(filelist,
                                 startChannel=edgetrim,
                                 endChannel=1024-edgetrim,
                                 outdir=OutputDirectory,
                                 flagSpike=True, spikeThresh=3,
                                 flagRMS=True,  plotTimeSeries=plotTimeSeries,
                                 flagRipple=True, rippleThresh=1.2,
                                 pixPerBeam=4.0,
                                 rmsThresh=1.1,
                                 robust=False,
                                 blorder=scanblorder,
                                 plotsubdir='timeseries/',
                                 windowStrategy='cubemask',
                                 maskfile=maskfile,
                                 outname=filename, **kwargs)
                
            else:
                warnings.warn("Mask file not found. Proceeding to grid without mask file")

                gbtpipe.griddata(filelist,
                                 startChannel=edgetrim,
                                 endChannel=1024-edgetrim,
                                 outdir=OutputDirectory,
                                 flagSpike=True, spikeThresh=3,
                                 flagRMS=True,  plotTimeSeries=plotTimeSeries,
                                 flagRipple=True, rippleThresh=1.2,
                                 pixPerBeam=4.0,
                                 rmsThresh=1.1,
                                 robust=False,
                                 blorder=scanblorder,
                                 plotsubdir='timeseries/',
                                 outname=filename, **kwargs)
        if os.path.exists(maskfile):
            warnings.warn("Processing with mask {0}".format(maskfile))
            postprocess.cleansplit(filename + '.fits',
                                   spectralSetup=setup,
                                   HanningLoops=1, blorder=posblorder,
                                   spatialSmooth=1.3,
                                   maskfile=maskfile,
                                   **kwargs)
        else:
            warnings.warn("Processing without mask")
            postprocess.cleansplit(filename + '.fits',
                                   spectralSetup=setup,
                                   HanningLoops=1, blorder=posblorder,
                                   spatialSmooth=1.3,
                                   **kwargs)
