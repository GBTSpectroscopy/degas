
from degas import analysis_stack
import os

dataDir = os.path.join(os.environ['ANALYSISDIR'],'IR5_regrid')



analysis_stack.stack('CO', ##line
                     os.path.join(dataDir,'NGC2903_HCN_rebase3_smooth1.3_hanning1_smooth.fits'), ##cube -- is this already read in via a spectralCube command?
                     vtype, # vtype
                     os.path.join(datadir,'NGC2903_12CO_peakInt_regrid.fits'), ## basemap
                     basetype, ## basetype
                     bintype, ## bintype
                     sfrmap=None, ## sfrmap
                     weightmap= None, ##None
                     database=os.path.join(os.environ['SCRIPTDIR'],'degas_base.fits') #database
                     datadir='') #datadir
                     
                     
