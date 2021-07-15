from degas import pipeline
import os

gallist = ['IC0342',
           'NGC0337',
           'NGC2146',
           'NGC2903',
           'NGC3147',
           'NGC3521',
           'NGC3631',
           'NGC4030',
           'NGC4038',
           'NGC4258',
           'NGC4321',
           'NGC4414',
           'NGC4501',
           'NGC4535',
           'NGC4569',
           'NGC5055',
           'NGC6946'
           'NGC7331',
           'NGC5248']
gallist = gallist[-2:]
# gallist = ['NGC2146', 'IC0342', 'NGC2903']

degasdir = '/mnt/bigdata/erosolow/surveys/DEGAS/'
from degas import catalogs
catalogs.updateLogs('ObservationLog.csv')
pipeline.reduceAll(release='QA0', galaxyList=gallist,
                   OffType='PCA')
