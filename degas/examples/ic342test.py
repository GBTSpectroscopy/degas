from degas import pipeline
import os
from degas.gridding import gridGalaxy
from degas import catalogs

ppo = False
gallist = ['IC0342']
degasdir = '/mnt/bigdata/erosolow/surveys/DEGAS/'
datadir='/mnt/bigdata/erosolow/surveys/DEGAS/'

#catalogs.updateLogs('ObservationLog.csv')
#pipeline.reduceAll(release='QA0', galaxyList=gallist,
#                   OffType='PCA')

gridGalaxy(galaxy='IC0342', setup='12CO',
           release='QA0', datadir=datadir,
           PostprocOnly=ppo)

gridGalaxy(galaxy='IC0342', setup='HCN_HCO+',
           release='QA0', datadir=datadir, PostprocOnly=ppo)


gridGalaxy(galaxy='IC0342', setup='13CO_C18O',
           release='QA0', datadir=datadir, PostprocOnly=ppo)



HCNgals = [
    'NGC4038',
    'NGC2146',
    'NGC6946',
    'NGC7331',
    'NGC5248',
    'NGC2903',
    'NGC4321',
    'NGC5055',
    'NGC4501',
    'NGC3147',
    'NGC3521',
    'NGC4414',
    'NGC0337',
    'NGC3631',
    'NGC4030',
    'NGC4258',
    'NGC4535',
    'NGC4569',
]

HCNgals=['IC0342']

COgals = [
    'NGC4038',
    'NGC2146',
    'NGC7331',
    'NGC2903',
    'NGC4321',
    'NGC5055',
    'NGC4501',
    'NGC3147',
    'NGC0337',
    'NGC4569',
    'NGC3521',
    'NGC3631',
    'NGC4030',
    'NGC4258',
    'NGC4414',
    'NGC4535',
    'IC0342',
]



# gridGalaxy(galaxy='NGC5055', setup='13CO_C18O', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC5055', setup='HCN_HCO+', release='QA0', datadir=datadir)



#gridGalaxy(galaxy='NGC7331', setup='HCN_HCO+',
#           release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC6946', setup='HCN_HCO+', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC4569', setup='HCN_HCO+', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC4569', setup='13CO_C18O', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC4501', setup='13CO_C18O', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC4501', setup='HCN_HCO+', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC4414', setup='13CO_C18O', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC4414', setup='HCN_HCO+', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC4321', setup='HCN_HCO+', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC4321', setup='13CO_C18O', release='QA0', datadir=datadir)

# # gridGalaxy(galaxy='NGC4038', setup='13CO_C18O', release='QA0', datadir=datadir)
# # gridGalaxy(galaxy='NGC4038', setup='HCN_HCO+', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC3521', setup='HCN_HCO+', release='QA0', datadir=datadir)


# gridGalaxy(galaxy='NGC2903', setup='HCN_HCO+', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC2903', setup='13CO_C18O', release='QA0', datadir=datadir)

# gridGalaxy(galaxy='NGC2146', setup='HCN_HCO+', release='QA0', datadir=datadir)
# gridGalaxy(galaxy='NGC2146', setup='13CO_C18O', release='QA0', datadir=datadir)

# # gridGalaxy(galaxy='IC0342', setup='13CO_C18O', release='QA0', datadir=datadir)
# # gridGalaxy(galaxy='IC0342', setup='HCN_HCO+', release='QA0', datadir=datadir)

