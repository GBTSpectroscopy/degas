from degas.gridding import gridGalaxy

datadir='/mnt/bigdata/erosolow/surveys/DEGAS/'

datadir = '/mnt/bigdata/erosolow/surveys/DEGAS/'
ppo = True
# gridGalaxy(galaxy='IC0342', setup='12CO',
#            release='QA0', datadir=datadir,
#            PostprocOnly=ppo)

gallist = [#'IC0342',
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
           'NGC6946']
gallist = gallist[-2:]
HCNgals = gallist

# HCNgals = ['NGC2903', 'NGC2146', 'IC0342']
# HCNgals = ['IC0342']
for gal in HCNgals:
    gridGalaxy(galaxy=gal, setup='HCN_HCO+',
               release='QA0', datadir=datadir, PostprocOnly=ppo)


COgals = gallist
# COgals = ['NGC2903', 'NGC2146', 'IC0342']
# COgals = ['NGC2146']
# COgals = ['IC0342']
for gal in COgals:
     gridGalaxy(galaxy=gal, setup='13CO_C18O',
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

