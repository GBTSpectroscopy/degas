import os
import degas

os.chdir(os.environ['DEGASDIR'])

# getting the observational log file
degas.catalogs.updateLogs()  # technically only needs to be done if it's changed.

# reduce all galaxies
degas.reduceAll(outputDir=os.environ['DEGASDIR'])

# reduce a specific galaxy and overwrite
degas.reduceAll(outputDir=os.environ['DEGASDIR'],galaxyList=['IC0342'],overwriter=True)
