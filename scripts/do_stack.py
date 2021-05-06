from degas.analysis_stack import makeSampleTable
import os

release = 'IR5p1'

scriptDir = os.environ['SCRIPTDIR']
regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

#myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['NGC2903'])

results = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1')
