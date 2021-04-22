from degas.analysis_stack import makeResultsFITSTable
import os

release = 'IR5p1'

scriptDir = os.environ['SCRIPTDIR']
regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

#myresults = makeResultsFITSTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['NGC2146'])

results = makeResultsFITSTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1')
