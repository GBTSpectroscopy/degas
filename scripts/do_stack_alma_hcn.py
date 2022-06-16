from degas.analysis_stack import makeSampleTable, pruneSampleTable
import os

release = 'alma_hcn'

scriptDir = os.environ['SCRIPTDIR']
regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)
#outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='stack_'+release, release='DR1',sourceList=['NGC2903','NGC4321'])

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_mom1.fits','stack_'+release+'_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))



