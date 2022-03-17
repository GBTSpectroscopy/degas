from degas.analysis_stack import makeSampleTable, pruneSampleTable
import os

release = 'IR6p1'
#release = 'IR6p1'

scriptDir = os.environ['SCRIPTDIR']
regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid_33as')
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+"_33as")
#outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='stack_'+release+"_33as", release='DR1')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_33as_mom1.fits','stack_'+release+'_33as_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))



