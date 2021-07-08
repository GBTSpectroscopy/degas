from degas.analysis_stack import makeSampleTable
import os

release = 'IR5p1'

scriptDir = os.environ['SCRIPTDIR']
regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

#from astropy.table import Table
#degas = Table.read(os.path.join(scriptDir,'degas_base.fits'))
#idx = degas['DR1'] == 1
#mylist = degas[idx]['NAME'].tolist()
#mylist.remove('NGC4038')

# testing out removing NGC4038 and it's bad  units for now.
#myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=mylist)


#myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['NGC0337'])

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1')
