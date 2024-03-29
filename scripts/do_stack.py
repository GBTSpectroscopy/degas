from degas.analysis_stack import makeSampleTable, pruneSampleTable
import os

release = 'IR6p1'

scriptDir = os.environ['SCRIPTDIR']
regridDir = os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')


# stack using mom1
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release)

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='stack_'+release, release='DR1')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_mom1.fits','stack_'+release+'_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

# stack using peak velocity

outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_peakVelocity')

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='peakVelocity',outname='stack_'+release, release='DR1')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_peakVelocity.fits','stack_'+release+'_peakVelocity_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

# stack using sigmaSFR
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_spatialR21')

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='stack_'+release+'_spatialR21', release='DR1', R21='sigmaSFR')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_spatialR21_mom1.fits','stack_'+release+'_spatialR21_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

# Stack using sigmaSFR and peak velocity
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+'_spatialR21_peakVelocity')

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='peakVelocity',outname='stack_'+release+'_spatialR21', release='DR1',R21='sigmaSFR')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_spatialR21_peakVelocity.fits','stack_'+release+'_spatialR21_peakVelocity_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

# stack using 24micron LTIR
outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_'+release+"_L24micron")

myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',ltir='single',outname='stack_'+release+'_L24micron', release='DR1')

myresults2 = pruneSampleTable(outDir,'stack_'+release+'_L24micron_mom1.fits','stack_'+release+'_L24micron_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

#----------------------------------------------------------------------


## stuff for testing.


#from astropy.table import Table
#degas = Table.read(os.path.join(scriptDir,'degas_base.fits'))
#idx = degas['DR1'] == 1
#mylist = degas[idx]['NAME'].tolist()
#mylist.remove('NGC4038')

# testing out removing NGC4038 and it's bad  units for now.
#myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=mylist)

# from degas.analysis_stack import makeSampleTable, pruneSampleTable
# import os
    
# release = 'IR6p1'

# scriptDir = os.environ['SCRIPTDIR']
# regridDir=os.path.join(os.environ['ANALYSISDIR'],release+'_regrid')
# outDir = os.path.join(os.environ['ANALYSISDIR'],'stack_test')

# myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['NGC3521','NGC4569','NGC5055','NGC6946'])

# myresults2 = pruneSampleTable(outDir,'test_mom1.fits','test_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

#sourceList=['NGC2903', 'NGC6946', 'NGC3521']) # no HI, HI, HI

# #myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['NGC2903','NGC3631'])
# #myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='test', release='DR1',sourceList=['IC0342'])

# myresults = makeSampleTable(regridDir, outDir, scriptDir, vtype='mom1',outname='stack_'+release, release='DR1')

# myresults2 = pruneSampleTable(outDir,'stack_'+release+'_mom1.fits','stack_'+release+'_mom1_pruned.fits',overrideFile=os.path.join(scriptDir,'manualOverrides.csv'))

