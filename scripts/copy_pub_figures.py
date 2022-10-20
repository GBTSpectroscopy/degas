import os
import shutil
import glob

baseDir = '/lustre/cv/users/akepley/degas'
regridDir = os.path.join(baseDir,'IR6p1_regrid')
mom0Dir = os.path.join(baseDir,'') # need to define this once I get the moment 0 images I want.
stackDir = os.path.join(baseDir,'stack_IR6p1')

paperDir = os.path.join(baseDir,'plots')

# copy sample overview plots -- DO THIS NEED TO HAPPEN?

# copy stack 


# copy stacks plots
filelist = glob.glob(os.path.join(stackDir,'stack_plots_pruned',"*_radius_stacks.pdf"))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

# copy stack trends

# copy empire degas comparison
filelist = glob.glob(os.path.join(baseDir,'empire_degas_comp','*_HCN_stack_compare.pdf'))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

filelist = glob.glob(os.path.join(baseDir,'empire_degas_comp','*_HCOp_stack_compare.pdf'))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

filelist = glob.glob(os.path.join(baseDir,'empire_degas_comp','*_13CO_stack_compare.pdf'))
for myfile in filelist:
    shutil.copy(myfile, paperDir)



# copy fdense and sfedense plots
filelist = glob.glob(os.path.join(baseDir,'sfe_fdense_trends','fdense_*.pdf'))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

filelist = glob.glob(os.path.join(baseDir,'sfe_fdense_trends','sfe_dense_*.pdf'))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

# copy fdense and sfedense fit coeff plots
filelist = glob.glob(os.path.join(baseDir,'sfe_fdense_trends','ratio_*hist*.pdf'))
for myfile in filelist:
    shutil.copy(myfile,paperDir)


# copy correlation plots
filelist = glob.glob(os.path.join(baseDir,"sfe_fdense_trends","*correlation_coeffs.pdf"))
for myfile in filelist:
    shutil.copy(myfile, paperDir)

filelist = glob.glob(os.path.join(baseDir,"sfe_fdense_trends","*correlation_coeffs_spatialR21.pdf"))
for myfile in filelist:
    shutil.copy(myfile, paperDir)
