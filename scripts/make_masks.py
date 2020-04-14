import os
from degas.masking import cubemask
#import degas
from astropy.table import Table
import glob

# set desired mask parameters
peak_cut = 5.0
low_cut = 3.0

# set up the relevant directories
maskDir = os.environ['MASKDIR'] 
otherDataDir = os.environ['OTHERDATA']
scriptDir = os.environ['SCRIPTDIR']

# get list of galaxies in degas DR1
degas_table = Table.read(os.path.join(scriptDir,"degas_base.fits"))

idx_dr1 = degas_table['DR1'] == 1

# Extract list of galaxies via fancy list comprehension

# heracles
heracles_list =  [os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'heracles','*gauss15.fits'))]

# bima song
bima_list = [ os.path.basename(image).split('_')[0] for image in glob.glob(os.path.join(otherDataDir,'bima_song','*gauss15_fixed.fits'))]

# OVRO
other_list1 = [ os.path.basename(image).split('.')[0].upper() for image in glob.glob(os.path.join(otherDataDir,'temp_co','*co.cmmsk.fits'))]

# NRO
other_list2 = [ os.path.basename(image)[0:5].replace('N','NGC') for image in glob.glob(os.path.join(otherDataDir,'temp_co','N*ICO.FITS'))]

## then there's two misc( NGC4030, ngc4038), which are one-offs.

for galaxy in degas_table[idx_dr1]:

    if galaxy['NAME'] in heracles_list:
        cubemask(os.path.join(otherDataDir,'heracles',
                              galaxy['NAME']+'_heracles_gauss15.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)
    elif galaxy['NAME'] in bima_list:
        cubemask(os.path.join(otherDataDir,'bima_song',
                              galaxy['NAME']+'_bima_gauss15_fixed.fits'),
                 galaxy['NAME']+'_mask.fits',
                 outDir=maskDir,
                 peak_cut=peak_cut,low_cut=low_cut)
    elif galaxy['NAME'] in other_list1:
        ## cmmsk -- OVRO -- original images need to be smoothed
        print("in otherlist1")
    elif galaxy['NAME'] in other_list2:
        ## check on smoothing here
        print("in otherlist2")
    elif galaxy['NAME'] == 'NGC4030':
        ## check on smoothing here
        print("in assorted")
    elif galaxy['NAME'] == 'NGC4038':
        ## check on smoothing here
        print("in assorted")
    else:
        print(galaxy['NAME']+" doesn't appear to have ancillary CO data.")


    

#cubemask('/lustre/cv/users/akepley/degas/ancillary_data/heracles/NGC0337_heracles_gauss15.fits','test.fits',peak_cut=5.0,low_cut=3.0,outDir=maskDir)


