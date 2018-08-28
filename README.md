This package builds data cubes for ARGUS observations of dense gas tracers.

## Installation

The package expects an Anaconda-like user-based python environment which you can install on many machines from [continuum.io](https://www.continuum.io/downloads).  In addition to many of the standard astronomical packages, this software currently requires a few beta packages for GBT-specific context.

Install the `degas` package:
```
pip install https://github.com/low-sky/degas/archive/master.zip
```


## Running the Pipeline

To execute the pipeline, we run the `calscans` data to calibrate the scans.  In the example file, the data are in the directory `TGBT15A_901_34.raw.vegas` running from scan 82 to scan 105.  The reference scans refer to the first scan of a vane calibration pair.

```
import degas
degas.calscans('TGBT15A_901_34.raw.vegas', start=82, stop=105, refscans=[80])
```
The pipeline creates a set of calibrated FITS files which can then be gridded into cubes using the gridding example given in the `e2e.py` file in the `examples` directory.

## Local Installation of the DEGAS pipeline

You can run the DEGAS calibration pipeline on your own computer instead of on a GB computer. You will need:
* The DEGAS data tree
* The GBT Weather database
* A local copy of this repository

You can organize this however you want. A suggested structure might be a main degas directory with the DEGAS rawdata, GBT Weather database, and code underneath it.

_Note: if you are copying between NRAO sites, be careful! The bwlimit given below is greater than the bandwidth of the internal links. Contact Amanda before trying to do this._

### DEGAS Data tree

Copy over the DEGAS data tree to your local computer.  In practice, this is easiest by doing a push from the GB machines to your computer using `rsync`.  Specifically, you need to grab the `rawdata` directory from our space on lustre. If you can't see lustre, make sure you are on one of the computers connected to lustre (list here: http://www.gb.nrao.edu/pubcomputing/public.shtml. You generally want the ones at the bottom (newton, planck, etc)).

```
rsync -ahv --bwlimit=8000 /lustre/pipeline/scratch/DEGAS/rawdata username@machine:/path/to/destination/DEGAS/.
```


### Weather Database

Copy over the GBT weather database to your machine.  It is located at `/users/rmaddale/Weather/ArchiveCoeffs/` on the GB machines.  You only need the files that start with `Coeffs`.  Make a directory on your local machine.  Let's call it `/home/username/GBTWeather/` in this example. 

```
rsync -ahv --bwlimit=8000 /users/rmaddale/Weather/ArchiveCoeffs/Coeffs* username@machine:/home/username/GBTWeather
```

Again, if you can't see lustre, make sure you are on one of the computers connected to lustre (list here: http://www.gb.nrao.edu/pubcomputing/public.shtml. You generally want the ones at the bottom (newton, planck, etc)).


### Install and configure the DEGAS repository

1. Clone this repository and install it locally by running
`python setup.py install`
from inside the `degas` directory. Bonus hint: To find the path for the repository, click on the green button in the top right hand side of the page that says "clone or download". Type "git clone" on the command line and copy and paste the location displayed when you click the green button to the command line.

2. Set environment variables that point to the DEGAS data and the GBT Weather.  The variable DEGASDIR should be pointed to the parent directory containing the raw file.  The variable GBTWEATHER points to the GBT Weather archive.

```
export DEGASDIR=/path/to/destination/DEGAS/
export GBTWEATHER=/home/username/GBTWeather/
```

3. Try it out.  You can check that things are working by trying the calibration pipeline. From `ipython` prompt, you can try:

```
import os
import degas
os.chdir(os.environ['DEGASDIR'])
degas.catalogs.updateLogs()
degas.reduceAll(outputDir=os.environ['DEGASDIR'])
```
