This package builds data cubes for ARGUS observations of dense gas tracers.

## Installation

The package expects and Anaconda like user-based python environment which you can install on many machines from [continuum.io](https://www.continuum.io/downloads).  In addition to many of the standard astronomical packages, this software currently requires a few beta packages for GBT-specific context.

First, install the GBT pipeline package:
```
pip install https://github.com/low-sky/gbtpipe/archive/master.zip
```

Next, install the Keystone Survey's gridder utilities:
```
pip install https://github.com/GBTAmmoniaSurvey/keystone/archive/master.zip
```

Finally, install the `degas` package:
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

