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
