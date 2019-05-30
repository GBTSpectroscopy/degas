This is some guidance on how the DEGAS gridding pipeline works.  There are three main steps in the DEGAS calibration.

  1. _Scan Calibration_ -- This takes raw ARGUS data and performs the ON-OFF calibration steps and is run by `gbtpipe.ArgusCal`, specifically the `calscans()` routine.  This routine is managed in a pipeline by the `degas.pipeline.wrapper()` function which parses observation logs and combines those with selection commands to determine the batches of scans to process for a given block of scans on the telescope.  The result of this program is a set of calibrated scans stored in a FITS bintable that track antenna temperature vs frequency with metadata such as RA, DEC and observation time.  The main user-contrallable parameters are:
      * `OffSelector` which how to identify sky spectra within the block of scans with parameters like `SpatialMask` or `RowEnds` as laid out in `gbtpipe.ArgusCal`.
      * `Offtype` -- This determines how to take a set of sky spectra and interpolate an OFF spectrum between them.  This includes `linefit` which fits a line to each the position vs brightness in each spectral channel and uses the line as an OFF or `median` which takes the median off all the sky values.
  2. _Gridding_  -- This step takes a block of spectra with known sky positions and averages them into a map using the On-the-fly mapping procedure.  Each position in the map is a weighted average of the nearby samples, with the weighting determined by the gridding kernel.  This functionality is run by `gbtpipe.Gridding` but the `degas.gridding` routines wrap the `gbtpipe` tools with the selected parameters to create individual data cubes.  The key parameters that control the quality of the data are all in `gbtpipe.Gridding.griddata`:
      * `templateHeader` -- `astropy.fits.Header` object that specifies the astrometry of the output data.
      * `startChannel`, `endChannel` -- Channels to select out of the ARGUS bandpass. This should probably exclude the edges of the band.
      * 
  3. _Post-Processing_ -- 
  
The current call to the gridding command looks like:

```
gbtpipe.griddata(filelist,
                 startChannel=edgetrim,  # edgetrim default to 100 channels
                 endChannel=1024-edgetrim, 
                 baselineRegion = [slice(edgetrim,
                                         edgetrim+basebuff,1),  # Basebuff set to 64
                                   slice(448,576,1),
                                   slice(1024-edgetrim-basebuff,
                                         1024-basebuff,1)],  # Region to fit a baseline PER SCAN in channel units
                 outdir=OutputDirectory,  # Where to stick the file
                 blorder=5,    # Scans are fit by 5th order polynomial
                 flagRMS=True,  # Flag scans with RMS higher than expected noise level from radiometer equation times a threshold
                 plotTimeSeries=plotTimeSeries,  # Generate QA plots
                 flagRipple=True, # Flag ripples in the spectrum
                 pixPerBeam=4.0,   # Pixels per beam
                 plotsubdir='timeseries',
                 outname=filename, **kwargs)
```
