This is some guidance on how the DEGAS gridding pipeline works.  There are three main steps in the DEGAS calibration.

  1. _Scan Calibration_ -- This takes raw ARGUS data and performs the ON-OFF calibration steps and is run by `gbtpipe.ArgusCal`, specifically the `calscans()` routine.  This routine is managed in a pipeline by the `degas.pipeline.wrapper()` function which parses observation logs and combines those with selection commands to determine the batches of scans to process for a given block of scans on the telescope.  The result of this program is a set of calibrated scans stored in a FITS bintable that track antenna temperature vs frequency with metadata such as RA, DEC and observation time.  The main user-contrallable parameters are:
      * `OffSelector` which how to identify sky spectra within the block of scans with parameters like `SpatialMask` or `RowEnds` as laid out in `gbtpipe.ArgusCal`.
      * `Offtype` -- This determines how to take a set of sky spectra and interpolate an OFF spectrum between them.  This includes `linefit` which fits a line to each the position vs brightness in each spectral channel and uses the line as an OFF or `median` which takes the median off all the sky values.
  2. _Gridding_  -- This step takes a block of spectra with known sky positions and averages them into a map using the On-the-fly mapping procedure.  Each position in the map is a weighted average of the nearby samples, with the weighting determined by the gridding kernel.  This functionality is run by `gbtpipe.Gridding` but the `degas.gridding` routines wrap the `gbtpipe` tools with the selected parameters to create individual data cubes.  The key parameters that control the quality of the data are all in `gbtpipe.Gridding.griddata`:
      * `templateHeader` -- `astropy.fits.Header` object that specifies the astrometry of the output data.
      * `startChannel`, `endChannel` -- Channels to select out of the ARGUS bandpass. This should probably exclude the edges of the band.
      * `doBaseline` -- perform per-scan baseline.  Currenty `True` and it probably should stay that way.
      * `baselineRegion` -- `np.slice` or `list` of slices that specify signal free portions of the spectrum.  Currently set to exclude the edge and a small section in the middle.
      * `blorder` -- Baseline order for fit of scan
      * `flagRMS` -- Boolean (currently `True`) that flags scans that have high RMS values compared to expectations from radiometer equation.  The flagging threshold is set by the `rmsThresh` keyword, which is currently set to 1.25.  RMS is calculated using a robust method that looks at channel-to-channel variations.
      * `flagSpike` -- This looks for channel-to channel jumps in the spectrum that are a factor of `spikeThresh` times the RMS level.  Currently, `spikeThresh` is set to 10.
      * `flagRipple` -- This measures the RMS of a spectrum post-baselineing and flags the spectrum if mean-sum-of-squares of the data is 2x the spectrum noise.  The 2x is hardwired and should probably be experimented with.  Spectra with ripples in them are systematically positive / negative over parts of the bandpass, which inflates the mean-sum-of-squares, prompting flagging.
  3. _Post-Processing_ -- Postprocessing is accomplished by the `degas.postprocess.cleansplit()` routine.  This rebaselines the data, per position (as opposed to per scan) and subtracts of lingering baselines.  The baseline window is specified by velocities of the galaxes in Velocity space, specifically the values `Vgalaxy` (the velocity of the galaxy) and `Vwindow` (the velocity width of the baseline window on each side of `Vgalaxy` so the width is really 2 times `Vwindow`.  The actual routine that calls this is `gbtpipe.Baseline.rebaseline()` which runs on FITS cubes and fits a polynomial with order `blorder`, currently 3.

The main thing to experiment with are (1) baseline order and regions on the per-scan and per-position steps, (2) flagging commands to reject spectra / channels that look bad.  The best option is probably to make calls to the `gbtpipe.griddata` call as below with a fixed filelist and look at the timeseries QA plots to verify functionality of baselining, and flagging.

The current call to the gridding command looks like:

```
gbtpipe.griddata(filelist, # List of files to process
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
