# Purpose: Generate moment maps

# Date         Programmer       Description of Changes
#----------------------------------------------------------------------
# 8/3/2020      A.A. Kepley     Original Code


from degas.products import makeMap
import os
import glob 
import aplpy
import matplotlib.pyplot as plt

releaseDir = os.environ['RELEASEDIR']

fitsList = glob.glob(os.path.join(releaseDir,"*.fits"))

for myfits in fitsList:
    # not doing any masking -- mostly for diagnostic purposes.
    makeMap(myfits,maptype='peakIntensity')
    makeMap(myfits,maptype='moment',order=0)


mom0List = glob.glob(os.path.join(releaseDir, "*_mom0.fits"))

for mom0 in mom0List:
    fig = aplpy.FITSFigure(mom0)
    fig.show_colorscale()
    fig.add_beam(color='gray')
    fig.add_colorbar()
    fig.set_title(os.path.basename(mom0))
    fig.save(mom0.replace(".fits",".png"))
    plt.close()

peakList = glob.glob(os.path.join(releaseDir, "*_peakInt.fits"))

for peak in peakList:
    fig = aplpy.FITSFigure(peak)
    fig.show_colorscale()
    fig.add_beam(color='gray')
    fig.add_colorbar()
    fig.set_title(os.path.basename(peak))
    fig.save(peak.replace(".fits",".png"))
    plt.close()
