# Purpose: Generate moment maps

# Date         Programmer       Description of Changes
#----------------------------------------------------------------------
# 8/3/2020      A.A. Kepley     Original Code
# 7/29/2021     A.A. Kepley     Updates for IR6p0


from degas.products import makeMap
import os
import glob 
import aplpy
import matplotlib.pyplot as plt

releaseDir = os.path.join(os.environ['ANALYSISDIR'],'IR6p1')

fitsList = glob.glob(os.path.join(releaseDir,"*hanning1.fits"))

mom0Dir = os.path.join(releaseDir,'mom0')
if not os.path.exists(mom0Dir):
    os.mkdir(mom0Dir)

peakDir = os.path.join(releaseDir,'peak')
if not os.path.exists(peakDir):
    os.mkdir(peakDir)

for myfits in fitsList:
    # not doing any masking -- mostly for diagnostic purposes.
    makeMap(myfits,peakDir,maptype='peakIntensity')
    makeMap(myfits,mom0Dir, maptype='moment',order=0)


mom0List = glob.glob(os.path.join(mom0Dir, "*_mom0.fits"))

for mom0 in mom0List:
    fig = aplpy.FITSFigure(mom0)
    fig.show_colorscale()
    fig.add_beam(color='gray')
    fig.add_colorbar()
    fig.set_title(os.path.basename(mom0))
    fig.save(mom0.replace(".fits",".png"))
    plt.close()

peakList = glob.glob(os.path.join(peakDir, "*_peakInt.fits"))

for peak in peakList:
    fig = aplpy.FITSFigure(peak)
    fig.show_colorscale()
    fig.add_beam(color='gray')
    fig.add_colorbar()
    fig.set_title(os.path.basename(peak))
    fig.save(peak.replace(".fits",".png"))
    plt.close()
