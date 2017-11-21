from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models, fitting
import gbtpipe
import degas
import numpy as np
from gbtpipe import Weather,Pipeutils, Calibration

# Adapted from Amanda Kepley's peak_fit.py.

def peak_fit(dcrfile, vaneCalScan=11, peakScan=13, feedNum=10):

# Input
# -------

# dcrfile = 'AGBT17A_304_02.raw.dcr.fits'
# vanecalScan = 11
#peakScan = 13
# feedNum = 10 # Feed 10 is the feed that tyypically points. Note this is the 1-based feed numbering, not the zero-based GBT numbering.

pu = Pipeutils()
weather = Weather()
cal = Calibration()

# open the dcr file and get the data
# ----------------------------------

peak = fits.open(dcrfile)
dcr = peak[1]

# Calculate Tcal from the vanecal scan
# -------------------------------------

# This is partially based on the degas code to calculate tsys, but
# with all the fancy stuff thrown out so I can experiment.

# hack to get weather calcs to work. needs to have log directory already existing.
log = gbtpipe.Logging('', 'gbtpipeline')

# doing the actual calculation
twarm = np.mean(dcr.data['TWARM'][(dcr.data['SCAN'] == vanecalScan) & (dcr.data['FEED'] == feedNum)]+273.15)
tbg = 2.725 #K

timestamps = dcr.data['DATE-OBS'][(dcr.data['SCAN'] == vanecalScan) & (dcr.data['FEED'] == feedNum)]
mjd = np.mean(np.array([pu.dateToMjd(stamp) for stamp in timestamps]))
avgfreq = np.mean(dcr.data['OBSFREQ'][(dcr.data['SCAN'] == vanecalScan) & (dcr.data['FEED'] == feedNum)])

tatm = weather.retrieve_Tatm(mjd, avgfreq,  log=log, forcecalc=True)
zenithtau = weather.retrieve_zenith_opacity(mjd,avgfreq,log=log, forcecalc=True)
elevation = np.mean(dcr.data['ELEVATIO'][(dcr.data['SCAN'] == vanecalScan) & (dcr.data['FEED'] == feedNum)])

tau = cal.elevation_adjusted_opacity(zenithtau, elevation)
tcal = (tatm - tbg) + (twarm - tatm) * np.exp(tau) 
vaneCounts = np.mean(dcr.data['DATA'][(dcr.data['SCAN'] == vanecalScan) & (dcr.data['FEED'] == feedNum)])
skyCounts = np.mean(dcr.data['DATA'][(dcr.data['SCAN'] == vanecalScan+1) & (dcr.data['FEED'] == feedNum)])

# Get the Peak scans
# ------------------

# peak scans are in groups of four: 2 azimuth scans and then 2
# elevation scans. 
peak1counts = dcr.data['DATA'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]

peak1ra = dcr.data['CRVAL2'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]
peak1dec = dcr.data['CRVAL3'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]

peak1ra0 = dcr.data['TRGTLONG'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]
peak1dec0 = dcr.data['TRGTLAT'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]

peak1coord =  dcr.data['RADESYS'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]

# Calculating the offsets.
# ------------------------

peak1path = SkyCoord(ra=peak1ra,dec=peak1dec,frame=peak1coord[0].lower(),unit='deg')
peak1center =  SkyCoord(ra=peak1ra0,dec=peak1dec0,frame=peak1coord[0].lower(),unit='deg')
peak1offset = peak1center.separation(peak1path)
minidx = peak1offset.arcmin.argmin()

peak1offset[0:minidx] = -1.0*peak1offset[0:minidx]
plt.plot(peak1offset.arcmin,peak1counts)

# Applying Calculate Tsys for the Calibration Scans
# --------------------------------------------------

## TODO -- how is this different from the spectral line case?? I think
## it's just that we have no off, so we're scaling absolutely. I've
## got a guess below.

avgfreq = np.mean(dcr.data['OBSFREQ'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)])

zenithtau = weather.retrieve_zenith_opacity(mjd,avgfreq,log=log, forcecalc=True)
elevation = np.mean(dcr.data['ELEVATIO'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)])
tau = cal.elevation_adjusted_opacity(zenithtau, elevation)

peak1Ta = (tcal/(vaneCounts-skyCounts)) * peak1counts  
plt.plot(peak1offset.arcmin,peak1Ta)

# Fitting a model to the peak
# -----------------------------

#g_init = models.Gaussian1D(amplitude=2000,mean=0.0,stddev=0.2) + models.Polynomial1D(1)
g_init = models.Gaussian1D(amplitude=(max(peak1Ta) - median(peak1Ta)),mean=0.0,stddev=10.0/60.0) + models.Polynomial1D(2,c0=0)
fit_g = fitting.LevMarLSQFitter()

# Fitting with a large offset in the polynominal causes havoc, so
# subtract a baseline from the data. This is a trick Carl uses a lot.
rescaleFactor = median(peak1Ta)
g = fit_g(g_init,peak1offset.arcmin,peak1Ta - rescaleFactor)
plt.plot(peak1offset.arcmin,peak1Ta)
plt.plot(peak1offset.arcmin,g(peak1offset.arcmin)+rescaleFactor)

# Get the source flux density from the ALMA data base
# ---------------------------------------------------

srcName = dcr.data['OBJECT'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)][0]
timestamps = dcr.data['DATE-OBS'][(dcr.data['SCAN'] == peakScan) & (dcr.data['FEED'] == feedNum)]

# Code between --- from Brian Kent based on beta code in ALMA pipeline. Used with permission. 
#---
from caldb import fluxservice
from dateutil.parser import parse

scanDate = parse(timestamps[0])

scanDateStr = scanDate.strftime("%d-%B-%Y")

table = fluxservice(frequency=str(avgfreq), sourcename="J"+srcName, date=scanDateStr)

Snu = float(table['fluxdensity'])
print "Snu=",Snu
#---

# Calculating the efficiencies
# ----------------------------

# using equations from Ron Maddelena's memo and David Frayer's guide to Argus

d = 100e2 #m

#eta_a = 0.352 * g.amplitude_0 * np.exp(tau) / Snu

# Since it's Ta* that comes out of the above equations, I think I can drop the np.exp(tau) factor.
eta_a = 0.352 * g.amplitude_0 / Snu
print 'eta_a', eta_a

avgwave = (avgfreq * u.Hz).to(u.cm,equivalencies=u.spectral())

# Here I'm assuming the x and y directions are the
# same. Straightforward to do with other direction of peak scan.
fwhm = 2.0 * sqrt(2.0*math.log(2.0))*g.stddev_0 /60.0 # degrees

eta_mb = (0.8899 * (math.radians(fwhm))**2 * d**2 / avgwave.value**2) * eta_a
print 'eta_mb', eta_mb
