# This makes decadal pdfs for data that's in single-file format as opposed to decadal-file format.
# This will be combined with plot_decadalpdf_script after completion.
# Also, this makes SST plots

import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np
from glob import glob
import time

import difference_over_time as dot
import plotting_tools as plt_tls
import process_glens_fun as pgf

# Inputs
cntrlFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/control.001.SST.r90x45.shift.annual.nc'
fdbckFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/feedback.001.SST.r90x45.shift.annual.nc'

baselineFlag = 1
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210225_formeeting/'
saveName = 'pdf_SSTglobal_cntrlfdbck_TESTREGION'
dpi_val = 300

levOfInt = 500 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'area' #aspirational

cntrlIntToPlot = [2010,2020,2030,2040,2050,2090]
fdbckIntToPlot = [2020,2030,2040,2050,2090]
timePeriod = 10 #number of years, i.e. 10 = decade

cntrlDset = xr.open_dataset(cntrlFile)
fdbckDset = xr.open_dataset(fdbckFile)

if regionToPlot == 'global':
    cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    dataKey = pgf.discover_data_var(cntrlDset)
    cntrlDarr = cntrlDsetLoi[dataKey]
    cntrlDarrMnSpc = cntrlDarr.mean(dim=['lat','lon'])

    fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    fdbckDarr = fdbckDsetLoi[dataKey]
    fdbckDarrMnSpc = fdbckDarr.mean(dim=['lat','lon'])
else:
    cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    latCondition = (cntrlDsetLoi['lat']>19.5) & (cntrlDsetLoi['lat']<30)
    cntrlDsetLoiLat = cntrlDsetLoi == latCondition
    lonCondition = (cntrlDsetLoi['lon']>-97) & (cntrlDsetLoi['lon']<-80)
    cntrlDarr = cntrlDsetLoiLat == lonCondition
    dataKey = pgf.discover_data_var(cntrlDset)
    # cntrlDsetLoiBetweenLat = cntrlDsetLoi.se(lat>19.5)

    cntrlDarrMnSpc = cntrlDarr.mean(dim=['lat','lon'])
