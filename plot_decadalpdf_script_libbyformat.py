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
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210311_regiontesting/'
saveName = 'pdf_SST_cntrlfdbck_TESTPOINT'
dpi_val = 300

levOfInt = 500 #hPa
# latOfInt = np.array([19.5,30])
latOfInt = -28.6
lonOfInt = 112
quantileForFig = 0.66
regionToPlot = 'point' #aspirational

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

elif regionToPlot == 'regional':
    cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt)
    dataKey = pgf.discover_data_var(cntrlDset)
    cntrlDarrLoi = cntrlDsetLoi[dataKey]
    lats = cntrlDsetLoi['lat'] #feedback and control are on same grid, fortunately
    lons = cntrlDsetLoi['lon']
    latMask = (lats>19.5) & (lats<30)
    lonMask = (lons>263) & (lons<280)
    cntrlDarrLoiAoi = cntrlDarrLoi[:,latMask,lonMask]
    cntrlDarrMnSpc = cntrlDarrLoiAoi.mean(dim=['lat','lon'])

    fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    dataKey = pgf.discover_data_var(fdbckDset)
    fdbckDarrLoi = fdbckDsetLoi[dataKey]
    fdbckDarrLoiAoi = fdbckDarrLoi[:,latMask,lonMask]
    fdbckDarrMnSpc = fdbckDarrLoiAoi.mean(dim=['lat','lon'])

elif regionToPlot == 'point':
    cntrlDsetLoiPoi = cntrlDset.sel(z_t=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
    dataKey = pgf.discover_data_var(cntrlDset)
    cntrlDarrLoiPoi = cntrlDsetLoiPoi[dataKey]
    cntrlDarrMnSpc = cntrlDarrLoiPoi #need better variable name

    fdbckDsetLoiPoi = fdbckDset.sel(z_t=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
    dataKey = pgf.discover_data_var(fdbckDset)
    fdbckDarrLoiPoi = fdbckDsetLoiPoi[dataKey]
    fdbckDarrMnSpc = fdbckDarrLoiPoi

else:
    print("Invalid region!")


# take 2010-2019 average and remove this
# baselineMeanToRmv = dot.average_over_years(cntrlDarrMnSpc,2010,2019)
# cntrlDarrMnSpcNorm = cntrlDarrMnSpc - baselineMeanToRmv
# fdbckDarrMnSpcNorm = fdbckDarrMnSpc - baselineMeanToRmv

cntrlActive = cntrlDarrMnSpc
fdbckActive = fdbckDarrMnSpc

cntrlYears = cntrlActive['time'].dt.year.data
cntrlHandlesToPlot = list()
cntrlHandlesToPlot = pgf.extract_doi(cntrlIntToPlot, cntrlYears, timePeriod, cntrlActive, cntrlHandlesToPlot)

fdbckYears = fdbckActive['time'].dt.year.data
fdbckHandlesToPlot = list()
fdbckHandlesToPlot = pgf.extract_doi(fdbckIntToPlot, fdbckYears, timePeriod, fdbckActive, fdbckHandlesToPlot)

handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot

colorsToPlot = plt_tls.select_colors(baselineFlag,len(cntrlIntToPlot)-1,len(fdbckIntToPlot))
labelsToPlot = list(['2010-2019 Baseline'])
labelsToPlot = plt_tls.generate_labels(labelsToPlot, cntrlIntToPlot, timePeriod, 'RCP8.5')
labelsToPlot = plt_tls.generate_labels(labelsToPlot, fdbckIntToPlot, timePeriod, 'SAI')
labelsToPlot.append('(-28.6,112) SST PDFs in GLENS')

print(colorsToPlot)
plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, savePath, saveName)

print('Completed!')
