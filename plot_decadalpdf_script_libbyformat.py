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
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210301_regiontesting/'
saveName = 'pdf_SSTglobal_cntrlfdbck_TESTREGION'
dpi_val = 300

levOfInt = 500 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

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
    print("Not TODAY")
    # cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    # dataKey = pgf.discover_data_var(cntrlDset)
    # cntrlDsetLoi = cntrlDsetLoi[dataKey]
    # print(cntrlDsetLoi)
    # latCondition = np.where((cntrlDsetLoi['lat']>19.5) & (cntrlDsetLoi['lat']<30))[0]
    # print(latCondition)
    # cntrlDsetLoiLat = cntrlDsetLoi[latCondition]
    # lonCondition = np.where((cntrlDsetLoi['lon']>-97) & (cntrlDsetLoi['lon']<-80))[0]
    # cntrlDarr = cntrlDsetLoiLat[lonCondition]
    # cntrlDarrMnSpc = cntrlDarr.mean(dim=['lat','lon'])
    #
    # fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
    # dataKey = pgf.discover_data_var(fdbckDset)
    # fdbckDsetLoi = fdbckDsetLoi[dataKey]
    # print(fdbckDsetLoi)
    # # latCondition2 = np.where((fdbckDsetLoi['lat']>19.5) & (fdbckDsetLoi['lat']<30))[0]
    # # print(latCondition2)
    # import sys
    # sys.exit()
    # fdbckDsetLoiLat = fdbckDsetLoi[latCondition]
    # # lonCondition = np.where((fdbckDsetLoi['lon']>-97) & (fdbckDsetLoi['lon']<-80))
    # fdbckDarr = fdbckDsetLoiLat[lonCondition]
    # fdbckDarrMnSpc = fdbckDarr.mean(dim=['lat','lon'])
    # print(fdbckDarrMnSpc)
# print(cntrlDarrMnSpc)

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
labelsToPlot.append('Global SST PDFs in GLENS')

print(colorsToPlot)
plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, savePath, saveName)

print('Completed!')
