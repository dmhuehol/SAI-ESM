# This makes timespan-avg pdfs for data that's in single-file format as opposed to decadal-file format.
# This is to test the ability to produce ENSEMBLE MEANS
# This will be combined with plot_decadalpdf_script after completion.
# Also, this makes SST plots

# latOfInt/lonOfInt can be a np.array([min,max]), or a single value. If
# regionToPlot is set to 'global', latOfInt/lonOfInt are ignored.

import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sn
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np
import scipy.stats as stats
from glob import glob
import time

import difference_over_time as dot
import plotting_tools as plt_tls
import process_glens_fun as pgf
import region_library as rlib

# Inputs
cntrlFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/control.001.SST.r90x45.shift.annual.nc'
fdbckFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/feedback.001.SST.r90x45.shift.annual.nc'

cntrlFile2 = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/control.002.SST.r90x45.shift.annual.nc'
fdbckFile2 = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/feedback.002.SST.r90x45.shift.annual.nc'

baselineFlag = 0
regionToPlot = 'regional' #'global' 'regional' 'point'
regOfInt = rlib.Nino34()
levOfInt = 500 #z_t for SST, often hPa for other data
latOfInt = regOfInt['regLats']#np.array([-35,-22])
lonOfInt = regOfInt['regLons']#np.array([108,115])
cntrlIntToPlot = [2020,2045,2070]#[2020,2030,2040,2050,2090]#[2020,2050]
fdbckIntToPlot = [2020,2045,2070]#[2020,2030,2040,2050,2090]#[2020,2050]
timePeriod = 30 #number of years, i.e. 10 = decade

plotStyle = 'kde' #'kde' or 'hist'
areaAvgBool = False
titleStr = regOfInt['regStr'] + ' SST PDFs ensemble mean members 1&2'
# titleStr = 'Gulf of Mexico SST PDFs in GLENS' #use when region is set manually
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210318_regionrefinement/ens/'
saveName = 'pdf_' + plotStyle + '_SST_cntrlfdbck_' + regOfInt['regSaveStr'] + '_30yr_ENSMEAN12'
# saveName = 'pdf_hist_SST_cntrlfdbck_REGIONHERE_30yr_MEANTEST' #use when region is set manually
dpiVal = 400

cntrlDset = xr.open_dataset(cntrlFile)
cntrlDset2 = xr.open_dataset(cntrlFile2)
fdbckDset = xr.open_dataset(fdbckFile)
fdbckDset2 = xr.open_dataset(fdbckFile2)

if regionToPlot == 'global':
    cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is equivalent to level
    dataKey = pgf.discover_data_var(cntrlDset)
    cntrlDarr = cntrlDsetLoi[dataKey]

    fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt)
    fdbckDarr = fdbckDsetLoi[dataKey]

    if areaAvgBool:
        print("Spatially averaging across globe")
        cntrlDarrMnSpc = cntrlDarr.mean(dim=['lat','lon'])
        fdbckDarrMnSpc = fdbckDarr.mean(dim=['lat','lon'])
    else:
        print("No spatial average applied")
        cntrlDarrMnSpc = cntrlDarr.stack(spc=("lat","lon"))
        fdbckDarrMnSpc = fdbckDarr.stack(spc=("lat","lon"))

elif regionToPlot == 'regional':
    cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is equivalent to level
    dataKey = pgf.discover_data_var(cntrlDset)
    cntrlDarrLoi = cntrlDsetLoi[dataKey]
    lats = cntrlDsetLoi['lat'] #feedback and control are on same grid, fortunately
    lons = cntrlDsetLoi['lon']
    latMask = (lats>latOfInt[0]) & (lats<latOfInt[1])
    lonMask = (lons>lonOfInt[0]) & (lons<lonOfInt[1])
    cntrlDarrLoiAoi = cntrlDarrLoi[:,latMask,lonMask]

    fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt)
    dataKey = pgf.discover_data_var(fdbckDset)
    fdbckDarrLoi = fdbckDsetLoi[dataKey]
    fdbckDarrLoiAoi = fdbckDarrLoi[:,latMask,lonMask]

    cntrlDsetLoi2 = cntrlDset2.sel(z_t=levOfInt) #z_t is equivalent to level
    cntrlDarrLoi2 = cntrlDsetLoi2[dataKey]
    cntrlDarrLoiAoi2 = cntrlDarrLoi2[:,latMask,lonMask]

    fdbckDsetLoi2 = fdbckDset2.sel(z_t=levOfInt)
    fdbckDarrLoi2 = fdbckDsetLoi2[dataKey]
    fdbckDarrLoiAoi2 = fdbckDarrLoi2[:,latMask,lonMask]

    if areaAvgBool:
        print("Spatially averaging across region")
        cntrlDarrMnSpc = cntrlDarrLoiAoi.mean(dim=['lat','lon'])
        fdbckDarrMnSpc = fdbckDarrLoiAoi.mean(dim=['lat','lon'])
    else:
        print("No spatial average applied")
        cntrlDarrMnSpc = cntrlDarrLoiAoi.stack(spc=("lat","lon"))
        fdbckDarrMnSpc = fdbckDarrLoiAoi.stack(spc=("lat","lon"))

        cntrlDarrMnSpc2 = cntrlDarrLoiAoi2.stack(spc=("lat","lon"))
        fdbckDarrMnSpc2 = fdbckDarrLoiAoi2.stack(spc=("lat","lon"))

        # ens mean
        cntrlDarrMnSpc = xr.concat([cntrlDarrMnSpc, cntrlDarrMnSpc2], dim="mem")
        fdbckDarrMnSpc = xr.concat([fdbckDarrMnSpc, fdbckDarrMnSpc2], dim="mem")
        cntrlDarrMnSpc = cntrlDarrMnSpc.mean(dim="mem")
        fdbckDarrMnSpc = fdbckDarrMnSpc.mean(dim="mem")
        print(cntrlDarrMnSpc)
        # cntrlDarrMnSpc = np.mean(cntrlDarrMnSpc,cntrlDarrMnSpc2) #ens mean
        # fdbckDarrMnSpc = np.mean(fdbckDarrMnSpc,fdbckDarrMnSpc2)

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

# Remove 2010-2019 average
# baselineMeanToRmv = dot.average_over_years(cntrlDarrMnSpc,2010,2019)
# cntrlDarrMnSpcNorm = cntrlDarrMnSpc - baselineMeanToRmv
# fdbckDarrMnSpcNorm = fdbckDarrMnSpc - baselineMeanToRmv

iqr = stats.iqr(cntrlDarrMnSpc)
# binwidth = 2*iqr*(10 ** -1/3) # the Freedman-Diaconis rule
binwidth = 0.2 #the Let's Not Overthink This rule
print(binwidth)

cntrlActive = cntrlDarrMnSpc
fdbckActive = fdbckDarrMnSpc

# Extract the decades of interest from the control and feedback datasets
cntrlYears = cntrlActive['time'].dt.year.data
cntrlHandlesToPlot = list()
cntrlHandlesToPlot = pgf.extract_doi(cntrlIntToPlot, cntrlYears, timePeriod, cntrlActive, cntrlHandlesToPlot)
fdbckYears = fdbckActive['time'].dt.year.data
fdbckHandlesToPlot = list()
fdbckHandlesToPlot = pgf.extract_doi(fdbckIntToPlot, fdbckYears, timePeriod, fdbckActive, fdbckHandlesToPlot)
handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot

# If not applying a spatial average, flatten data so dimensions don't confuse plotting code
if ~areaAvgBool:
    for ind, h in enumerate(handlesToPlot):
        handlesToPlot[ind] = h.data.flatten()

# Generate colors and strings for plots and filenames
if baselineFlag:
    colorsToPlot = plt_tls.select_colors(baselineFlag,len(cntrlIntToPlot)-1,len(fdbckIntToPlot))
else:
    colorsToPlot = plt_tls.select_colors(baselineFlag,len(cntrlIntToPlot),len(fdbckIntToPlot))
if baselineFlag:
    labelsToPlot = list(['2010-2019 Baseline'])
else:
    labelsToPlot = list()
labelsToPlot = plt_tls.generate_labels(labelsToPlot, cntrlIntToPlot, timePeriod, 'RCP8.5')
labelsToPlot = plt_tls.generate_labels(labelsToPlot, fdbckIntToPlot, timePeriod, 'SAI')
labelsToPlot.append(titleStr)

print(colorsToPlot) # For troubleshooting

# Make KDE or histograms
if plotStyle == 'kde':
    plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, savePath, saveName, dpiVal)
elif plotStyle == 'hist':
    plt_tls.plot_pdf_hist(handlesToPlot, colorsToPlot, labelsToPlot, savePath, saveName, binwidth, dpiVal)
elif plotStyle == 'step':
    plt_tls.plot_pdf_step(handlesToPlot, colorsToPlot, labelsToPlot, savePath, saveName, binwidth, dpiVal)
else:
    print("Invalid plot style") #TODO: make this a formal error message

print('Completed!')
