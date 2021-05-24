# This makes timespan-avg pdfs for data that's in single-file format as opposed to decadal-file format.
# This will be combined with plot_decadalpdf_script after completion.
# Also, this makes SST plots

# latOfInt/lonOfInt can be a np.array([min,max]), or a single value. If
# regionToPlot is set to 'global', latOfInt/lonOfInt are ignored.
from icecream import ic
import sys

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
import fun_convert_unit as fcu

# Inputs
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

baselineFlag = 0
levOfInt = 'stratosphere' #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
regionToPlot = rlib.NoLandLatitude() #'global', rlib.Place(), [latN,lonE360]
areaAvgBool = False
cntrlIntToPlot = [2020,2090]#[2020,2030,2040,2050,2090]#[2020,2050]
fdbckIntToPlot = [2020,2090]#[2020,2030,2040,2050,2090]#[2020,2050]
timePeriod = 10 #number of years, i.e. 10 = decade

# regionToPlot = 'regional' #'global' 'regional' 'point'
# regOfInt = rlib.EasternEurope()
# levOfInt = 1000 #z_t for SST, often hPa for other data
# latOfInt = regOfInt['regLats']#np.array([-35,-22])
# lonOfInt = regOfInt['regLons']#np.array([108,115])

plotStyle = 'step' #'kde' or 'hist' or 'step'
# titleStr = regOfInt['regStr'] + ' T PDFs in GLENS'
titleStr = 'NoLandLat60S stratosphere ozone PDFs in GLENS' #use when region is set manually
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210524_4pPdfReg/'
# saveName = 'pdf_' + plotStyle + '_T_cntrlfdbck_' + regOfInt['regSaveStr'] + '_10yr_TESTINSET'
saveName = 'pdf_step_stratosphere_O3_cntrlfdbck_NoLandLat60S_10yr_nospcavg' #use when region is set manually
dpiVal = 400

glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
#dataKey = '' #Override automatic variable discovery here
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

# Obtain levels
glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, levOfInt)
glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, levOfInt)

# Deal with area
glensCntrlAoi, locStr, locTitleStr = pgf.manage_area(glensCntrlLoi, regionToPlot, areaAvgBool)
glensFdbckAoi, locStr, locTitleStr = pgf.manage_area(glensFdbckLoi, regionToPlot, areaAvgBool)

# Remove 2010-2019 average
# baselineMeanToRmv = dot.average_over_years(cntrlDarrMnSpc,2010,2019)
# cntrlDarrMnSpcNorm = cntrlDarrMnSpc - baselineMeanToRmv
# fdbckDarrMnSpcNorm = fdbckDarrMnSpc - baselineMeanToRmv

# Unit conversion
cntrlToPlot = fcu.molmol_to_ppm(glensCntrlAoi)
fdbckToPlot = fcu.molmol_to_ppm(glensFdbckAoi)

iqr = stats.iqr(glensCntrlAoi)
# binwidth = 2*iqr*(10 ** -1/3) # the Freedman-Diaconis rule
binwidth = 1 #the Let's Not Overthink This rule
ic(binwidth)

cntrlActive = cntrlToPlot
fdbckActive = fdbckToPlot

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
