''' plot_pdf_script
Plot pdfs for RCP8.5 ("Control") and SAI ("Feedback") values for a GLENS output
variable. Three formats are available: a kernel density estimate, a histogram,
or a step plot.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

import sys

import xarray as xr
xr.set_options(keep_attrs=True)
import matplotlib.pyplot as plt
import seaborn as sn
import cmocean
import numpy as np
import scipy.stats as stats

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

baselineFlag = 0 #Include 2010-2019 period from RCP8.5 "Control" scenario
levOfInt = 'total' #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
regionToPlot = 'global' #'global', rlib.Place(), [latN,lonE360]
areaAvgBool = False #Apply spatial average over region or include every point individually
cntrlPoi = [2020,2090]#[2020,2030,2040,2050,2090]#[2020,2050]
fdbckPoi = [2020,2090]#[2020,2030,2040,2050,2090]#[2020,2050]
timePeriod = 10 #number of years, i.e. 10 = decade

plotStyle = 'kde' #'kde' or 'hist' or 'step'
savePath = '/Users/dhueholt/Documents/GLENS_fig/20210525_github/'
savePrfx = 'pdf_' + plotStyle
dpiVal = 400

# Open data
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
# baselineMeanToRmv = dot.average_over_years(glensCntrlAoi,2010,2019)
# glensCntrlAoi = glensCntrlAoi - baselineMeanToRmv
# glensFdbckAoi = glensFdbckAoi - baselineMeanToRmv

# Unit conversion
cntrlToPlot = fcu.molmol_to_ppm(glensCntrlAoi)
fdbckToPlot = fcu.molmol_to_ppm(glensFdbckAoi)

iqr = stats.iqr(glensCntrlAoi)
# binwidth = 2*iqr*(10 ** -1/3) # the Freedman-Diaconis rule
binwidth = 10 #the Let's Not Overthink This rule

# Extract the decades of interest from the control and feedback datasets
cntrlYears = cntrlToPlot['time'].dt.year.data
cntrlHandlesToPlot = list()
cntrlHandlesToPlot = pgf.extract_doi(cntrlPoi, cntrlYears, timePeriod, cntrlToPlot, cntrlHandlesToPlot)
fdbckYears = fdbckToPlot['time'].dt.year.data
fdbckHandlesToPlot = list()
fdbckHandlesToPlot = pgf.extract_doi(fdbckPoi, fdbckYears, timePeriod, fdbckToPlot, fdbckHandlesToPlot)
handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot

# If not applying a spatial average, flatten data so dimensions don't confuse plotting code
if ~areaAvgBool:
    for ind, h in enumerate(handlesToPlot):
        handlesToPlot[ind] = h.data.flatten()

# Generate colors and strings for plots and filenames
if baselineFlag:
    colorsToPlot = plt_tls.select_colors(baselineFlag,len(cntrlPoi)-1,len(fdbckPoi))
else:
    colorsToPlot = plt_tls.select_colors(baselineFlag,len(cntrlPoi),len(fdbckPoi))
if baselineFlag:
    labelsToPlot = list(['2010-2019 Baseline'])
else:
    labelsToPlot = list()
labelsToPlot = plt_tls.generate_labels(labelsToPlot, cntrlPoi, timePeriod, 'RCP8.5')
labelsToPlot = plt_tls.generate_labels(labelsToPlot, fdbckPoi, timePeriod, 'SAI')
varStr = glensDarrFdbck.long_name
varSave = varStr.replace(" ","")
levStr = pgf.make_level_string(cntrlToPlot, levOfInt)
timeStr = str(timePeriod) + 'yr'
titleStr = varStr + ' ' + levStr + ' ' + locTitleStr
labelsToPlot.append(titleStr)
if areaAvgBool:
    spcStr = 'spcavg'
else:
    spcStr = 'nospcavg'
unit = cntrlToPlot.attrs['units']
saveName = savePath + savePrfx + '_' + timeStr + '_' + varSave + '_' + levStr + '_' + locStr + '_' + spcStr

# Make kde, histograms, or step plots
if plotStyle == 'kde':
    plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, unit, saveName, dpiVal)
elif plotStyle == 'hist':
    plt_tls.plot_pdf_hist(handlesToPlot, colorsToPlot, labelsToPlot, unit, saveName, binwidth, dpiVal)
elif plotStyle == 'step':
    plt_tls.plot_pdf_step(handlesToPlot, colorsToPlot, labelsToPlot, unit, saveName, binwidth, dpiVal)
else:
    sys.exit('Invalid plot style')
