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

import difference_over_time as dot
import plotting_tools as plt_tls
import process_glens_fun as pgf

# Inputs
cntrlFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/control.001.SST.r90x45.shift.annual.nc'
fdbckFile = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/feedback.001.SST.r90x45.shift.annual.nc'

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210211_formeeting/'
saveName = 'pdf_SSTglobal_cntrlfdbck'
dpi_val = 300

levOfInt = 500 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

# request periods as input e.g. [2030,2039,2090,2099]

cntrlDset = xr.open_dataset(cntrlFile)
cntrlDsetLoi = cntrlDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
dataKey = pgf.discover_data_var(cntrlDset)
cntrlDarr = cntrlDsetLoi[dataKey]
cntrlDarrMnSpc = cntrlDarr.mean(dim=['lat','lon'])

fdbckDset = xr.open_dataset(fdbckFile)
fdbckDsetLoi = fdbckDset.sel(z_t=levOfInt) #z_t is the equivalent to level, I guess?
fdbckDarr = fdbckDsetLoi[dataKey]
fdbckDarrMnSpc = fdbckDarr.mean(dim=['lat','lon'])
# print(cntrlDarrMnSpc)

# take 2010-2019 average and remove this
baselineMeanToRmv = dot.average_over_years(cntrlDarrMnSpc,2010,2019)
cntrlDarrMnSpcNorm = cntrlDarrMnSpc - baselineMeanToRmv
# print(cntrlDarrMnSpcNorm['time'].dt.year.data)
baselineToPlot = cntrlDarrMnSpcNorm[0:10]
cntrl20502059ToPlot = cntrlDarrMnSpcNorm[40:50]
cntrl20902099ToPlot = cntrlDarrMnSpcNorm[80:90]

fdbckDarrMnSpcNorm = fdbckDarrMnSpc - baselineMeanToRmv
fdbck20502059ToPlot = fdbckDarrMnSpcNorm[30:40]
fdbck20902099ToPlot = fdbckDarrMnSpcNorm[70:80]
#
# handlesToPlot = [baselineNormDarr.data]
#
# for c,cfile in enumerate(cntrlIn):
#     actvCntrl = xr.open_dataset(cfile)
#     actvCntrlLoi = actvCntrl.sel(lev=levOfInt)
#     actvCntrlDarr = actvCntrlLoi[dataKey]
#     actvCntrlDarrMnSpc = actvCntrlDarr.mean(dim=['lat','lon'])
#     actvNormCntrlDarr = actvCntrlDarrMnSpc - baselineMeanToRmv
#     handlesToPlot.append(actvNormCntrlDarr.data)
#     #for each file, repeat the above process
#
# # Handle feedback + control separately; they may not have the same number of files
# for f,ffile in enumerate(fdbckIn):
#     actvFdbck = xr.open_dataset(ffile)
#     actvFdbckLoi = actvFdbck.sel(lev=levOfInt)
#     actvFdbckDarr = actvFdbckLoi[dataKey]
#     actvFdbckDarrMnSpc = actvFdbckDarr.mean(dim=['lat','lon'])
#     actvNormFdbckDarr = actvFdbckDarrMnSpc - baselineMeanToRmv
#     handlesToPlot.append(actvNormFdbckDarr.data)
#
#
# # print(handlesToPlot)
#
# plt_tls.plot_pdf_kdeplot(handlesToPlot, ['slateblue','rosybrown','lightcoral','firebrick','orchid','purple'], ['2010-2019 Baseline','2020-2029 RCP8.5', '2080-2089 RCP8.5', '2090-2099 RCP8.5', '2080-2089 SAI', '2090-2099 SAI', 'Global temperature PDFs in GLENS'], savePath, saveName)

plt_tls.plot_pdf_kdeplot((baselineToPlot, cntrl20502059ToPlot, cntrl20902099ToPlot, fdbck20502059ToPlot, fdbck20902099ToPlot), ['slateblue','lightcoral','firebrick','orchid','purple'], ['2010-2019 Baseline', '2050-2059 RCP8.5', '2090-2099 RCP8.5', '2050-2059 SAI', '2090-2099 SAI', 'Global SST PDFs in GLENS'], savePath, saveName)

print('Completed!')
