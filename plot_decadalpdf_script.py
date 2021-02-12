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
cntrlPath = '/Users/dhueholt/Documents/GLENS_data/control/'
fdbckPath = '/Users/dhueholt/Documents/GLENS_data/feedback/'
cntrlIn = glob(cntrlPath + '*.nc')
fdbckIn = glob(fdbckPath + '*.nc')

baselinePath = '/Users/dhueholt/Documents/GLENS_data/'
baselineFile = glob(baselinePath + '*.nc')
print(baselineFile)

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210211_formeeting/'
saveName = 'pdf_Tglobal_cntrlfdbck'
dpi_val = 300

levOfInt = 1000 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

baselineDs = xr.open_dataset(baselineFile[0])
baselineDsLoi = baselineDs.sel(lev=levOfInt)
dataKey = pgf.discover_data_var(baselineDsLoi)
baselineDarr = baselineDsLoi[dataKey]
baselineDarrMnSpc = baselineDarr.mean(dim=['lat','lon'])
baselineMeanToRmv = baselineDarrMnSpc.mean(dim='time')
baselineNormDarr = baselineDarrMnSpc - baselineMeanToRmv

handlesToPlot = [baselineNormDarr.data]

for c,cfile in enumerate(cntrlIn):
    actvCntrl = xr.open_dataset(cfile)
    actvCntrlLoi = actvCntrl.sel(lev=levOfInt)
    actvCntrlDarr = actvCntrlLoi[dataKey]
    actvCntrlDarrMnSpc = actvCntrlDarr.mean(dim=['lat','lon'])
    actvNormCntrlDarr = actvCntrlDarrMnSpc - baselineMeanToRmv
    handlesToPlot.append(actvNormCntrlDarr.data)
    #for each file, repeat the above process

# Handle feedback + control separately; they may not have the same number of files
for f,ffile in enumerate(fdbckIn):
    actvFdbck = xr.open_dataset(ffile)
    actvFdbckLoi = actvFdbck.sel(lev=levOfInt)
    actvFdbckDarr = actvFdbckLoi[dataKey]
    actvFdbckDarrMnSpc = actvFdbckDarr.mean(dim=['lat','lon'])
    actvNormFdbckDarr = actvFdbckDarrMnSpc - baselineMeanToRmv
    handlesToPlot.append(actvNormFdbckDarr.data)


# print(handlesToPlot)

plt_tls.plot_pdf_kdeplot(handlesToPlot, ['slateblue','rosybrown','lightcoral','firebrick','orchid','purple'], ['2010-2019 Baseline','2020-2029 RCP8.5', '2080-2089 RCP8.5', '2090-2099 RCP8.5', '2080-2089 SAI', '2090-2099 SAI', 'Global temperature PDFs in GLENS'], savePath, saveName)


print('Completed!')
