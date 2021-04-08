# Plot GLENS data on a globe given location or region (eventually)

import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np
import sys

import difference_over_time as dot
import process_glens_fun as pgf
import plotting_tools as plt_tls

# Inputs
filenameCntrl = 'annual_control.nc'#'b.e15.B5505C5WCCML45BGCR.f09_g16.control.001.cam.h0.T.201001-201912.nc'
filenameFdbck = 'annual_feedback.nc'#'b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.001.cam.h0.T.209001-209912.nc'
dataPath = '/Users/dhueholt/Documents/GLENS_data/attemptMap/'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210401_dataformat/'
dpi_val = 400

startInt = [2010,2019]
finalInt = [2090,2099]
levOfInt = 1000 #hPa
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

# Control
glensDsetCntrl = xr.open_dataset(cntrlPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey] #HERE IT IS
print(glensDarrCntrl)
toiStart = dot.average_over_years(glensDarrCntrl,startInt[0],startInt[1]) # 2010-2019 is baseline, injection begins 2020
toiEndCntrl = dot.average_over_years(glensDarrCntrl,finalInt[0],finalInt[1])
# diffToiCntrl = toiEndCntrl - toiStart
#format: lev,lat,lon

# Feedback
glensDsetFdbck = xr.open_dataset(fdbckPath)
glensDarrFdbck = glensDsetFdbck[dataKey]
toiEndFdbck = dot.average_over_years(glensDarrFdbck,finalInt[0],finalInt[1])
# diffToiFdbck = toiEndFdbck - toiStart
diffToiFdbck = toiEndFdbck - toiEndCntrl

# Calculate
diffToiFdbck = diffToiFdbck.sel(lev=levOfInt)
# diffToiLevCntrl = diffToiCntrl.sel(lev=levOfInt)
# diffToiLevFdbck = diffToiFdbck.sel(lev=levOfInt)

# diffToiLevFdbckCntrl = diffToiLevFdbck - diffToiLevCntrl

# diffToiEndFdbckCntrl = toiEndFdbck - toiEndCntrl
# diffToiEndLevFdbckCntrl = diffToiEndFdbckCntrl.sel(lev=levOfInt)
#
# diffToiEndLevFdbckCntrlAbs = np.abs(diffToiEndLevFdbckCntrl)
# diffToiEndLevFdbckCntrlAbsNorm = diffToiEndLevFdbckCntrlAbs / np.max(diffToiEndLevFdbckCntrlAbs)
# quantCut = diffToiEndLevFdbckCntrlAbsNorm.quantile(quantileForFig)
# diffToiEndLevFdbckCntrlQ = diffToiEndLevFdbckCntrlAbsNorm
# diffToiEndLevFdbckCntrlQ.data[diffToiEndLevFdbckCntrlQ.data < quantCut.data] = np.nan
#
# diffToiEndFdbckCntrl = toiEndFdbck - toiEndCntrl
# diffToiEndLevFdbckCntrl = diffToiEndFdbckCntrl.sel(lev=levOfInt)

# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

# Make title text
firstDcd = str(startInt[0]) + '-' + str(startInt[1])
lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
cntrlStr = 'RCP8.5'
# fdbckStr = 'SAI'
tempStr = 'SAI - RCP8.5'
levStr = str(levOfInt) + 'mb'
varStr = glensDarrCntrl.long_name
quantileStr = str(quantileForFig)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(1,1,1,projection=mapProj) #nrow ncol index
cmap = cmocean.cm.haline
cmapSeq = cm.Purples
minVal = -8#-diffToiFdbck.quantile(0.97).data
maxVal = 1#diffToiFdbck.quantile(0.97).data

plt_tls.drawOnGlobe(ax, diffToiFdbck, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
# plt.title(lastDcd + ' - ' + firstDcd + ' ' + tempStr + ' ' + levStr + ' ' + varStr)
plt.title(lastDcd + ' ' + tempStr + ' ' + levStr + ' ' + varStr)

prfxStr = '1p_FdbckCntrl_'
saveStr = prfxStr + dataKey + str(levOfInt) + '_' + regionToPlot + '_' + str(startInt[0]) + str(startInt[1]) + '_' + str(finalInt[0]) + str(finalInt[1])
savename = savePath + saveStr + '.png'
plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')

print('Completed!')
