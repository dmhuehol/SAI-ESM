# Plot GLENS data on a globe given location or region (eventually)

from icecream import ic
import sys

import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

import difference_over_time as dot
import process_glens_fun as pgf
import plotting_tools as plt_tls

# Inputs
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

startInt = [2020,2029]
finalInt = [2090,2099]
levOfInt = 1000 #hPa
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210519_regionsAndOzone/'
savePrfx = 'INPROGRESS_globe_1p_FdbckCntrl_'
dpi_val = 400

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

# Obtain levels
glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, levOfInt)
glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, levOfInt)

# Average over years but also apply level mask (TODO: REFACTOR)
toiStart = dot.average_over_years(glensCntrlLoi,startInt[0],startInt[1]) # 2010-2019 is baseline, injection begins 2020
toiEndCntrl = dot.average_over_years(glensCntrlLoi,finalInt[0],finalInt[1])
toiEndFdbck = dot.average_over_years(glensFdbckLoi,finalInt[0],finalInt[1])
diffToiFdbck =  toiEndFdbck - toiEndCntrl

# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
# Make title text
firstDcd = str(startInt[0]) + '-' + str(startInt[1])
lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
cntrlStr = 'RCP8.5'
# fdbckStr = 'SAI'
tempStr = 'RCP8.5 - SAI'
levStr = str(levOfInt) + 'mb'
varStr = glensDarrCntrl.long_name
quantileStr = str(quantileForFig)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(1,1,1,projection=mapProj) #nrow ncol index
cmap = cmocean.cm.delta#cmocean.cm.haline
cmapSeq = cm.Purples
minVal = -diffToiFdbck.quantile(0.99).data
maxVal = diffToiFdbck.quantile(0.99).data

ic(diffToiFdbck)
plt_tls.drawOnGlobe(ax, diffToiFdbck, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' ' + tempStr + ' ' + levStr + ' ' + varStr)

saveStr = savePrfx + dataKey + str(levOfInt) + '_' + regionToPlot + '_' + str(startInt[0]) + str(startInt[1]) + '_' + str(finalInt[0]) + str(finalInt[1])
savename = savePath + saveStr + '.png'
plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')

print('Completed!')
