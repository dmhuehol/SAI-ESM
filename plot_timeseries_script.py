# Plot timeseries for a GLENS dataset given location or region (eventually)

import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

import difference_over_time as dot
import plotting_tools as plt_tls
import process_glens_fun as pgf

# Inputs
filenameCntrl = 'control.001.U.r90x45.shift.annual.nc'
filenameFdbck = 'feedback.001.U.r90x45.shift.annual.nc'
dataPath = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210204_timeseries/'
saveFile = 'timeseries_testU_'
saveName = savePath + saveFile
dpi_val = 300

levOfInt = 1000 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

# Control
glensDsetCntrl = xr.open_dataset(cntrlPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]

# Feedback
glensDsetFdbck = xr.open_dataset(fdbckPath)
glensDarrFdbck = glensDsetFdbck[dataKey]

# Obtain matching periods
# Format: time,lev,lat,lon
bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)

glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
cntrlToPlot = glensCntrlPoi.sel(lev=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]
fdbckToPlot = glensFdbckPoi.sel(lev=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
# glensFdbckToPlot = glensFdbckLoi.where((glensFdbckLoi.lat == -65.45) & (glensFdbckLoi.lon == 148))

# Plotting
yStr = glensDarrFdbck.units
varStr = glensDarrFdbck.long_name
startStr = str(bndDct['strtYrMtch'])
endStr = str(bndDct['endYrMtch'])
latStr = str(np.round_(cntrlToPlot.lat.data,decimals=2))
lonStr = str(np.round_(cntrlToPlot.lon.data,decimals=2))
locStr = latStr + '_' + lonStr
locTitleStr = '(' + latStr + ',' + lonStr + ')'

plt.figure(figsize=(12,2.73*2))
# ax = plt.subplot(2,2,1) #nrow ncol index
plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
plt.legend()
plt.ylabel(yStr)
plt.title(varStr + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)
# plt.show()

plt.savefig(saveName + locStr + '.png',dpi=dpi_val,bbox_inches='tight')

print('Completed!')
