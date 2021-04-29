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
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210429_refactoringAndOzone/'
saveFile = 'timeseries_testO3_'
saveName = savePath + saveFile
dpi_val = 400

# levOfInt = 1000 #hPa
latOfInt = 34
lonOfInt = 282
regionToPlot = 'point' #aspirational

# Control
glensDsetCntrl = xr.open_dataset(cntrlPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]

# Feedback
glensDsetFdbck = xr.open_dataset(fdbckPath)
glensDarrFdbck = glensDsetFdbck[dataKey]

bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)
glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
print(glensCntrlPoi['lev'])
glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]
levOfInt = glensCntrlPoi['lev'].data[len(glensCntrlPoi['lev'].data)-1]
print(levOfInt)

if regionToPlot == 'point':
    # glensCntrlPoi = glensCntrlPoi.sum(dim='lev')#
    cntrlToPlot = glensCntrlPoi.sel(lev=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
    # cntrlToPlot = glensCntrlPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")
    # glensFdbckPoi = glensFdbckPoi.sum(dim='lev')
    fdbckToPlot = glensFdbckPoi.sel(lev=levOfInt, lat=latOfInt, lon=lonOfInt, method="nearest")
    # fdbckToPlot = glensFdbckPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")

# Plotting
yStr = glensDarrFdbck.units
varStr = glensDarrFdbck.long_name
startStr = str(bndDct['strtYrMtch'])
endStr = str(bndDct['endYrMtch'])
levStr = str(np.round_(levOfInt,decimals=1))
latStr = str(np.round_(cntrlToPlot.lat.data,decimals=2))
lonStr = str(np.round_(cntrlToPlot.lon.data,decimals=2))
locStr = latStr + '_' + lonStr
locTitleStr = '(' + latStr + ',' + lonStr + ')'

plt.figure()
plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
plt.legend()
plt.ylabel(yStr)
# plt.title('TCO' + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)
plt.autoscale(enable=True, axis='x', tight=True)
plt.title(varStr + ' ' + levStr + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)
plt.savefig(saveName + locStr + '_' + levStr + '.png',dpi=dpi_val,bbox_inches='tight')

print("Completed! :D")
