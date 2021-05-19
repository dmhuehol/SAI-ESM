# Make timeseries from GLENS output showing progression of both RCP8.5 ("Control)
# and SAI ("Feedback") for a variable.
#
# Written by Daniel Hueholt | May 2021
# Graduate Research Assistant at Colorado State University

from icecream import ic
import sys

import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

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

levOfInt = 'stratosphere' #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
regionToPlot = 'regional' #'global', 'regional', 'point'
regOfInt = rlib.EasternEurope()
latOfInt = regOfInt['regLats']#34
lonOfInt = regOfInt['regLons']#282

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210519_regionsAndOzone/'
saveFile = 'JUSTNOW_timeseries_testO3_'
saveName = savePath + saveFile
dpi_val = 400

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)
glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
ic(glensCntrlPoi['lev'])
glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]

# Obtain levels
glensCntrlPoi = pgf.obtain_levels(glensCntrlPoi, levOfInt)
glensFdbckPoi = pgf.obtain_levels(glensFdbckPoi, levOfInt)

# Deal with area (potentially break off to new function)
if regionToPlot == 'global':
    ic('global')
    latWeights = np.cos(np.deg2rad(glensCntrlPoi['lat']))
    glensCntrlPoiWght = glensCntrlPoi.weighted(latWeights)
    glensFdbckPoiWght = glensFdbckPoi.weighted(latWeights)
    cntrlToPlot = glensCntrlPoiWght.mean(dim=['lat','lon'])
    fdbckToPlot = glensFdbckPoiWght.mean(dim=['lat','lon'])
elif regionToPlot == 'regional':
    ic('regional')
    lats = glensCntrlPoi['lat'] #feedback and control are on same grid, fortunately
    lons = glensCntrlPoi['lon']
    latMask = (lats>latOfInt[0]) & (lats<latOfInt[1])
    lonMask = (lons>lonOfInt[0]) & (lons<lonOfInt[1])
    cntrlBoxMask = glensCntrlPoi[:,latMask,lonMask]
    fdbckBoxMask = glensFdbckPoi[:,latMask,lonMask]
    cntrlToPlot = cntrlBoxMask.mean(dim=['lat','lon'])
    fdbckToPlot = fdbckBoxMask.mean(dim=['lat','lon'])
elif regionToPlot == 'point':
    ic('point')
    cntrlToPlot = glensCntrlPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")
    fdbckToPlot = glensFdbckPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")
else:
    sys.exit('Input error! Check value for regionToPlot.')

# Unit conversion
cntrlToPlot = fcu.molmol_to_ppb(cntrlToPlot)
fdbckToPlot = fcu.molmol_to_ppb(fdbckToPlot)

# Plotting
yStr = cntrlToPlot.units
varStr = glensDarrFdbck.long_name
startStr = str(bndDct['strtYrMtch'])
endStr = str(bndDct['endYrMtch'])
levStr = pgf.make_level_string(glensCntrlPoi, levOfInt)

if regionToPlot == 'global':
    locStr = 'global'
    locTitleStr = 'global'
elif regionToPlot == 'regional':
    locStr = regOfInt['regSaveStr']
    locTitleStr = regOfInt['regStr']
elif regionToPlot == 'point':
    latStr = str(np.round_(cntrlToPlot.lat.data,decimals=2))
    lonStr = str(np.round_(cntrlToPlot.lon.data,decimals=2))
    locStr = latStr + '_' + lonStr
    locTitleStr = '(' + latStr + ',' + lonStr + ')'
ic(levStr, locStr)

# Make timeseries
plt.figure()
plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
plt.legend()
plt.ylabel(yStr)
plt.autoscale(enable=True, axis='x', tight=True)
plt.title(varStr + ' ' + levStr + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)
plt.savefig(saveName + locStr + '_' + levStr + '.png',dpi=dpi_val,bbox_inches='tight')

print("Completed! :D")
