
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

# Inputs
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

levOfInt = 'total'#5.96 * 10**(-6) #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
latOfInt = 34
lonOfInt = 282
regionToPlot = 'point' #aspirational

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210518_tsRfctrng/'
saveFile = 'timeseries_testO3_'
saveName = savePath + saveFile
dpi_val = 400

glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)

dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)
glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
ic(glensCntrlPoi['lev'])
glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]
levs = glensCntrlPoi['lev'].data
# levMask = levs > levActive
# ic(levOfInt)

if levOfInt == 'total':
    glensCntrlPoi = glensCntrlPoi.sum(dim='lev')
    glensFdbckPoi = glensFdbckPoi.sum(dim='lev')
elif levOfInt == 'troposphere':
    indTpause = pgf.find_closest_level(glensCntrlPoi, 200, levName='lev') #simple split on 200hPa for now
    levMask = levs > indTpause
    glensCntrlPoi = glensCntrlPoi[:,levMask,:,:]
    glensCntrlPoi = glensCntrlPoi.sum(dim='lev')
    glensFdbckPoi = glensFdbckPoi[:,levMask,:,:]
    glensFdbckPoi = glensFdbckPoi.sum(dim='lev')
elif levOfInt == 'stratosphere':
    indTpause = pgf.find_closest_level(glensCntrlPoi, 200, levName='lev') #simple split on 200hPa for now
    levMask = levs <= indTpause
    glensCntrlPoi = glensCntrlPoi[:,levMask,:,:]
    glensCntrlPoi = glensCntrlPoi.sum(dim='lev')
    glensFdbckPoi = glensFdbckPoi[:,levMask,:,:]
    glensFdbckPoi = glensFdbckPoi.sum(dim='lev')
elif np.size(levOfInt)==2:
    indHghrPres = pgf.find_closest_level(glensCntrlPoi, max(levOfInt), levName='lev')
    indLowerPres = pgf.find_closest_level(glensCntrlPoi, min(levOfInt), levName='lev')
    levOfInt = [levs[indLowerPres],levs[indHghrPres]]
    glensCntrlPoi = glensCntrlPoi[:,indLowerPres:indHghrPres+1,:,:]
    glensCntrlPoi = glensCntrlPoi.sum(dim='lev')
    glensFdbckPoi = glensFdbckPoi[:,indLowerPres:indHghrPres+1,:,:]
    glensFdbckPoi = glensFdbckPoi.sum(dim='lev')
else:
    indClosest = pgf.find_closest_level(glensCntrlPoi, levOfInt, levName='lev')
    levActive = glensCntrlPoi['lev'].data[indClosest]
    glensCntrlPoi = glensCntrlPoi.sel(lev=levActive)
    glensFdbckPoi = glensFdbckPoi.sel(lev=levActive)

if regionToPlot == 'point':
    cntrlToPlot = glensCntrlPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")
    fdbckToPlot = glensFdbckPoi.sel(lat=latOfInt, lon=lonOfInt, method="nearest")

# Plotting
yStr = glensDarrFdbck.units
varStr = glensDarrFdbck.long_name
startStr = str(bndDct['strtYrMtch'])
endStr = str(bndDct['endYrMtch'])
if isinstance(levOfInt,str):
    levStr = levOfInt
elif np.size(levOfInt)==2:
    if (np.round_(levOfInt[0],decimals=1)==0) | (np.round_(levOfInt[1],decimals=1)==0):
        levStr = str(np.round_(levOfInt,decimals=6))
    else:
        levStr = str(np.round_(levOfInt,decimals=1))
elif np.round_(levActive,decimals=1) == 0:
    levStr = str(np.round_(levActive,decimals=6))
else:
    levStr = str(np.round_(levActive,decimals=1))
ic(levStr)
latStr = str(np.round_(cntrlToPlot.lat.data,decimals=2))
lonStr = str(np.round_(cntrlToPlot.lon.data,decimals=2))
locStr = latStr + '_' + lonStr
locTitleStr = '(' + latStr + ',' + lonStr + ')'

plt.figure()
plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
plt.legend()
plt.ylabel(yStr)
plt.autoscale(enable=True, axis='x', tight=True)
plt.title(varStr + ' ' + levStr + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)
plt.savefig(saveName + locStr + '_' + levStr + '.png',dpi=dpi_val,bbox_inches='tight')

print("Completed! :D")
