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

levOfInt = 1000 #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
regionToPlot = [-70,282] #'global', rlib.Place(), [latN,lonE360]

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210521_rfctrngArea4pPdf/'
saveFile = 'NEWSTRING_timeseries_O3_'
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
cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlPoi, regionToPlot)
fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckPoi, regionToPlot)

# Unit conversion
cntrlToPlot = fcu.molmol_to_ppb(cntrlToPlot)
fdbckToPlot = fcu.molmol_to_ppb(fdbckToPlot)

# Plotting
yStr = cntrlToPlot.units
varStr = glensDarrFdbck.long_name
startStr = str(bndDct['strtYrMtch'])
endStr = str(bndDct['endYrMtch'])
levStr = pgf.make_level_string(glensCntrlPoi, levOfInt)
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
ic(saveName + locStr + '_' + levStr + '.png')

print("Completed! :D")
