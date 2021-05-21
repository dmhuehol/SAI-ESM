"""
Plot difference between RCP8.5 ("Control") and SAI ("Feedback") values for a
GLENS output variable on a globe. Equal Earth map projection used by default.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
"""

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
import fun_convert_unit as fcu

# Inputs
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

startInt = [2010,2019]
finalInt = [2090,2099]
levOfInt = 1000 #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210520_ozoneAndRfctrng/'
savePrfx = 'globe_1p_FdbckCntrl_'
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

# Average over years
toiStart = dot.average_over_years(glensCntrlLoi, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
toiEndCntrl = dot.average_over_years(glensCntrlLoi, finalInt[0], finalInt[1])
toiEndFdbck = dot.average_over_years(glensFdbckLoi, finalInt[0], finalInt[1])
diffToiFdbck =  toiEndCntrl - toiEndFdbck

# Unit conversion
diffToiFdbckPlot = fcu.molmol_to_ppb(diffToiFdbck)

# Plotting
firstDcd = str(startInt[0]) + '-' + str(startInt[1])
lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
sceneStr = 'RCP8.5 - SAI'
levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)
varStr = glensDarrCntrl.long_name

CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
plt.figure(figsize=(12, 2.73*2))
ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
cmap = cmocean.cm.delta
minVal = -diffToiFdbckPlot.quantile(0.99).data
# minVal = -7 #Override automatic colorbar minimum here
maxVal = diffToiFdbckPlot.quantile(0.99).data
# maxVal = 7 #Override automatic colorbar maximum here

plt_tls.drawOnGlobe(ax, diffToiFdbckPlot, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' ' + sceneStr + ' ' + levStr + ' ' + varStr)
# plt.title("2010-2019 Baseline - 2090-2099 SAI [50 0] ozone") #Override automatic title generation here
saveStr = savePrfx + dataKey + '_' + levStr + '_' + lastDcd
# saveStr = 'globe_1p_FdbckCntrl_O3_[50 0]_C2010-2019_F2090-2099'#Override automatic filename generation here
savename = savePath + saveStr + '.png'
plt.savefig(savename, dpi=dpi_val, bbox_inches='tight')
ic(savename)

print('Completed!')
