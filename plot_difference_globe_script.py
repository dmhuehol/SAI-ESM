''' plot_difference_globe_script
Plot differences between RCP8.5 ("Control") and SAI ("Feedback") values for a
GLENS output variable on a 4-panel globe.

Equal Earth map projection used by default.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import xarray as xr
xr.set_options(keep_attrs=True)
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
levOfInt = 200 #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
quantileOfInt = 0.67

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210525_github/'
savePrfx = 'globe_4p_FdbckCntrl_'
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

# Unit conversion
glensCntrlLoi = fcu.molmol_to_ppm(glensCntrlLoi)
glensFdbckLoi = fcu.molmol_to_ppm(glensFdbckLoi)

# Average over years
toiStart = dot.average_over_years(glensCntrlLoi, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
toiEndCntrl = dot.average_over_years(glensCntrlLoi, finalInt[0], finalInt[1])
toiEndFdbck = dot.average_over_years(glensFdbckLoi, finalInt[0], finalInt[1])

# Calculate 4-panel values
diffToiCntrl = toiEndCntrl - toiStart
diffToiFdbck = toiEndFdbck - toiStart
diffEndCntrlFdbck = toiEndCntrl - toiEndFdbck
diffEndCntrlFdbckAbsNormQ = pgf.isolate_change_quantile(diffEndCntrlFdbck, quantileOfInt)

# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

# Make title text
firstDcd = str(startInt[0]) + '-' + str(startInt[1])
lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
cntrlStr = 'RCP8.5'
fdbckStr = 'SAI'
levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)
varStr = glensDarrCntrl.long_name
quantileStr = str(quantileOfInt)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
cmap = cmocean.cm.delta
cmapSeq = cmocean.cm.dense
minVal = -diffToiCntrl.quantile(0.99).data
maxVal = diffToiCntrl.quantile(0.99).data

plt_tls.drawOnGlobe(ax, diffToiCntrl, glensDarrCntrl.lat, glensDarrCntrl.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + levStr + ' ' + varStr)

ax2 = plt.subplot(2,2,2,projection=mapProj)
plt_tls.drawOnGlobe(ax2, diffToiFdbck, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + levStr + ' ' + varStr)

ax3 = plt.subplot(2,2,3,projection=mapProj)
plt_tls.drawOnGlobe(ax3, diffEndCntrlFdbck, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + levStr + ' ' + varStr)

ax4 = plt.subplot(2,2,4,projection=mapProj)
plt_tls.drawOnGlobe(ax4, diffEndCntrlFdbckAbsNormQ, glensDarrFdbck.lat, glensDarrFdbck.lon, cmapSeq, vmin=0, vmax=1, cbarBool=True, fastBool=True, extent='max')
plt.title(lastDcd + ' ' + fdbckStr + ' - ' + cntrlStr + ' ' + levStr + ' ' + '|' + 'norm' + '\u0394' + varStr + '|' + '>' + quantileStr + 'Q')

saveStr = savePrfx + dataKey + '_' + str(levOfInt) + '_' + str(startInt[0]) + str(startInt[1]) + '_' + str(finalInt[0]) + str(finalInt[1])
savename = savePath + saveStr + '.png'
plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')
ic(savename)

print('Completed!')
