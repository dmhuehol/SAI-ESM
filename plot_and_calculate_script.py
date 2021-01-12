import xarray as xr
import difference_over_time as dot
import plotting_tools as plt_tls
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

filenameCntrl = 'control.001.T.r90x45.shift.annual.nc'
glensDatasetCntrl = dot.import_glens_dataset(filenameCntrl)
startInt = [2020,2030]
finalInt = [2090,2099]
toiStartCntrl = dot.average_over_years(glensDatasetCntrl,startInt[0],startInt[1])
toiEndCntrl = dot.average_over_years(glensDatasetCntrl,finalInt[0],finalInt[1])
tempDiffOverToiCntrl = toiEndCntrl - toiStartCntrl
tempDiffOverToi1000Cntrl = tempDiffOverToiCntrl[0,:,:]
# glensEndToi is formatted: lev,lat,lon
# tempDiff1000Cntrl = tempDiffCntrl[0,:,:]
# glensDatasetCntrl, tempDiffCntrl = dot.simple_diff_calc(glensDatasetCntrl)

filenameFdbck = 'feedback.001.T.r90x45.shift.annual.nc'
glensDatasetFdbck = dot.import_glens_dataset(filenameFdbck)
toiStartFdbck = dot.average_over_years(glensDatasetFdbck,startInt[0],startInt[1])
toiEndFdbck = dot.average_over_years(glensDatasetFdbck,finalInt[0],finalInt[1])
tempDiffOverToiFdbck = toiEndFdbck - toiStartFdbck
tempDiffOverToi1000Fdbck = tempDiffOverToiFdbck[0,:,:]
# glensDatasetFdbck, tempDiffFdbck = dot.simple_diff_calc(glensDatasetFdbck)
# tempDiff1000Fdbck = tempDiffFdbck[0,:,:]

tempDiffOverToi1000FdbckCntrl = tempDiffOverToi1000Fdbck - tempDiffOverToi1000Cntrl

tempDiffOverToiEndFdbckCntrl = toiEndFdbck - toiEndCntrl
tempDiffOverToiEnd1000FdbckCntrl = tempDiffOverToiEndFdbckCntrl[0,:,:]
# print(tempDiffOverToiEnd1000FdbckCntrl.equals(tempDiffOverToi1000FdbckCntrl))

quantCut = tempDiffOverToiEnd1000FdbckCntrl.quantile(0.33)
tempDiffOverToiEnd1000FdbckCntrlQ = tempDiffOverToiEnd1000FdbckCntrl
tempDiffOverToiEnd1000FdbckCntrlQ.data[tempDiffOverToiEnd1000FdbckCntrlQ.data > quantCut.data] = np.nan
# print(tempDiffOverToiEnd1000FdbckCntrlQ.data)

tempDiffOverToiEndFdbckCntrl = toiEndFdbck - toiEndCntrl
tempDiffOverToiEnd1000FdbckCntrl = tempDiffOverToiEndFdbckCntrl[0,:,:]
# print(tempDiffOverToiEnd1000FdbckCntrl.equals(tempDiffOverToi1000FdbckCntrl))

# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
cmap = cmocean.cm.curl
minVal = -8
maxVal = 8
# print(maxVal)
plt_tls.drawOnGlobe(ax, tempDiffOverToi1000Cntrl, glensDatasetCntrl.lat, glensDatasetCntrl.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2090-2099 - 2020-2030 CONTROL (RCP8.5) 1000mb temp')

ax2 = plt.subplot(2,2,2,projection=mapProj)
plt_tls.drawOnGlobe(ax2, tempDiffOverToi1000Fdbck, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2090-2099 - 2020-2030 FEEDBACK (SAI) 1000mb temp')

ax3 = plt.subplot(2,2,3,projection=mapProj)
plt_tls.drawOnGlobe(ax3, tempDiffOverToiEnd1000FdbckCntrlQ, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2090-2099 FDBCK - CNTRL 1000mb temp 33Q< only')

ax4 = plt.subplot(2,2,4,projection=mapProj)
plt_tls.drawOnGlobe(ax4, tempDiffOverToiEnd1000FdbckCntrl, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2090-2099 FDBCK - CNTRL 1000mb temp')
plt.show()

print('cats')
