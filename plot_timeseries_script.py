import xarray as xr
import difference_over_time as dot
import plotting_tools as plt_tls
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

# Inputs
filenameCntrl = 'control.001.T.r90x45.shift.annual.nc'
filenameFdbck = 'feedback.001.T.r90x45.shift.annual.nc'
dataPath = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck
saveStr = '4p_FdbckControl_T_global_20102020_20202029'
savePath = '/Users/dhueholt/Documents/GLENS_fig/'
savename = savePath + saveStr + '.png'
dpi_val = 300

# Control
glensDatasetCntrl = xr.open_dataset(cntrlPath)
startInt = [2010,2020]
finalInt = [2020,2029]
toiStart = dot.average_over_years(glensDatasetCntrl,startInt[0],startInt[1]) # 2010-2020 is baseline, injection begins 2020
toiEndCntrl = dot.average_over_years(glensDatasetCntrl,finalInt[0],finalInt[1])
tempDiffOverToiCntrl = toiEndCntrl - toiStart
tempDiffOverToi1000Cntrl = tempDiffOverToiCntrl[0,:,:]
# glensEndToi is formatted: lev,lat,lon

# Feedback
glensDatasetFdbck = xr.open_dataset(fdbckPath)
toiEndFdbck = dot.average_over_years(glensDatasetFdbck,finalInt[0],finalInt[1])
tempDiffOverToiFdbck = toiEndFdbck - toiStart
tempDiffOverToi1000Fdbck = tempDiffOverToiFdbck[0,:,:]

# Calculate
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


glensCntrlPoi = glensDatasetCntrl.T[10:90,0,6,37]
glensFdbckPoi = glensDatasetFdbck.T[:,0,6,37]
datasetYears = glensDatasetFdbck['time'].dt.year.data
# print(glensPoi)
# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

plt.figure(figsize=(12,2.73*2))
# ax = plt.subplot(2,2,1) #nrow ncol index
plt.plot(datasetYears,glensCntrlPoi)
plt.plot(datasetYears,glensFdbckPoi)
plt.show()
# cmap = cmocean.cm.curl
# minVal = -8
# maxVal = 8
# # print(maxVal)
# plt_tls.drawOnGlobe(ax, tempDiffOverToi1000Cntrl, glensDatasetCntrl.lat, glensDatasetCntrl.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
# plt.title('2020-2029 - 2010-2020 CONTROL (RCP8.5) 1000mb temp')
#
# ax2 = plt.subplot(2,2,2,projection=mapProj)
# plt_tls.drawOnGlobe(ax2, tempDiffOverToi1000Fdbck, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
# plt.title('2020-2029 - 2010-2020 FEEDBACK (SAI) 1000mb temp')
#
# ax3 = plt.subplot(2,2,3,projection=mapProj)
# plt_tls.drawOnGlobe(ax3, tempDiffOverToiEnd1000FdbckCntrlQ, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
# plt.title('2020-2029 FDBCK - CNTRL 1000mb temp 33Q< only')
#
# ax4 = plt.subplot(2,2,4,projection=mapProj)
# plt_tls.drawOnGlobe(ax4, tempDiffOverToiEnd1000FdbckCntrl, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
# plt.title('2020-2029 FDBCK - CNTRL 1000mb temp')
# plt.show()

plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')

print('Completed!')
