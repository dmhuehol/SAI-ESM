import xarray as xr
import difference_over_time as dot
import plotting_tools as plt_tls
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean

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

CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(1,2,1,projection=mapProj)
cmap = cmocean.cm.curl
minVal = -10 #tempDiff.min()
maxVal = 10 #tempDiff.max()
# print(maxVal)
plt_tls.drawOnGlobe(ax, tempDiffOverToi1000Cntrl, glensDatasetCntrl.lat, glensDatasetCntrl.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2099-2010 difference CONTROL (RCP8.5) 1000mb temp')

ax2 = plt.subplot(1,2,2,projection=mapProj)
plt_tls.drawOnGlobe(ax2, tempDiffOverToi1000Fdbck, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2099-2010 difference FEEDBACK (SAI) 1000mb temp')
plt.show()

print('cats')
