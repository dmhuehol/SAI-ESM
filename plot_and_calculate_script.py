import xarray as xr
import difference_over_time as dot
import plotting_tools as plt_tls
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean

filenameCntrl = 'control.001.T.r90x45.shift.annual.nc'
glensDatasetCntrl, tempDiffCntrl = dot.simple_diff_calc(filenameCntrl)
tempDiff1000Cntrl = tempDiffCntrl[0,:,:]

filenameFdbck = 'feedback.001.T.r90x45.shift.annual.nc'
glensDatasetFdbck, tempDiffFdbck = dot.simple_diff_calc(filenameFdbck)
tempDiff1000Fdbck = tempDiffFdbck[0,:,:]

CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(1,2,1,projection=mapProj)
cmap = cmocean.cm.curl
minVal = -10 #tempDiff.min()
maxVal = 10 #tempDiff.max()
# print(maxVal)
plt_tls.drawOnGlobe(ax, tempDiff1000Cntrl, glensDatasetCntrl.lat, glensDatasetCntrl.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2099-2010 difference CONTROL (RCP8.5) 1000mb temp')

ax2 = plt.subplot(1,2,2,projection=mapProj)
plt_tls.drawOnGlobe(ax2, tempDiff1000Fdbck, glensDatasetFdbck.lat, glensDatasetFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title('2099-2010 difference FEEDBACK (SAI) 1000mb temp')
plt.show()

print('cats')
