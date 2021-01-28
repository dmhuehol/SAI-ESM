import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np

import plotting_tools as plt_tls

filePath = '/Users/dhueholt/Documents/GLENS_data/'
filename = 'b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.001.cam.h0.T.209001-209912.nc'
infile = filePath + filename

glensDataset = xr.open_dataset(infile)

# time, lev, lat, lon
looksee = glensDataset.T[3,3,:,:]
lookC = looksee - 273.15
lookCr = lookC
lookCr[look]

figTime = str(glensDataset.T[3,3,3,3].time.data)
figLevel = str(glensDataset.T[3,3,3,3].lev.data)
titleText = figTime + '    ' + figLevel

# Plotting
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

plt.figure(figsize=(12,2.73*2))
ax = plt.subplot(1,1,1,projection=mapProj) #nrow ncol index
cmap = cmocean.cm.curl
minVal = -39
maxVal = 16
# print(maxVal)
plt_tls.drawOnGlobe(ax, lookC, glensDataset.lat, glensDataset.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
plt.title(titleText)
plt.show()

print("cats")
