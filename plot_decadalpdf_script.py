import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import numpy as np
from glob import glob

import difference_over_time as dot
import plotting_tools as plt_tls
import process_glens_fun as pgf

# Inputs
cntrlPath = '/Users/dhueholt/Documents/GLENS_data/control/'
fdbckPath = '/Users/dhueholt/Documents/GLENS_data/feedback/'
cntrlIn = glob(cntrlPath + '*.nc')
fdbckIn = glob(fdbckPath + '*.nc')

baselinePath = '/Users/dhueholt/Documents/GLENS_data/'
baselineFile = glob(baselinePath + '*.nc')
print(baselineFile)

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210204_pdf/'
saveName = 'pdf_test_'
dpi_val = 300

levOfInt = 1000 #hPa
latOfInt = 34
lonOfInt = -78
quantileForFig = 0.66
regionToPlot = 'global' #aspirational

baselineDs = xr.open_dataset(baselineFile[0])
baselineDsLoi = baselineDs.sel(lev=levOfInt)
dataKey = pgf.discover_data_var(baselineDsLoi)
baselineDarr = baselineDsLoi[dataKey]

baselineDarrMn = baselineDarr.mean(dim='time')
print(np.shape(baselineDarrMn))
bstdev = np.std(baselineDarrMn)
bBinWdth = 0.2*bstdev
bBins = np.arange(int(np.min(baselineDarrMn)), int(np.max(baselineDarrMn)), 0.2)
# print(bBins)
plt_tls.plot_pdf_kdeplot(baselineDarrMn.data, 'slateblue', bBins, 'baseline', savePath, saveName)


print('Completed!')
