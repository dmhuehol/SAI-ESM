from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)
import numpy as np
import xarray as xr

p_out = '/Users/dhueholt/Documents/ecology_fig/20240212_posterAndFinalFigs/'
s = 850
e = 869
name_out = str(s) + str(e) + '_lastma_timeseries.png'
p_lastma = '/Users/dhueholt/Documents/ecology_data/annual_2mTemp/past1000_002_TREFHT_0850-1849_annual.nc'
ds_lastma = xr.open_dataset(p_lastma)
da_lastma = ds_lastma.TREFHT - 273.15
da_lastma.attrs['units'] = '\u00b0C'
lat_weights = np.cos(np.deg2rad(da_lastma.lat))
da_lastma_wght = da_lastma.weighted(lat_weights)
da_lastma_regmn = da_lastma_wght.mean(dim=('lat', 'lon'))
np_lastma = da_lastma_regmn.data
years = da_lastma.time.dt.year

plt.figure()

plt.rcParams.update({'font.family': 'Red Hat Display'})
plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.size': 14})

plt.plot(years, np_lastma, linewidth=2.3, color='#3F6593')
plt.xlim(s, e)
plt.ylim(14,15)
plt.title('CESM2(WACCMma) Last Millennium 2m temp global mean')
plt.ylabel('\u00b0C')
plt.savefig(p_out + name_out, dpi=400)
