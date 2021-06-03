''' plot_zonal_mean
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor
incididunt ut labore et dolore magna aliqua. Ut sem nulla pharetra diam sit amet
 nisl. Sed blandit libero volutpat sed cras ornare arcu dui vivamus.

Written by Daniel Hueholt | June 2021
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
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_solin/'
filenameCntrl = 'control_003_SOLIN_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_SOLIN_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

startInt = [2010,2019]
finalInt = [2090,2099]

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210602_OzoneAndRefinements/'
savePrfx = '1zonmn_FdbckCntrl_'
dpi_val = 400

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

def zonal_mean():

    # Time mean
    glensCntrlTmAvg = glensDarrCntrl.mean(dim='time')
    glensFdbckTmAvg = glensDarrFdbck.mean(dim='time')

    # Zonal mean
    cntrlToPlot = glensCntrlTmAvg.mean(dim='lon')
    fdbckToPlot = glensFdbckTmAvg.mean(dim='lon')

    # Plotting
    yStr = cntrlToPlot.units
    varStr = glensDarrFdbck.long_name
    lats = cntrlToPlot['lat']
    # startStr = str(bndDct['strtYrMtch'])
    # endStr = str(bndDct['endYrMtch'])
    # levStr = pgf.make_level_string(glensCntrlPoi, levOfInt)

    plt.figure()
    plt.plot(lats, cntrlToPlot.data, color='#DF8C20', label='RCP8.5')
    # plt.plot(lats, fdbckToPlot.data, color='#20DFCC', label='SAI')
    plt.legend()
    plt.ylabel(yStr)
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.title('Zonal mean, time mean insolation')
    plt.savefig(savePath + savePrfx + '.png', dpi=dpi_val)
    ic(savePrfx)

zonal_mean()
