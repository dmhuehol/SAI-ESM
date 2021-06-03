''' plot_zonal_mean
Makes figures related to zonal mean values, such as zonal mean, time mean
figures or zonal mean-height contour plots.

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
import cmasher
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
levOfInt = 'stratosphere'

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210603_ncarAndZonal/zonal/'
savePrfx = 'zonmn_FdbckCntrl_'
dpi_val = 400

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

def zonal_mean():
    ''' Plot zonal mean, time mean of a given variable '''
    # Obtain levels
    try:
        glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, levOfInt)
        glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, levOfInt)
        levCheck = 1
    except:
        glensCntrlLoi = glensDarrCntrl
        glensFdbckLoi = glensDarrFdbck
        levCheck = 0

    # Zonal mean
    glensCntrlZnAvg = glensCntrlLoi.mean(dim='lon')
    glensFdbckZnAvg = glensFdbckLoi.mean(dim='lon')

    # Unit conversion
    glensCntrlZnAvg = fcu.molmol_to_ppm(glensCntrlZnAvg)
    glensFdbckZnAvg = fcu.molmol_to_ppm(glensFdbckZnAvg)

    # Time mean
    toiStart = dot.average_over_years(glensCntrlZnAvg, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlZnAvg, finalInt[0], finalInt[1])
    toiEndFdbck = dot.average_over_years(glensFdbckZnAvg, finalInt[0], finalInt[1])

    # Plotting
    yStr = glensCntrlZnAvg.units
    varStr = glensDarrFdbck.long_name
    lats = glensCntrlZnAvg['lat']

    plt.figure()
    plt.plot(lats, toiStart.data, color='#767676', label='Baseline 2010-2019')
    plt.plot(lats, toiEndCntrl.data, color='#DF8C20', label='RCP8.5 2090-2099')
    plt.plot(lats, toiEndFdbck.data, color='#20DFCC', label='SAI 2090-2099')
    plt.legend()
    plt.ylabel(yStr)
    plt.xlabel('deg N')
    plt.autoscale(enable=True, axis='x', tight=True)
    if levCheck:
        plt.title('Zonal mean time mean' + ' ' + str(levOfInt) + ' ' + varStr)
        saveStr = savePath + savePrfx + dataKey + '_' + str(levOfInt)
        plt.savefig(saveStr + '.png', dpi=dpi_val)
    else:
        plt.title('Zonal mean time mean' + ' ' + varStr)
        saveStr = savePath + savePrfx + dataKey
        plt.savefig(saveStr + '.png', dpi=dpi_val)

    ic(savePrfx)
    ic(saveStr)

def zonal_mean_height():
    ''' Plot contours of the zonal and time mean of a variable against height.
    Note that this kind of plot only applies when a variable has a level
    component, and cannot be used for those that have no level dimension.
    '''

    # Zonal mean
    glensCntrlZnAvg = glensDarrCntrl.mean(dim='lon')
    glensFdbckZnAvg = glensDarrFdbck.mean(dim='lon')

    # Unit conversion
    glensCntrlZnAvg = fcu.molmol_to_ppm(glensCntrlZnAvg)
    glensFdbckZnAvg = fcu.molmol_to_ppm(glensFdbckZnAvg)

    # Time mean
    toiStart = dot.average_over_years(glensCntrlZnAvg, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlZnAvg, finalInt[0], finalInt[1])
    toiEndFdbck = dot.average_over_years(glensFdbckZnAvg, finalInt[0], finalInt[1])

    # Calculate
    diffEndCntrlFdbck = toiEndCntrl.data - toiEndFdbck.data

    # Plotting
    savePrfx = 'zonmnhght_FdbckCntrl_2090-2099_'
    yStr = glensCntrlZnAvg['lev'].units
    varStr = glensDarrFdbck.long_name
    lats = glensCntrlZnAvg['lat'].data
    levs = glensCntrlZnAvg['lev'].data

    plt.figure()
    # conSet = plt.contour(lats,levs, toiStart.data,levels=np.arange(0,15,2), colors='black', linewidths=0.75)
    # conSet = plt.contour(lats,levs, toiEndFdbck.data,levels=np.arange(0,15,2), colors='black', linewidths=0.75)
    # conSet = plt.contour(lats,levs, toiEndCntrl.data,levels=np.arange(0,15,2), colors='black', linewidths=0.75)
    conSet = plt.contour(lats,levs, diffEndCntrlFdbck,levels=np.arange(-1,1,0.5), colors='white', linewidths=0.75)
    plt.clabel(conSet,fmt='%1.1f')
    # plt.contourf(lats,levs, toiStart.data, levels=np.arange(0,15), cmap='Oranges')
    plt.contourf(lats,levs, diffEndCntrlFdbck, levels=np.arange(-2,2,0.5), cmap=cmasher.cm.iceburn)
    # plt.contourf(lats,levs, toiEndCntrl.data, levels=np.arange(0,15), cmap='Reds')
    # plt.contourf(lats,levs, toiEndFdbck.data, levels=np.arange(0,15), cmap='Purples')

    # plt.legend()
    plt.ylabel(yStr)
    plt.xlabel('deg N')
    if glensCntrlZnAvg['lev'].positive == 'down':
        plt.gca().invert_yaxis()
    plt.autoscale(enable=True, axis='x', tight=True)

    plt.yscale("log") #Log y-axis
    plt.ylim(1000,0.1)
    plt.gca().set_yticks([1000,100,10,1,0.1])
    plt.gca().set_yticklabels(['1000','100','10','1','0.1'])

    plt.title('RCP8.5 2090-2099 - SAI 2090-2099: Zonal mean time mean' + ' ' + varStr, fontsize=11)
    saveStr = savePath + savePrfx + dataKey
    plt.savefig(saveStr + '.png', dpi=dpi_val)

    ic(savePrfx)
    ic(saveStr)


zonal_mean_height()
