''' plot_difference_globes
Functions to plot differences between RCP8.5 ("Control") and SAI ("Feedback")
scenarios for a GLENS output variable on a 4-panel globe.
plot_basic_difference_globe: 4 panels showing different plots wrt scenario
plot_vertical_difference_globe: 4 panels showing different plots wrt height

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

savePath = '/Users/dhueholt/Documents/GLENS_fig/20210527_threeMoods/globeFunctions/'
savePrfx = 'globe_4p_FdbckCntrl_'
dpi_val = 400

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

def plot_basic_difference_globe():
    ic(dir())
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

def plot_single_basic_difference_globe():
    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, levOfInt)
    glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, levOfInt)

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, finalInt[0], finalInt[1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, finalInt[0], finalInt[1])
    diffToiFdbck =  toiEndCntrl - toiEndFdbck

    # Unit conversion
    diffToiFdbckPlot = fcu.molmol_to_ppm(diffToiFdbck)

    # Plotting
    firstDcd = str(startInt[0]) + '-' + str(startInt[1])
    lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
    sceneStr = 'RCP8.5 - SAI'
    levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)
    varStr = glensDarrCntrl.long_name

    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.delta
    minVal = -diffToiFdbckPlot.quantile(0.99).data
    # minVal = -7 #Override automatic colorbar minimum here
    maxVal = diffToiFdbckPlot.quantile(0.99).data
    # maxVal = 7 #Override automatic colorbar maximum here

    plt_tls.drawOnGlobe(ax, diffToiFdbckPlot, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
    plt.title(lastDcd + ' ' + sceneStr + ' ' + levStr + ' ' + varStr)
    # plt.title("2010-2019 Baseline - 2090-2099 SAI [50 0] ozone") #Override automatic title generation here
    saveStr = savePrfx + dataKey + '_' + levStr + '_' + lastDcd
    # saveStr = 'globe_1p_FdbckCntrl_O3_[50 0]_C2010-2019_F2090-2099'#Override automatic filename generation here
    savename = savePath + saveStr + '.png'
    plt.savefig(savename, dpi=dpi_val, bbox_inches='tight')
    ic(savename)

    print('Completed!')

def plot_vertical_difference_globe():

    # Unit conversion
    glensDarrCntrlUnit = fcu.molmol_to_ppm(glensDarrCntrl)
    glensDarrFdbckUnit = fcu.molmol_to_ppm(glensDarrFdbck)

    # Average over years
    toiEndCntrl = dot.average_over_years(glensDarrCntrlUnit, finalInt[0], finalInt[1])
    toiEndFdbck = dot.average_over_years(glensDarrFdbckUnit, finalInt[0], finalInt[1])

    # Obtain levels
    toiEndCntrlTotal = pgf.obtain_levels(toiEndCntrl, 'total')
    toiEndFdbckTotal = pgf.obtain_levels(toiEndFdbck, 'total')
    toiEndCntrlTrop = pgf.obtain_levels(toiEndCntrl, 'troposphere')
    toiEndFdbckTrop = pgf.obtain_levels(toiEndFdbck, 'troposphere')
    toiEndCntrlStrat = pgf.obtain_levels(toiEndCntrl, 'stratosphere')
    toiEndFdbckStrat = pgf.obtain_levels(toiEndFdbck, 'stratosphere')
    toiEndCntrlLowStrat = pgf.obtain_levels(toiEndCntrl, [250,50])
    toiEndFdbckLowStrat = pgf.obtain_levels(toiEndFdbck, [250,50])

    # Calculate 4-panel values
    diffToiCntrlFdbckTotal = toiEndCntrlTotal - toiEndFdbckTotal
    diffToiCntrlFdbckTrop = toiEndCntrlTrop - toiEndFdbckTrop
    diffToiCntrlFdbckLowStrat = toiEndCntrlLowStrat - toiEndFdbckLowStrat
    diffToiCntrlFdbckStrat = toiEndCntrlStrat - toiEndFdbckStrat

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

    # Make title text
    firstDcd = str(startInt[0]) + '-' + str(startInt[1])
    lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
    cntrlStr = 'RCP8.5'
    fdbckStr = 'SAI'
    # levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)
    varStr = glensDarrCntrl.long_name
    quantileStr = str(quantileOfInt)

    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.delta
    cmapSeq = cmocean.cm.dense
    # minVal = -diffToiCntrl.quantile(0.99).data
    # maxVal = diffToiCntrl.quantile(0.99).data

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensDarrCntrl.lat, glensDarrCntrl.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total ozone concentration')
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ozone concentration')

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title('[250,50] mb ozone concentration')

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckStrat.quantile(0.99), vmax=diffToiCntrlFdbckStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ozone concentration')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    # saveStr = savePrfx + dataKey + '_' + str(levOfInt) + '_' + str(startInt[0]) + str(startInt[1]) + '_' + str(finalInt[0]) + str(finalInt[1])
    saveStr = 'vertical_globe_test'
    savename = savePath + saveStr + '.png'
    plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')
    ic(savename)

    print('Completed!')

def plot_vertical_baseline_difference_globe():

    # Unit conversion
    glensDarrCntrlUnit = fcu.molmol_to_ppm(glensDarrCntrl)
    glensDarrFdbckUnit = fcu.molmol_to_ppm(glensDarrFdbck)

    # Average over years
    toiStart = dot.average_over_years(glensDarrCntrlUnit, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
    toiEndFdbck = dot.average_over_years(glensDarrFdbckUnit, finalInt[0], finalInt[1])

    # Obtain levels
    toiStartCntrlTotal = pgf.obtain_levels(toiStart, 'total')
    toiEndFdbckTotal = pgf.obtain_levels(toiEndFdbck, 'total')
    toiStartCntrlTrop = pgf.obtain_levels(toiStart, 'troposphere')
    toiEndFdbckTrop = pgf.obtain_levels(toiEndFdbck, 'troposphere')
    toiStartCntrlStrat = pgf.obtain_levels(toiStart, 'stratosphere')
    toiEndFdbckStrat = pgf.obtain_levels(toiEndFdbck, 'stratosphere')
    toiStartCntrlLowStrat = pgf.obtain_levels(toiStart, [250,50])
    toiEndFdbckLowStrat = pgf.obtain_levels(toiEndFdbck, [250,50])

    # Calculate 4-panel values
    diffToiCntrlFdbckTotal = toiStartCntrlTotal - toiEndFdbckTotal
    diffToiCntrlFdbckTrop = toiStartCntrlTrop - toiEndFdbckTrop
    diffToiCntrlFdbckLowStrat = toiStartCntrlLowStrat - toiEndFdbckLowStrat
    diffToiCntrlFdbckStrat = toiStartCntrlStrat - toiEndFdbckStrat

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

    # Make title text
    firstDcd = str(startInt[0]) + '-' + str(startInt[1])
    lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
    cntrlStr = 'RCP8.5'
    fdbckStr = 'SAI'
    # levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)
    varStr = glensDarrCntrl.long_name
    quantileStr = str(quantileOfInt)

    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.delta
    cmapSeq = cmocean.cm.dense
    # minVal = -diffToiCntrl.quantile(0.99).data
    # maxVal = diffToiCntrl.quantile(0.99).data

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensDarrCntrl.lat, glensDarrCntrl.lon, cmap, vmin=-15, vmax=15, cbarBool=True, fastBool=True, extent='max')
    plt.title('Total ozone concentration')
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-0.5, vmax=0.5, cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ozone concentration')

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-2, vmax=2, cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title('[250,50] mb ozone concentration')

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-15, vmax=15, cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ozone concentration')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    # saveStr = savePrfx + dataKey + '_' + str(levOfInt) + '_' + str(startInt[0]) + str(startInt[1]) + '_' + str(finalInt[0]) + str(finalInt[1])
    saveStr = 'vertical_baseline_globe_test'
    savename = savePath + saveStr + '.png'
    plt.savefig(savename,dpi=dpi_val,bbox_inches='tight')
    ic(savename)

    print('Completed!')
