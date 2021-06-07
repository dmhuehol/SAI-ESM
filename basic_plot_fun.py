''' basic_plot_fun
Contains the basic plotting functions for GLENS output variables:
difference globes, timeseries, and pdfs.

--DIFFERENCE GLOBES--
Functions to plot differences between RCP8.5 ("Control") and SAI ("Feedback")
scenarios for a GLENS output variable on a 4-panel globe. Equal Earth map
projection used by default.
plot_basic_difference_globe: 4 panels showing different plots wrt scenario
plot_single_basic_difference_globe: 1 panel plot, flexible
plot_vertical_difference_globe: 4 panels showing different plots wrt height (RCP - SAI)
plot_vertical_baseline_difference_globe: 4 panels showing different plots wrt height (BASELINE - SAI)

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

# Example inputs
# dataDict = {
#     "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_o3/',
#     "fnameCntrl": 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc',
#     "fnameFdbck": 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
# }
#
# setDict = {
#     "startIntvl": [2010,2019],
#     "endIntvl": [2090,2099],
#     "levOfInt": 200, #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels,
#     "regOfInt": 'global',  #'global', rlib.Place(), [latN,lonE360]
#     "quantileOfInt": 0.67
# }
#
# outDict = {
#     "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210604_refactoringZonalAndMore/',
#     "dpiVal": 400
# }

## DIFFERENCE GLOBES

def plot_basic_difference_globe(dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe
        (1) change over time for RCP8.5 ("Control")
        (2) change over time for SAI/GEO8.5 ("Feedback")
        (3) difference between RCP8.5 and SAI/GEO8.5 by end interval
        (4) show only where normalized values of (3) are above a given quantile
    '''
    # Open data
    glensDarrCntrl, glensDarrFdbck, dataKey = pgf.open_data(dataDict)

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, setDict["levOfInt"])

    # Unit conversion
    glensCntrlLoi = fcu.molmol_to_ppm(glensCntrlLoi)
    glensFdbckLoi = fcu.molmol_to_ppm(glensFdbckLoi)

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])

    # Calculate 4-panel values
    diffToiCntrl = toiEndCntrl - toiStart
    diffToiFdbck = toiEndFdbck - toiStart
    diffEndCntrlFdbck = toiEndCntrl - toiEndFdbck
    diffEndCntrlFdbckAbsNormQ = pgf.isolate_change_quantile(diffEndCntrlFdbck, setDict["quantileOfInt"])

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

    # Make title text
    firstDcd = str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1])
    lastDcd = str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1])
    cntrlStr = 'RCP8.5'
    fdbckStr = 'SAI'
    levStr = pgf.make_level_string(glensCntrlLoi, setDict["levOfInt"])
    varStr = glensDarrCntrl.long_name
    quantileStr = str(setDict["quantileOfInt"])

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

    savePrfx = 'globe_4p_FdbckCntrl_' #Easy modification for unique filename
    saveStr = savePrfx + dataKey + '_' + str(setDict["levOfInt"]) + '_' + str(setDict["startIntvl"][0]) + str(setDict["startIntvl"][1]) + '_' + str(setDict["endIntvl"][0]) + str(setDict["endIntvl"][1])
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename,dpi=outDict["dpiVal"],bbox_inches='tight')
    ic(savename)

def plot_single_basic_difference_globe(dataDict, setDict, outDict):
    ''' Plot 1 panel difference globe '''
    # Open data
    glensDarrCntrl, glensDarrFdbck, dataKey = pgf.open_data(dataDict)

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, setDict["levOfInt"])

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    diffToiFdbck =  toiEndCntrl - toiEndFdbck

    # Unit conversion
    diffToiFdbckPlot = fcu.molmol_to_ppm(diffToiFdbck)

    # Plotting
    firstDcd = str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1])
    lastDcd = str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1])
    sceneStr = 'RCP8.5 - SAI'
    levStr = pgf.make_level_string(glensCntrlLoi, setDict["levOfInt"])
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

    savePrfx = 'globe_1p_FdbckCntrl_' #Easy modification for unique filename
    saveStr = savePrfx + dataKey + '_' + levStr + '_' + lastDcd
    # saveStr = 'globe_1p_FdbckCntrl_O3_[50 0]_C2010-2019_F2090-2099' #Override automatic filename generation here
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    ic(savename)

def plot_vertical_difference_globe(dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between RCP8.5 and SAI/GEO8.5
    values at ending interval by level
        (1) Total
        (2) Troposphere
        (3) 250mb to 50mb (can be modified for any layer)
        (4) Stratosphere
    '''
    # Open data
    glensDarrCntrl, glensDarrFdbck, dataKey = pgf.open_data(dataDict)

    # Unit conversion
    glensDarrCntrlUnit = fcu.molmol_to_ppm(glensDarrCntrl)
    glensDarrFdbckUnit = fcu.molmol_to_ppm(glensDarrFdbck)

    # Average over years
    toiEndCntrl = dot.average_over_years(glensDarrCntrlUnit, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensDarrFdbckUnit, setDict["endIntvl"][0], setDict["endIntvl"][1])

    # Obtain levels
    toiEndCntrlTotal = pgf.obtain_levels(toiEndCntrl, 'total')
    toiEndFdbckTotal = pgf.obtain_levels(toiEndFdbck, 'total')
    toiEndCntrlTrop = pgf.obtain_levels(toiEndCntrl, 'troposphere')
    toiEndFdbckTrop = pgf.obtain_levels(toiEndFdbck, 'troposphere')
    toiEndCntrlStrat = pgf.obtain_levels(toiEndCntrl, 'stratosphere')
    toiEndFdbckStrat = pgf.obtain_levels(toiEndFdbck, 'stratosphere')
    layerToPlot = [250,50]
    toiEndCntrlLowStrat = pgf.obtain_levels(toiEndCntrl, layerToPlot)
    toiEndFdbckLowStrat = pgf.obtain_levels(toiEndFdbck, layerToPlot)

    # Calculate 4-panel values
    diffToiCntrlFdbckTotal = toiEndCntrlTotal - toiEndFdbckTotal
    diffToiCntrlFdbckTrop = toiEndCntrlTrop - toiEndFdbckTrop
    diffToiCntrlFdbckLowStrat = toiEndCntrlLowStrat - toiEndFdbckLowStrat
    diffToiCntrlFdbckStrat = toiEndCntrlStrat - toiEndFdbckStrat

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

    # Make title text
    firstDcd = str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1])
    lastDcd = str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1])
    cntrlStr = 'RCP8.5'
    fdbckStr = 'SAI'
    varStr = glensDarrCntrl.long_name
    quantileStr = str(setDict["quantileOfInt"])

    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.delta
    cmapSeq = cmocean.cm.dense
    # minVal = -diffToiCntrl.quantile(0.99).data
    # maxVal = diffToiCntrl.quantile(0.99).data

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensDarrCntrl.lat, glensDarrCntrl.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total ' + varStr)
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + varStr)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + varStr)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-diffToiCntrlFdbckStrat.quantile(0.99), vmax=diffToiCntrlFdbckStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + varStr)
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = 'globe_4p_vertical_FdbckCntrl_' #Modify manually for differentiation
    saveStr = savePrfx + dataKey + '_' + str(setDict["endIntvl"][0]) + str(setDict["endIntvl"][1])
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename,dpi=outDict["dpiVal"],bbox_inches='tight')
    ic(savename)

def plot_vertical_baseline_difference_globe(dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between baseline 2010-2019
     and SAI/GEO8.5 values at ending interval by level
        (1) Total
        (2) Troposphere
        (3) Input layer (250mb to 50mb by default)
        (4) Stratosphere
    '''
    # Open data
    glensDarrCntrl, glensDarrFdbck, dataKey = pgf.open_data(dataDict)

    # Unit conversion
    glensDarrCntrlUnit = fcu.molmol_to_ppm(glensDarrCntrl)
    glensDarrFdbckUnit = fcu.molmol_to_ppm(glensDarrFdbck)

    # Average over years
    toiStart = dot.average_over_years(glensDarrCntrlUnit, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndFdbck = dot.average_over_years(glensDarrFdbckUnit, setDict["endIntvl"][0], setDict["endIntvl"][1])

    # Obtain levels
    toiStartCntrlTotal = pgf.obtain_levels(toiStart, 'total')
    toiEndFdbckTotal = pgf.obtain_levels(toiEndFdbck, 'total')
    toiStartCntrlTrop = pgf.obtain_levels(toiStart, 'troposphere')
    toiEndFdbckTrop = pgf.obtain_levels(toiEndFdbck, 'troposphere')
    toiStartCntrlStrat = pgf.obtain_levels(toiStart, 'stratosphere')
    toiEndFdbckStrat = pgf.obtain_levels(toiEndFdbck, 'stratosphere')
    layerToPlot = [250,50]
    toiStartCntrlLowStrat = pgf.obtain_levels(toiStart, layerToPlot)
    toiEndFdbckLowStrat = pgf.obtain_levels(toiEndFdbck, layerToPlot)

    # Calculate 4-panel values
    diffToiCntrlFdbckTotal = toiStartCntrlTotal - toiEndFdbckTotal
    diffToiCntrlFdbckTrop = toiStartCntrlTrop - toiEndFdbckTrop
    diffToiCntrlFdbckLowStrat = toiStartCntrlLowStrat - toiEndFdbckLowStrat
    diffToiCntrlFdbckStrat = toiStartCntrlStrat - toiEndFdbckStrat

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)

    # Make title text
    firstDcd = str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1])
    lastDcd = str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1])
    cntrlStr = 'RCP8.5'
    fdbckStr = 'SAI'
    # levStr = pgf.make_level_string(glensCntrlLoi, setDict["levOfInt"])
    varStr = glensDarrCntrl.long_name
    quantileStr = str(setDict["quantileOfInt"])

    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.delta
    cmapSeq = cmocean.cm.dense
    # minVal = -diffToiCntrl.quantile(0.99).data
    # maxVal = diffToiCntrl.quantile(0.99).data

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensDarrCntrl.lat, glensDarrCntrl.lon, cmap, vmin=-15, vmax=15, cbarBool=True, fastBool=True, extent='max')
    plt.title('Total ' + varStr)
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-0.5, vmax=0.5, cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + varStr)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-2, vmax=2, cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + varStr)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=-15, vmax=15, cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + varStr)
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = 'globe_4p_vertical_baseline_' #Modify manually for differentiation
    saveStr = savePrfx + dataKey + '_' + str(setDict["levOfInt"]) + '_' + str(setDict["startIntvl"][0]) + str(setDict["startIntvl"][1]) + '_' + str(setDict["endIntvl"][0]) + str(setDict["endIntvl"][1])
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename,dpi=outDict["dpiVal"],bbox_inches='tight')
    ic(savename)

## TIMESERIES

def plot_timeseries(dataDict, setDict, outDict):
    ''' Make timeseries of GLENS output variable for RCP8.5 ("Control") and SAI/GEO8.5 ("Feedback") '''
    # Open data
    glensDarrCntrl, glensDarrFdbck, dataKey = pgf.open_data(dataDict)

    bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)
    glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
    glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]

    # Obtain levels
    glensCntrlPoi = pgf.obtain_levels(glensCntrlPoi, setDict["levOfInt"])
    glensFdbckPoi = pgf.obtain_levels(glensFdbckPoi, setDict["levOfInt"])

    # Deal with area
    cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlPoi, setDict["regOfInt"], areaAvgBool=True)
    fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckPoi, setDict["regOfInt"], areaAvgBool=True)

    # Unit conversion
    cntrlToPlot = fcu.molmol_to_ppm(cntrlToPlot)
    fdbckToPlot = fcu.molmol_to_ppm(fdbckToPlot)

    # Plotting
    yStr = cntrlToPlot.units
    varStr = glensDarrFdbck.long_name
    startStr = str(bndDct['strtYrMtch'])
    endStr = str(bndDct['endYrMtch'])
    levStr = pgf.make_level_string(glensCntrlPoi, setDict["levOfInt"])
    ic(levStr, locStr)

    # Make timeseries
    plt.figure()
    plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
    plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
    plt.legend()
    plt.ylabel(yStr)
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.title(varStr + ' ' + levStr + ': ' + startStr + '-' + endStr + ' ' + locTitleStr)

    savePrfx = 'timeseries_' #Modify manually for differentiation
    saveStr = savePrfx + locStr + '_' + levStr
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename,dpi=outDict["dpiVal"],bbox_inches='tight')
    ic(savename)
