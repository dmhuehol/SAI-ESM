''' basic_plot_fun
Contains the basic plotting functions for GLENS output variables:
difference globes, timeseries, and pdfs. The same three dictionary inputs (defining
the input files, plot settings, and output images, respectively) are used for
each function.

--DIFFERENCE GLOBES--
Functions to plot differences between RCP8.5 ("Control") and SAI/GEO8.5 ("Feedback")
scenarios for a GLENS output variable on a 4-panel globe. Equal Earth map
projection used by default.
plot_basic_difference_globe: 4 panels showing different plots wrt scenario
plot_single_basic_difference_globe: 1 panel plot, flexible
plot_vertical_difference_globe: 4 panels showing different plots wrt height (RCP - SAI)
plot_vertical_baseline_difference_globe: 4 panels showing different plots wrt height (BASELINE - SAI)

--TIMESERIES--
Make timeseries showing progression of both RCP8.5 ("Control") and SAI/GEO8.5
("Feedback") for a GLENS output variable.
plot_timeseries: plots a timeseries

--PDFs--
Plot pdfs for RCP8.5 ("Control") and SAI ("Feedback") values for a GLENS output
variable.
plot_pdf: plots the pdfs as histogram, step plot, or kde depending on input

Written by Daniel Hueholt | July 2021
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
import scipy.stats as stats
import cftime

import difference_over_time as dot
import process_glens_fun as pgf
import plotting_tools as plt_tls
import fun_convert_unit as fcu
import region_library as rlib

## GLOBAL VARIABLES
ensPrp = {
    "dscntntyYrs": [2030],
    "drc": [21,4],
    "drf": [21,21],
    "drg2f": [10,10],
    "drc6c": [4,5]
}

## DIFFERENCE GLOBES

def plot_basic_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe
        (1) change over time for RCP8.5 ("Control")
        (2) change over time for SAI/GEO8.5 ("Feedback")
        (3) difference between RCP8.5 and SAI/GEO8.5 by end interval
        (4) show only where normalized values of (3) are above a given quantile
    '''
    # Check ending years before running
    # bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    # if bndDct["endYrMtch"] < setDict["endIntvl"][0]:
    #     print("Requested interval is not in input! Cancelling 4-panel FdbckCntrl globe")
    #     return

    # Obtain levels
    toiStart = list()
    toiEnd = list()
    trackScn = list()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = pgf.obtain_levels(rDarr, setDict["levOfInt"])
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = dot.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiEndLp = dot.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiStart.append(toiStartLp)
        else:
            toiEndLp = dot.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
        trackScn.append(rlzLoi.attrs['scenario'])
        toiEnd.append(toiEndLp)

    # Calculate 4-panel values
    diffToiR85 = toiEnd[0] - toiStart[0]
    diffToiS245 = toiEnd[3] - toiStart[1]
    wrldAvrtdG12R85 = toiEnd[1] - toiEnd[0]
    wrldAvrtdG15S245 = toiEnd[2].data - toiEnd[3].data #Account for subtle format differences between CESM2-WACCM SSP2-4.5 and SCIRIS runs
    wrldAvrtdG15S245.attrs = toiEnd[2].attrs
    # scnrsCmprd = toiEnd[2] - toiEnd[1]

    panels = (diffToiR85, diffToiS245, wrldAvrtdG12R85, wrldAvrtdG15S245)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance
    cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    md = pgf.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    plt.suptitle(md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']), fontsize=10)
    lats = rlzList[0].lat
    lons = rlzList[1].lon

    plt_tls.drawOnGlobe(ax, panels[0], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['cntrlStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['cntrlStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, panels[1], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc6c'][1])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc6c'][0])+']' + ' ' + md['ssp245Str'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc6c'][0])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc6c'][0])+']' + ' ' + md['ssp245Str'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['ssp245Str'], fontsize=10)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, panels[2], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, panels[3], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStrG2'] + '[r'+str(ensPrp['drg2f'][1])+']' + ' - ' + md['ssp245Str'] + '[r'+str(ensPrp['drc'][1])+']' + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStrG2'] + '[r'+str(ensPrp['drg2f'][0])+']' + ' - ' + md['ssp245Str'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['fdbckStrG2'] + ' - ' + md['ssp245Str'] + ' ' + md['lstDcd'], fontsize=10)

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot 1 panel difference globe '''
    # Check ending years before running
    bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    if bndDct["endYrMtch"] < setDict["endIntvl"][0]:
        print("Requested interval is not in input! Cancelling single FdbckCntrl difference globe")
        return

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlRlz, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckRlz, setDict["levOfInt"])

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    diffToiFdbck = toiEndFdbck - toiEndCntrl

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = 'RdBu'#cmocean.cm.curl_r
    cbVals = plt_tls.find_widest_quantile(diffToiFdbck)
    # cbVals = [-7, 7] #Override automatic colorbar minimum here
    md = pgf.meta_book(setDict, dataDict, glensFdbckLoi, labelsToPlot=None)

    plt_tls.drawOnGlobe(ax, diffToiFdbck, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '-' + md['cntrlStr'] + ' ' + md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']))

    # plt.title("2010-2019 Baseline - 2090-2099 SAI [50 0] ozone") #Override automatic title generation here

    savePrfx = '' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_vertical_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between two scenario
    values (i.e. baseline - SAI[RCP]) at ending interval by level
        (1) Total
        (2) Troposphere
        (3) 250mb to 50mb (can be modified for any layer)
        (4) Stratosphere
    '''
    # Check ending years before running
    # bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    # if bndDct["endYrMtch"] < setDict["endIntvl"][0]:
    #     print("Requested interval is not in input! Cancelling 4-panel vertical difference globe.")
    #     return
    rlzRfr = rlzList[0]
    rlzFdbck = rlzList[1]
    # Average over years
    rlzToiEnd = list()
    for rc,rDarr in enumerate(rlzList):
        activeRlz = dot.average_over_years(rDarr, setDict["endIntvl"][0], setDict["endIntvl"][1])
        rlzToiEnd.append(activeRlz)

    # toiEndCntrl = dot.average_over_years(glensCntrlRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])
    # toiEndFdbck = dot.average_over_years(glensFdbckRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])

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
    plt.figure(figsize=(12,2.73*2))
    md = pgf.meta_book(setDict, dataDict, glensCntrlRlz, labelsToPlot=None)
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']', fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']', fontsize=10)
    else:
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '-' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']), fontsize=10)

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.curl_r

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total column ' + md['varStr'])
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + md['varStr'])

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + md['varStr'])

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckStrat.quantile(0.99), vmax=diffToiCntrlFdbckStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + md['varStr'])
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['vGl'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between baseline 2010-2019
     and SAI/GEO8.5 values at ending interval by level
        (1) Total
        (2) Troposphere
        (3) Input layer (250mb to 50mb by default)
        (4) Stratosphere
    '''

    # Average over years
    toiStart = dot.average_over_years(glensCntrlRlz, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndFdbck = dot.average_over_years(glensFdbckRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])

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

    # ic(diffToiCntrlFdbckStrat)
    # ic(np.min(diffToiCntrlFdbckStrat))
    # ic(np.max(diffToiCntrlFdbckStrat))
    # ic(diffToiCntrlFdbckStrat.quantile(0.01))
    # ic(-diffToiCntrlFdbckStrat.quantile(0.99))
    # sys.exit('STOP')

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    md = pgf.meta_book(setDict, dataDict, glensCntrlRlz, labelsToPlot=None)
    if (setDict["realization"] == 'mean'):
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ', fontsize=10)
    else:
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']))

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.curl_r

    plt_tls.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total column ' + md['varStr'])
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + md['varStr'])

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + md['varStr'])

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    v4 = plt_tls.find_widest_quantile(diffToiCntrlFdbckStrat)
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=v4[0], vmax=v4[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + md['varStr'])
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = ''
    saveStr = savePrfx + md['varSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['vGl'] + '_' + md['glbType']['bGl']
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

## TIMESERIES

def plot_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make timeseries of GLENS output variable for RCP8.5 ("Control") and SAI/GEO8.5 ("Feedback") '''
    # Set up data: Isolate time, level, and area of interest
    setYear = [2020, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzToi = rDarr.sel(time=timeSlice)
        rlzLoi = pgf.obtain_levels(rlzToi, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = pgf.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzToPlot.append(rlzAoi)

    plt.figure()
    md = pgf.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot=None)
    for rpc,rpv in enumerate(rlzToPlot):
        if 'GLENS1:Control' in rpv.scenario:
            activeColor = '#DF8C20'
            activeLabel = md['cntrlStr']
        elif 'GLENS1:Feedback' in rpv.scenario:
            activeColor = '#20DFCC'
            activeLabel = md['fdbckStr']
        elif 'GLENS2:Feedback' in rpv.scenario:
            activeColor = '#5C68FF'
            activeLabel = md['fdbckStrG2']
        elif 'CMIP6:Control' in rpv.scenario:
            activeColor = '#FED8B1'
            activeLabel = md['ssp245Str']
        else:
            sys.exit('Unknown scenario cannot be plotted!')
        yearsOfInt = rpv['time'].dt.year.data #Otherwise the x-axis will be the cftime object, which is ugly
        plt.plot(yearsOfInt, rpv.data, color=activeColor, label=activeLabel)

    b,t = plt.ylim()
    if (setDict["realization"] == 'mean'):
        plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed', label='RCP8.5 ens: 21 to 2030, 4 to 2095')
        leg = plt.legend()
        lText = leg.get_texts()
        # ic(l1,l2,l3) #troubleshooting if the size is changed for the wrong entry
        lText[len(lText)-1]._fontproperties = lText[len(lText)-2]._fontproperties.copy()
        lText[len(lText)-1].set_fontsize(7)
    else:
        plt.legend()
    plt.ylabel(md['unit'])
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.xlim(setYear[0],setYear[1])
    plt.title(md['varStr'] + ' ' + md['levStr'] + ': ' + md['strtStr'] + '-' + md['endStr'] + ' ' + locTitleStr  + ' ' + 'Ens ' + str(setDict['realization']))

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['strtStr'] + md['endStr'] + '_' + locStr + '_' + md['ensStr'] + '_' + md['pid']['ts']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

## PDFs

def plot_pdf(rlzList, dataDict, setDict, outDict):
    ''' Plot pdfs for RCP8.5 ("Control") and SAI ("Feedback") values for a GLENS output
    variable. Three formats are available: a kernel density estimate, a histogram,
    or a step plot.'''
    baselineFlag = True #True if plotting any data from before 2020 (during the "Baseline" period), False otherwise

    # Check for valid periods before running
    # bndDct = pgf.find_matching_year_bounds(rlzList)
    # for poic,poiv in enumerate(setDict["cntrlPoi"]):
    #     if bndDct["endYrMtch"] < poiv:
    #         del setDict["cntrlPoi"][poic]

    # Set up data
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = pgf.obtain_levels(rDarr, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = pgf.manage_area(rlzLoi, setDict["regOfInt"], setDict["areaAvgBool"])
        rlzToPlot.append(rlzAoi)

    iqr = stats.iqr(rlzToPlot[0],nan_policy='omit')
    binwidth = (2*iqr) / np.power(np.count_nonzero(~np.isnan(rlzToPlot[0])),1/3) # the Freedman-Diaconis rule (NaNs omitted as stackoverflow.com/a/21778195)
    if setDict["areaAvgBool"]:
        binwidth = binwidth/5 #binwidths need to be smaller for small N datasets
    # binwidth = 0.1 #the Let's Not Overthink This rule
    ic(iqr, binwidth)

    # Extract the decades of interest
    handlesToPlot = list()
    for scnData in rlzToPlot:
        periodsOfInt = scnData['time'].dt.year.data
        if 'GLENS1:Control' in scnData.scenario:
            cntrlHandlesToPlot = list()
            cntrlHandlesToPlot = pgf.extract_doi(setDict["cntrlPoi"], periodsOfInt, setDict["timePeriod"], scnData, cntrlHandlesToPlot)
        elif 'GLENS1:Feedback' in scnData.scenario:
            fdbckHandlesToPlot = list()
            fdbckHandlesToPlot = pgf.extract_doi(setDict["fdbckPoi"], periodsOfInt, setDict["timePeriod"], scnData, fdbckHandlesToPlot)
        elif 'GLENS2:Feedback' in scnData.scenario:
            glens2HandlesToPlot = list()
            glens2HandlesToPlot = pgf.extract_doi(setDict["glens2Poi"], periodsOfInt, setDict["timePeriod"], scnData, glens2HandlesToPlot)
        else:
            ic(scnData.scenario)
            #No sys.exit(), want to know what the error is if it fails here
    handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot + glens2HandlesToPlot

    # If not applying a spatial average, flatten data so dimensions don't confuse plotting code
    if ~setDict["areaAvgBool"]:
        for ind, h in enumerate(handlesToPlot):
            handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    if baselineFlag:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"])-1,len(setDict["fdbckPoi"]),len(setDict["glens2Poi"]))
    else:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"]),len(setDict["fdbckPoi"]),len(setDict["glens2Poi"]))
    labelsToPlot = list()
    labelsToPlot = plt_tls.generate_labels(labelsToPlot, setDict, ensPrp, baselineFlag)

    unit = rlzToPlot[0].attrs['units']
    md = pgf.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot)
    titleStr = md['varStr'] + ' ' + md['levStr'] + ' ' + locTitleStr + ' ' + 'Ens ' + str(setDict['realization'])
    labelsToPlot.append(titleStr)
    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['tmStr'] + '_' + locStr + '_' + md['ensStr'] + '_' + md['pid']['pdf'] + '_' + md['pdfStyle'] + '_' + md['spcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    # ic(colorsToPlot) # For troubleshooting

    # Make kde, histograms, or step plots
    try:
        if setDict["plotStyle"] == 'kde':
            plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, outDict["dpiVal"])
        elif setDict["plotStyle"] == 'hist':
            plt_tls.plot_pdf_hist(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
        elif setDict["plotStyle"] == 'step':
            plt_tls.plot_pdf_step(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
        else:
            sys.exit('Invalid plot style')
    except:
        ic('Failed on: ' + savePrfx + saveStr)
        pass
