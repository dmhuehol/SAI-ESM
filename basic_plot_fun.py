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
    "drf": [21,21]
}

## DIFFERENCE GLOBES

def plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe
        (1) change over time for RCP8.5 ("Control")
        (2) change over time for SAI/GEO8.5 ("Feedback")
        (3) difference between RCP8.5 and SAI/GEO8.5 by end interval
        (4) show only where normalized values of (3) are above a given quantile
    '''
    # Check ending years before running
    bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    if bndDct["endYrMtch"] < setDict["endIntvl"][0]:
        print("Output ends before end interval, cancelling 4-panel FdbckCntrl globe")
        return

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlRlz, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckRlz, setDict["levOfInt"])

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
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance
    cmapNorm = cmasher.iceburn_r
    cbVals = [-diffToiCntrl.quantile(0.99).data, diffToiCntrl.quantile(0.99).data]
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=glensCntrlLoi, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=glensFdbckLoi)

    plt_tls.drawOnGlobe(ax, diffToiCntrl, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['cntrlStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['cntrlStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    plt_tls.drawOnGlobe(ax2, diffToiFdbck, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drf'][1])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drf'][1])+']' + ' ' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + '[r'+str(ensPrp['drf'][0])+']' + ' - ' + md['frstDcd'] + '[r'+str(ensPrp['drf'][1])+']' + ' ' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    plt_tls.drawOnGlobe(ax3, diffEndCntrlFdbck, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + ' - ' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    plt_tls.drawOnGlobe(ax4, diffEndCntrlFdbckAbsNormQ, glensFdbckRlz.lat, glensFdbckRlz.lon, cmapNorm, vmin=-1., vmax=1., cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' ' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ' + md['levStr'] + ' ' + 'greatest third of norm. change', fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' ' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['levStr'] + ' ' + 'greatest third of norm. change', fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + ' - ' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + 'greatest third of norm. change', fontsize=10)

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
        print("Output ends before end interval, cancelling single FdbckCntrl difference globe")
        return

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlRlz, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckRlz, setDict["levOfInt"])

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    diffToiFdbck =  toiEndCntrl - toiEndFdbck

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance
    cbVals = [-diffToiFdbck.quantile(0.99).data, diffToiFdbck.quantile(0.99).data]
    # cbVals = [-7, 7] #Override automatic colorbar minimum here
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=glensCntrlLoi, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=glensFdbckLoi)

    plt_tls.drawOnGlobe(ax, diffToiFdbck, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '-' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']))


    # plt.title("2010-2019 Baseline - 2090-2099 SAI [50 0] ozone") #Override automatic title generation here

    savePrfx = '' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between RCP8.5 and SAI/GEO8.5
    values at ending interval by level
        (1) Total
        (2) Troposphere
        (3) 250mb to 50mb (can be modified for any layer)
        (4) Stratosphere
    '''
    # Check ending years before running
    bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    if bndDct["endYrMtch"] < setDict["endIntvl"][0]:
        print("Output ends before end interval, cancelling 4-panel vertical difference globe")
        return

    # Average over years
    toiEndCntrl = dot.average_over_years(glensCntrlRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])

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
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=False, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=glensCntrlRlz)
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']', fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']', fontsize=10)
    else:
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '-' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']), fontsize=10)

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance

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

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=False, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=glensCntrlRlz)
    if (setDict["realization"] == 'mean'):
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ', fontsize=10)
    else:
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']))

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance

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
    plt_tls.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckStrat.quantile(0.99), vmax=diffToiCntrlFdbckStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + md['varStr'])
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = ''
    saveStr = savePrfx + md['varSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['vGl'] + '_' + md['glbType']['bGl']
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

## TIMESERIES

def plot_timeseries(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Make timeseries of GLENS output variable for RCP8.5 ("Control") and SAI/GEO8.5 ("Feedback") '''

    setYear = [2020, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    glensCntrlPoi = glensCntrlRlz.sel(time=timeSlice)
    glensFdbckPoi = glensFdbckRlz.sel(time=timeSlice)
    # Uncomment below to automatically pick start/end years
    # bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    # glensCntrlPoi = glensCntrlRlz[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
    # glensFdbckPoi = glensFdbckRlz[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlPoi, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckPoi, setDict["levOfInt"])

    # Deal with area
    cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlLoi, setDict["regOfInt"], areaAvgBool=True)
    fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckLoi, setDict["regOfInt"], areaAvgBool=True)

    # Make timeseries
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=glensCntrlLoi, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=cntrlToPlot)
    plt.figure()
    yearsOfInt = glensCntrlPoi['time'].dt.year.data #bndDct['mtchYrs']
    plt.plot(yearsOfInt,cntrlToPlot.data,color='#DF8C20',label=md['cntrlStr'])
    plt.plot(yearsOfInt,fdbckToPlot.data,color='#20DFCC',label=md['fdbckStr'])
    b,t = plt.ylim()
    if (setDict["realization"] == 'mean'):
        plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed', label='RCP8.5 ens: 21 to 2030, 4 to 2095')
        leg = plt.legend()
        l1,l2,l3 = leg.get_texts()
        # ic(l1,l2,l3) #troubleshooting if the size is changed for the wrong entry
        l3._fontproperties = l2._fontproperties.copy()
        l3.set_fontsize(7)
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

def plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot pdfs for RCP8.5 ("Control") and SAI ("Feedback") values for a GLENS output
    variable. Three formats are available: a kernel density estimate, a histogram,
    or a step plot.'''
    baselineFlag = True #True if plotting any data from before 2020 (during the "Baseline" period), False otherwise

    # Check for valid periods before running
    bndDct = pgf.find_matching_year_bounds(glensCntrlRlz, glensFdbckRlz)
    for poic,poiv in enumerate(setDict["cntrlPoi"]):
        if bndDct["endYrMtch"] < poiv:
            del setDict["cntrlPoi"][poic]

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlRlz, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckRlz, setDict["levOfInt"])

    # Deal with area
    cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlLoi, setDict["regOfInt"], setDict["areaAvgBool"])
    fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckLoi, setDict["regOfInt"], setDict["areaAvgBool"])

    iqr = stats.iqr(cntrlToPlot,nan_policy='omit')
    binwidth = (2*iqr) / np.power(np.count_nonzero(~np.isnan(cntrlToPlot)),1/3) # the Freedman-Diaconis rule (NaNs omitted as stackoverflow.com/a/21778195)
    # binwidth = 0.1 #the Let's Not Overthink This rule
    ic(iqr, binwidth)

    # Extract the decades of interest from the control and feedback datasets
    cntrlYears = cntrlToPlot['time'].dt.year.data
    cntrlHandlesToPlot = list()
    cntrlHandlesToPlot = pgf.extract_doi(setDict["cntrlPoi"], cntrlYears, setDict["timePeriod"], cntrlToPlot, cntrlHandlesToPlot)
    fdbckYears = fdbckToPlot['time'].dt.year.data
    fdbckHandlesToPlot = list()
    fdbckHandlesToPlot = pgf.extract_doi(setDict["fdbckPoi"], fdbckYears, setDict["timePeriod"], fdbckToPlot, fdbckHandlesToPlot)
    handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot

    # If not applying a spatial average, flatten data so dimensions don't confuse plotting code
    if ~setDict["areaAvgBool"]:
        for ind, h in enumerate(handlesToPlot):
            handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    if baselineFlag:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"])-1,len(setDict["fdbckPoi"]))
    else:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"]),len(setDict["fdbckPoi"]))
    labelsToPlot = list()
    labelsToPlot = plt_tls.generate_labels(labelsToPlot, setDict, ensPrp, baselineFlag)

    unit = cntrlToPlot.attrs['units']
    md = pgf.meta_book(setDict, dataDict, labelsToPlot, glensCntrlLoi, glensFdbckRlz, cntrlToPlot)
    titleStr = md['varStr'] + ' ' + md['levStr'] + ' ' + locTitleStr + ' ' + 'Ens ' + str(setDict['realization'])
    labelsToPlot.append(titleStr)
    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['tmStr'] + '_' + locStr + '_' + md['ensStr'] + '_' + md['pid']['pdf'] + '_' + md['pdfStyle'] + '_' + md['spcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    # ic(colorsToPlot) # For troubleshooting

    # Make kde, histograms, or step plots
    if setDict["plotStyle"] == 'kde':
        plt_tls.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, outDict["dpiVal"])
    elif setDict["plotStyle"] == 'hist':
        plt_tls.plot_pdf_hist(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
    elif setDict["plotStyle"] == 'step':
        plt_tls.plot_pdf_step(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
    else:
        sys.exit('Invalid plot style')
