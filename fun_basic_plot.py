''' fun_basic_plot
Contains the basic plotting functions for GLENS/ARISE summary: difference
globes, timeseries, and pdfs. The same three dictionary inputs (defining the
input files, plot settings, and output, respectively) are used by each function.

--DIFFERENCE GLOBES--
Functions to plot "world avoided" scenarios for GLENS and ARISE on a 4-panel
globe. Equal Earth map projection used by default.
plot_basic_difference_globe: 4 panels showing different plots wrt scenario
plot_basic_difference_polar: 4 panels showing different plots wrt scenario in polar map projection
plot_single_basic_difference_globe: 1 panel plot, flexible
CURRENTLY NONFUNCTIONAL: plot_vertical_difference_globe: 4 panels showing different plots wrt height (RCP - SAI)
CURRENTLY NONFUNCTIONAL: plot_vertical_baseline_difference_globe: 4 panels showing different plots wrt height (BASELINE - SAI)

--TIMESERIES--
Make timeseries showing control and SAI scenarios for GLENS and ARISE.
plot_timeseries: plots a timeseries

--PDFs--
Plot pdfs for control and SAI scenarios for GLENS and ARISE.
plot_pdf: plots the pdfs as histogram, step plot, or kde depending on input

Written by Daniel Hueholt | October 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import cftime
import scipy.stats as stats
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn
import cmocean
import cmasher

import fun_process_data as fpd
import fun_plot_tools as fpt
import fun_convert_unit as fcu
import region_library as rlib

## GLOBAL VARIABLES
ensPrp = {
    "dscntntyYrs": [2030],
    "drc": [21,4], #GLENS Control
    "drf": [21,21], #GLENS Feedback
    "drari": [10,10], #ARISE
    "drs245": [5,5] #SSP2-4.5 Control
}

## DIFFERENCE GLOBES

def plot_basic_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe
        (1) change over time for RCP8.5 (GLENS control)
        (2) change over time for SSP2-4.5 (ARISE control)
        (3) diff between RCP8.5 and G1.2(8.5) for end interval (world avoided)
        (4) diff between SSP2-4.5 and G1.5(4.5) for end interval (world avoided)
    '''
    toiStart = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = fpd.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiStart[shrtScn] = toiStartLp
        else:
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
        toiEnd[shrtScn] = toiEndLp

    # Set up panels
    diffToiR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    wrldAvrtdG12R85 = toiEnd['G1.2(8.5)'] - toiEnd['RCP8.5']
    #CESM2-WACCM SSP2-4.5 uses CMIP6 variable names which are often different
    #than in ARISE; subtracting the two DataArrays directly results in nonsense
    wrldAvrtdG15S245 = toiEnd['G1.5(4.5)'].copy() #Duplicate pre-existing array to retain ARISE attributes and shape
    wrldAvrtdG15S245.data = toiEnd['G1.5(4.5)'].data - toiEnd['SSP2-4.5'].data #Subtract the data only
    # scnrsCmprd = toiEnd['G1.2(8.5)'] - toiEnd['G1.5(4.5)'] #Compare ARISE/GLENS CI scenarios USE WITH CAUTION: usually physically meaningless due to differences in model setup!

    panels = (diffToiR85, diffToiS245, wrldAvrtdG12R85, wrldAvrtdG15S245)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance
    # cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    cbVals = [-60,60] #Override automatic colorbar range here
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    plt.suptitle(md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']), fontsize=10)
    # plt.suptitle('2m temperature ens mean', fontsize=10) #Override automatic supertitle here
    lats = rlzList[0].lat
    lons = rlzList[0].lon

    fpt.drawOnGlobe(ax, panels[0], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    fpt.drawOnGlobe(ax2, panels[1], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, panels[2], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    fpt.drawOnGlobe(ax4, panels[3], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_glens_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe
        (1) change over time for RCP8.5 mid-century (GLENS control)
        (2) change over time for G1.2(8.5) mid-century (GLENS feedback)
        (3) diff between RCP8.5 and G1.2(8.5) for mid-century (world avoided)
        (4) diff between RCP8.5 and G1.2(8.5) for end of century (world avoided)
        Input endIntvl as 4 elements, i.e. [2041,2060,2076,2095]
    '''
    toiStart = dict()
    toiMid = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = fpd.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiMidLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][2], setDict["endIntvl"][3])
            toiStart[shrtScn] = toiStartLp
        else:
            toiMidLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][2], setDict["endIntvl"][3])
        toiMid[shrtScn] = toiMidLp
        toiEnd[shrtScn] = toiEndLp

    # Set up panels
    diffToiR85MdCn = toiMid['RCP8.5'] - toiStart['RCP8.5']
    diffToiG12R85MdCn = toiMid['G1.2(8.5)'] - toiStart['RCP8.5'] #Reference for G1.2(8.5) is RCP8.5
    wrldAvrtdG12R85MdCn = toiMid['G1.2(8.5)'] - toiMid['RCP8.5']
    wrldAvrtdG12R85EndCn = toiEnd['G1.2(8.5)'] - toiEnd['RCP8.5']

    panels = (diffToiR85MdCn, diffToiG12R85MdCn, wrldAvrtdG12R85MdCn, wrldAvrtdG12R85EndCn)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    # tropicalPal = seaborn.diverging_palette(133, 324, as_cmap=True)
    cmap = cmocean.cm.balance
    # cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    cbVals = [-70,70] #Override automatic colorbar range here
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    plt.suptitle(md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']), fontsize=10)
    # plt.suptitle('2m temperature ens mean', fontsize=10) #Override automatic supertitle here
    lats = rlzList[0].lat
    lons = rlzList[0].lon

    ic('Watch out! Titles are set manually for all four panels.')
    fpt.drawOnGlobe(ax, panels[0], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('2041-2060 - 2011-2030 RCP8.5')

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    # dimPalette = seaborn.diverging_palette(213,295, as_cmap=True)
    # fpt.drawOnGlobe(ax2, panels[1], lats, lons, dimPalette, vmin=-10, vmax=10, cbarBool=True, fastBool=True, extent='max')
    fpt.drawOnGlobe(ax2, panels[1], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('2041-2060 - 2011-2030 G1.2(8.5)')

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, panels[2], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('G1.2(8.5) - RCP8.5 2041-2060')

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    fpt.drawOnGlobe(ax4, panels[3], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('G1.2(8.5) - RCP8.5 2076-2095')

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + '2076-2095' + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['gfcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_arise_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 2-panel difference globe
        (1) diff between RCP8.5 and G1.2(8.5) for end interval (world avoided)
        (2) diff between SSP2-4.5 and G1.5(4.5) for end interval (world avoided)
    '''
    toiStart = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = fpd.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiStart[shrtScn] = toiStartLp
        else:
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
        toiEnd[shrtScn] = toiEndLp

    # Set up panels
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    #CESM2-WACCM SSP2-4.5 uses CMIP6 variable names which are often different
    #than in ARISE; subtracting the two DataArrays directly results in nonsense
    wrldAvrtdG15S245 = toiEnd['G1.5(4.5)'].copy() #Duplicate pre-existing array to retain ARISE attributes and shape
    wrldAvrtdG15S245.data = toiEnd['G1.5(4.5)'].data - toiEnd['SSP2-4.5'].data #Subtract the data only
    diffToiG15S245 = toiEnd['G1.5(4.5)'].copy()
    diffToiG15S245.data = toiEnd['G1.5(4.5)'].data - toiStart['SSP2-4.5'].data #Subtract the data only
    # scnrsCmprd = toiEnd['G1.2(8.5)'] - toiEnd['G1.5(4.5)'] #Compare ARISE/GLENS CI scenarios USE WITH CAUTION: usually physically meaningless due to differences in model setup!

    panels = (diffToiS245, wrldAvrtdG15S245, diffToiG15S245)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(3,1,1,projection=mapProj) #nrow ncol index
    tropicalPal = seaborn.diverging_palette(324, 133, as_cmap=True)
    cmap = tropicalPal#cmocean.cm.balance
    # cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    cbVals = [-0.1,0.1] #Override automatic colorbar range here
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    plt.suptitle(md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']), fontsize=10)
    # plt.suptitle('2m temperature ens mean', fontsize=10) #Override automatic supertitle here
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    # ic(lats, lons)

    ic('Warning! Titles are set manually')
    fpt.drawOnGlobe(ax, panels[0], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('2041-2060 - 2011-2030 SSP2-4.5')

    ax2 = plt.subplot(3,1,2,projection=mapProj)
    fpt.drawOnGlobe(ax2, panels[1], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('2041-2060 G1.5(4.5) - SSP2-4.5')

    ax3 = plt.subplot(3,1,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, panels[2], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('2041-2060 - 2011-2030 G1.5(4.5)')

    savePrfx = 'arise_'
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_single_basic_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel difference globe '''
    toiStart = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = fpd.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiStart[shrtScn] = toiStartLp
        else:
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
        toiEnd[shrtScn] = toiEndLp

    # Set up panels
    diffToiR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    wrldAvrtdG12R85 = toiEnd['G1.2(8.5)'] - toiEnd['RCP8.5']
    #CESM2-WACCM SSP2-4.5 uses CMIP6 variable names which are often different
    #than in ARISE; subtracting the two DataArrays directly results in nonsense
    wrldAvrtdG15S245 = toiEnd['G1.5(4.5)'].copy() #Duplicate pre-existing array to retain ARISE attributes and shape
    wrldAvrtdG15S245.data = toiEnd['G1.5(4.5)'].data - toiEnd['SSP2-4.5'].data #Subtract the data only
    # scnrsCmprd = toiEnd['G1.2(8.5)'] - toiEnd['G1.5(4.5)'] #Compare ARISE/GLENS CI scenarios USE WITH CAUTION: usually physically meaningless due to differences in model setup!
    blank = toiEnd['G1.5(4.5)'].copy()
    blank.data = toiEnd['G1.5(4.5)'] - toiEnd['G1.5(4.5)']

    panel = blank

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = 'RdBu'#cmocean.cm.curl_r
    # cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    cbVals = [-3,3] #Override automatic colorbar range here
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    lats = rlzList[0].lat
    lons = rlzList[0].lon

    fpt.drawOnGlobe(ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=False, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']' + ' - ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' ' + md['levStr'] + ' ' + md['varStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' ' + md['fdbckStr'] + '-' + md['cntrlStr'] + ' ' + md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']))

    plt.title(" ") #Override automatic title generation here

    savePrfx = '' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    # savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    savename = outDict["savePath"] + 'blankmap.eps'
    # plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.savefig(savename,format='eps')
    plt.close()
    ic(savename)

def plot_vertical_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference globe for difference between two scenario
    values (i.e. reference - SAI[RCP]) at ending interval by level
        (1) Total
        (2) Troposphere
        (3) 250mb to 50mb (can be modified for any layer)
        (4) Stratosphere
    '''
    rlzRfr = rlzList[0]
    rlzFdbck = rlzList[1]
    # Average over years
    rlzToiEnd = list()
    for rc,rDarr in enumerate(rlzList):
        activeRlz = fpd.average_over_years(rDarr, setDict["endIntvl"][0], setDict["endIntvl"][1])
        rlzToiEnd.append(activeRlz)

    # toiEndCntrl = fpd.average_over_years(glensCntrlRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])
    # toiEndFdbck = fpd.average_over_years(glensFdbckRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])

    # Obtain levels
    toiEndCntrlTotal = fpd.obtain_levels(toiEndCntrl, 'total')
    toiEndFdbckTotal = fpd.obtain_levels(toiEndFdbck, 'total')
    toiEndCntrlTrop = fpd.obtain_levels(toiEndCntrl, 'troposphere')
    toiEndFdbckTrop = fpd.obtain_levels(toiEndFdbck, 'troposphere')
    toiEndCntrlStrat = fpd.obtain_levels(toiEndCntrl, 'stratosphere')
    toiEndFdbckStrat = fpd.obtain_levels(toiEndFdbck, 'stratosphere')
    layerToPlot = [250,50]
    toiEndCntrlLowStrat = fpd.obtain_levels(toiEndCntrl, layerToPlot)
    toiEndFdbckLowStrat = fpd.obtain_levels(toiEndFdbck, layerToPlot)

    # Calculate 4-panel values
    diffToiCntrlFdbckTotal = toiEndCntrlTotal - toiEndFdbckTotal
    diffToiCntrlFdbckTrop = toiEndCntrlTrop - toiEndFdbckTrop
    diffToiCntrlFdbckLowStrat = toiEndCntrlLowStrat - toiEndFdbckLowStrat
    diffToiCntrlFdbckStrat = toiEndCntrlStrat - toiEndFdbckStrat

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12,2.73*2))
    md = fpd.meta_book(setDict, dataDict, glensCntrlRlz, labelsToPlot=None)
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][1])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']', fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][0])+']', fontsize=10)
    else:
        plt.suptitle(md['lstDcd'] + ' ' + md['cntrlStr'] + '-' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']), fontsize=10)

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.curl_r

    fpt.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total column ' + md['varStr'])
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    fpt.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + md['varStr'])

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + md['varStr'])

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    fpt.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensCntrlRlz.lat, glensCntrlRlz.lon, cmap, vmin=-diffToiCntrlFdbckStrat.quantile(0.99), vmax=diffToiCntrlFdbckStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
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
    toiStart = fpd.average_over_years(glensCntrlRlz, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndFdbck = fpd.average_over_years(glensFdbckRlz, setDict["endIntvl"][0], setDict["endIntvl"][1])

    # Obtain levels
    toiStartCntrlTotal = fpd.obtain_levels(toiStart, 'total')
    toiEndFdbckTotal = fpd.obtain_levels(toiEndFdbck, 'total')
    toiStartCntrlTrop = fpd.obtain_levels(toiStart, 'troposphere')
    toiEndFdbckTrop = fpd.obtain_levels(toiEndFdbck, 'troposphere')
    toiStartCntrlStrat = fpd.obtain_levels(toiStart, 'stratosphere')
    toiEndFdbckStrat = fpd.obtain_levels(toiEndFdbck, 'stratosphere')
    layerToPlot = [250,50]
    toiStartCntrlLowStrat = fpd.obtain_levels(toiStart, layerToPlot)
    toiEndFdbckLowStrat = fpd.obtain_levels(toiEndFdbck, layerToPlot)

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
    md = fpd.meta_book(setDict, dataDict, glensCntrlRlz, labelsToPlot=None)
    if (setDict["realization"] == 'mean'):
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + '[r'+str(ensPrp['drc'][0])+']' + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + '[r'+str(ensPrp['drf'][1])+']' + ' ', fontsize=10)
    else:
        plt.suptitle(md['frstDcd'] + ' ' + md['cntrlStr'] + ' - ' + md['lstDcd'] + ' ' + md['fdbckStr'] + ' ' + 'Ens: ' + str(setDict['realization']))

    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.curl_r

    fpt.drawOnGlobe(ax, diffToiCntrlFdbckTotal, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckTotal.quantile(0.99), vmax=diffToiCntrlFdbckTotal.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    plt.title('Total column ' + md['varStr'])
    # plt.title('lastDcd + ' - ' + firstDcd + ' ' + cntrlStr + ' ' + 't'otal' + ' ' + varStr)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    fpt.drawOnGlobe(ax2, diffToiCntrlFdbckTrop, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckTrop.quantile(0.99), vmax=diffToiCntrlFdbckTrop.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' - ' + firstDcd + ' ' + fdbckStr + ' ' + 'troposphere' + ' ' + varStr)
    plt.title('Troposphere ' + md['varStr'])

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, diffToiCntrlFdbckLowStrat, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=-diffToiCntrlFdbckLowStrat.quantile(0.99), vmax=diffToiCntrlFdbckLowStrat.quantile(0.99), cbarBool=True, fastBool=True, extent='max')
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + 'stratosphere' + ' ' + varStr)
    plt.title(str(layerToPlot) + ' mb' + ' ' + md['varStr'])

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    v4 = fpt.find_widest_quantile(diffToiCntrlFdbckStrat)
    fpt.drawOnGlobe(ax4, diffToiCntrlFdbckStrat, glensFdbckRlz.lat, glensFdbckRlz.lon, cmap, vmin=v4[0], vmax=v4[1], cbarBool=True, fastBool=True, extent='max')
    plt.title('Stratosphere ' + md['varStr'])
    # plt.title(lastDcd + ' ' + cntrlStr + ' - ' + fdbckStr + ' ' + '[250,50]' + ' ' + varStr)

    savePrfx = ''
    saveStr = savePrfx + md['varSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['vGl'] + '_' + md['glbType']['bGl']
    savename = outDict["savePath"] + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_basic_difference_polar(rlzList, dataDict, setDict, outDict):
    ''' Plot 4-panel difference map with polar projection
        (1) change over time for RCP8.5 (GLENS control)
        (2) change over time for SSP2-4.5 (ARISE control)
        (3) diff between RCP8.5 and G1.2(8.5) for end interval (world avoided)
        (4) diff between SSP2-4.5 and G1.5(4.5) for end interval (world avoided)
    '''
    toiStart = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            toiStartLp = fpd.average_over_years(rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
            toiStart[shrtScn] = toiStartLp
        else:
            toiEndLp = fpd.average_over_years(rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
        toiEnd[shrtScn] = toiEndLp

    # Set up panels
    diffToiR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    wrldAvrtdG12R85 = toiEnd['G1.2(8.5)'] - toiEnd['RCP8.5']
    #CESM2-WACCM SSP2-4.5 uses CMIP6 variable names which are often different
    #than in ARISE; subtracting the two DataArrays directly results in nonsense
    wrldAvrtdG15S245 = toiEnd['G1.5(4.5)'].copy() #Duplicate pre-existing array to retain ARISE attributes and shape
    wrldAvrtdG15S245.data = toiEnd['G1.5(4.5)'].data - toiEnd['SSP2-4.5'].data #Subtract the data only
    # scnrsCmprd = toiEnd['G1.2(8.5)'] - toiEnd['G1.5(4.5)'] #Compare ARISE/GLENS CI scenarios USE WITH CAUTION: usually physically meaningless due to differences in model setup!

    panels = (diffToiR85, diffToiS245, wrldAvrtdG12R85, wrldAvrtdG15S245)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.Orthographic(0, 90)#N: (0,90) S: (180,-90)
    savePrfx = 'NPOLE_'
    # mapProj = cartopy.crs.Orthographic(180, -90)
    # savePrfx = 'SPOLE_'
    plt.figure(figsize=(12,2.73*2))
    ax = plt.subplot(2,2,1,projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.curl_r
    # cbVals = [-panels[0].quantile(0.75).data, panels[0].quantile(0.75).data]
    cbVals = [-1,1] #Override automatic colorbar range here
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    # plt.suptitle(md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']), fontsize=10)
    plt.suptitle('Feb ice thickness ens mean', fontsize=10) #Override automatic supertitle here
    lats = rlzList[0].lat
    lons = rlzList[0].lon

    fpt.drawOnGlobe(ax, panels[0], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['cntrlStr'], fontsize=10)

    ax2 = plt.subplot(2,2,2,projection=mapProj)
    fpt.drawOnGlobe(ax2, panels[1], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)
    else:
        plt.title(md['lstDcd'] + ' - ' + md['frstDcd'] + ' ' + md['s245Cntrl'], fontsize=10)

    ax3 = plt.subplot(2,2,3,projection=mapProj)
    fpt.drawOnGlobe(ax3, panels[2], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['fdbckStr'] + ' - ' + md['cntrlStr'] + ' ' + md['lstDcd'], fontsize=10)

    ax4 = plt.subplot(2,2,4,projection=mapProj)
    fpt.drawOnGlobe(ax4, panels[3], lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    if (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] > ensPrp['dscntntyYrs'][0]):
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)
    elif (setDict["realization"] == 'mean') & (setDict["endIntvl"][0] < ensPrp['dscntntyYrs'][0]):
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)
    else:
        plt.title(md['ariseStr'] + ' - ' + md['s245Cntrl'] + ' ' + md['lstDcd'], fontsize=10)

    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

## TIMESERIES

def plot_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make a simple timeseries of output variable '''
    # Set up data: Isolate time, level, and area of interest
    setYear = [2010, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzToi = rDarr.sel(time=timeSlice)
        rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzToPlot.append(rlzAoi)

    plt.figure()
    md = fpd.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot=None)
    for rpc,rpv in enumerate(rlzToPlot):
        if 'GLENS:Control' in rpv.scenario:
            activeColor = '#D93636'
            activeLabel = md['cntrlStr']
        elif 'GLENS:Feedback' in rpv.scenario:
            activeColor = '#8346C1'
            activeLabel = md['fdbckStr']
        elif 'ARISE:Feedback' in rpv.scenario:
            activeColor = '#12D0B2'
            activeLabel = md['ariseStr']
        elif 'ARISE:Control' in rpv.scenario:
            activeColor = '#F8A53D'
            activeLabel = md['s245Cntrl']
        else:
            sys.exit('Unknown scenario cannot be plotted!')
        yearsOfInt = rpv['time'].dt.year.data #Otherwise the x-axis will be the cftime object, which is ugly
        plt.plot(yearsOfInt, rpv.data, color=activeColor, label=activeLabel)

    b,t = plt.ylim()
    if (setDict["realization"] == 'mean'):
        plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed')
        leg = plt.legend()
        # lText = leg.get_texts()
        # lText[len(lText)-1]._fontproperties = lText[len(lText)-2]._fontproperties.copy() #Change text size for ensemble discontinuity label
        # lText[len(lText)-1].set_fontsize(7)
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
    ''' Plot pdfs for an output variable. Three formats are available: a kernel
    density estimate, a histogram, or a step plot.'''
    rfrncFlag = True #True if plotting any data from before 2020 (during the reference period), False otherwise

    # Set up data
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], setDict["areaAvgBool"])
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
        if 'GLENS:Control' in scnData.scenario:
            cntrlHandlesToPlot = list()
            cntrlHandlesToPlot = fpd.extract_intvl(setDict["cntrlPoi"], periodsOfInt, setDict["timePeriod"], scnData, cntrlHandlesToPlot)
        elif 'GLENS:Feedback' in scnData.scenario:
            fdbckHandlesToPlot = list()
            fdbckHandlesToPlot = fpd.extract_intvl(setDict["fdbckPoi"], periodsOfInt, setDict["timePeriod"], scnData, fdbckHandlesToPlot)
        elif 'ARISE:Feedback' in scnData.scenario:
            ariseHandlesToPlot = list()
            ariseHandlesToPlot = fpd.extract_intvl(setDict["arisePoi"], periodsOfInt, setDict["timePeriod"], scnData, ariseHandlesToPlot)
        elif 'ARISE:Control' in scnData.scenario:
            s245CntrlHandlesToPlot = list()
            s245CntrlHandlesToPlot = fpd.extract_intvl(setDict["s245CntrlPoi"], periodsOfInt, setDict["timePeriod"], scnData, s245CntrlHandlesToPlot)
        else:
            ic(scnData.scenario)
            #No sys.exit(), want to know what the error is if it fails here
    handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot + ariseHandlesToPlot + s245CntrlHandlesToPlot

    # If not applying a spatial average, flatten data so dimensions don't confuse plotting code
    if ~setDict["areaAvgBool"]:
        for ind, h in enumerate(handlesToPlot):
            handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    if rfrncFlag:
        colorsToPlot = fpt.select_colors(rfrncFlag,len(setDict["cntrlPoi"])-1,len(setDict["fdbckPoi"]),len(setDict["arisePoi"]),len(setDict["s245CntrlPoi"]))
    else:
        colorsToPlot = fpt.select_colors(rfrncFlag,len(setDict["cntrlPoi"]),len(setDict["fdbckPoi"]),len(setDict["arisePoi"]),len(setDict["s245CntrlPoi"]))
    labelsToPlot = list()
    labelsToPlot = fpt.generate_labels(labelsToPlot, setDict, ensPrp, rfrncFlag)

    unit = rlzToPlot[0].attrs['units']
    md = fpd.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot)
    titleStr = md['varStr'] + ' ' + md['levStr'] + ' ' + locTitleStr + ' ' + 'Ens ' + str(setDict['realization'])
    labelsToPlot.append(titleStr)
    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['tmStr'] + '_' + locStr + '_' + md['ensStr'] + '_' + md['pid']['pdf'] + '_' + md['pdfStyle'] + '_' + md['spcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    # ic(colorsToPlot) # For troubleshooting

    # Make kde, histograms, or step plots
    try:
        if setDict["plotStyle"] == 'kde':
            fpt.plot_pdf_kdeplot(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, outDict["dpiVal"])
        elif setDict["plotStyle"] == 'hist':
            fpt.plot_pdf_hist(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
        elif setDict["plotStyle"] == 'step':
            fpt.plot_pdf_step(handlesToPlot, colorsToPlot, labelsToPlot, unit, savename, binwidth, outDict["dpiVal"])
        else:
            sys.exit('Invalid plot style')
    except:
        ic('Failed on: ' + savePrfx + saveStr)
        pass
