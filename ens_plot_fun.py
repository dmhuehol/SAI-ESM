''' ens_plot_fun
Contains functions to plot data with ensemble visualizations.
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

import process_glens_fun as pgf
import plotting_tools as plt_tls
import fun_convert_unit as fcu
import region_library as rlib

## GLOBAL VARIABLES
ensPrp = {
    "dscntntyYrs": [2030],
    "drc": [21,4], #GLENS Control
    "drf": [21,21], #GLENS Feedback
    "drsci": [10,10], #SCIRIS
    "drs245": [4,5] #SSP2-4.5 Control
}

def plot_ens_spaghetti_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make a simple timeseries of output variable. Ensemble members are
    visualized in a familiar, basic spaghetti plot. '''
    # Set up data: Isolate time, level, and area of interest
    setYear = [2020, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzToi = rDarr.sel(time=timeSlice)
        rlzLoi = pgf.obtain_levels(rlzToi, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = pgf.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzToPlot.append(rlzAoi)

    # Make timeseries
    plt.figure()
    for rsc,rsDarr in enumerate(rlzList):
        for rc in rsDarr['realization'].data:
            dataToPlot, locStr, locTitleStr = pgf.manage_area(rsDarr[rc,:], setDict["regOfInt"], areaAvgBool=True)
            md = pgf.meta_book(setDict, dataDict, dataToPlot, labelsToPlot=None)
            if 'GLENS:Control' in dataToPlot.scenario:
                activeColor = '#D93636'
                activeLabel = md['cntrlStr']
            elif 'GLENS:Feedback' in dataToPlot.scenario:
                activeColor = '#8346C1'
                activeLabel = md['fdbckStr']
            elif 'SCIRIS:Feedback' in dataToPlot.scenario:
                activeColor = '#12D0B2'
                activeLabel = md['scirisStr']
            elif 'SCIRIS:Control' in dataToPlot.scenario:
                activeColor = '#F8A53D'
                activeLabel = md['s245Cntrl']
            else:
                sys.exit('Unknown scenario cannot be plotted!')
            yearsOfInt = dataToPlot['time'].dt.year.data #bndDct['mtchYrs']
            if rc==len(rsDarr['realization'].data)-1:
                plt.plot(yearsOfInt,dataToPlot.data,color=activeColor,label=activeLabel)
            else:
                plt.plot(yearsOfInt,dataToPlot.data,color=activeColor,linewidth=0.3)

    b,t = plt.ylim()
    plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed')
    leg = plt.legend()
    plt.ylabel(md['unit'])
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.xlim(setYear[0],setYear[1])
    plt.title(md['varStr'] + ' ' + md['levStr'] + ': ' + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'Ens ' + str(setDict['realization']))

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['pid']['ts']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_ens_pdf(rlzList, dataDict, setDict, outDict):
    ''' Plot pdfs for an output variable. Three formats are available: a kernel
    density estimate, a histogram, or a step plot. PDFs are made with respect
    to ensemble variability, i.e. the distribution consists of the values given
    by the different ensemble members with or without spatial averaging. Note
    this requires a sizeable ensemble to give useful results if spatial
    averaging is applied!
    '''
    baselineFlag = True #True if plotting any data from before 2020 (during the "Baseline" period), False otherwise

    # Set up data
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = pgf.obtain_levels(rDarr, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = pgf.manage_area(rlzLoi, setDict["regOfInt"], not setDict["dimOfVrblty"]['spcBool'])
        rlzToPlot.append(rlzAoi)

    mnRlzInd = len(rlzToPlot[0]['realization'])-1
    rlzForStats = rlzToPlot[0].sel(realization=mnRlzInd)
    iqr = stats.iqr(rlzForStats,nan_policy='omit')
    binwidth = (2*iqr) / np.power(np.count_nonzero(~np.isnan(rlzForStats)),1/3) # the Freedman-Diaconis rule (NaNs omitted as stackoverflow.com/a/21778195)
    if ~setDict["dimOfVrblty"]['spcBool']:
        binwidth = binwidth/5 #binwidths need to be smaller for small N datasets
    # binwidth = 0.1 #the Let's Not Overthink This rule
    ic(iqr, binwidth)

    # Extract the decades of interest
    handlesToPlot = list()
    for scnData in rlzToPlot:
        periodsOfInt = scnData['time'].dt.year.data
        if 'GLENS:Control' in scnData.scenario:
            cntrlHandlesToPlot = list()
            cntrlHandlesToPlot = pgf.extract_intvl(setDict["cntrlPoi"], periodsOfInt, setDict["timePeriod"], scnData, cntrlHandlesToPlot)
        elif 'GLENS:Feedback' in scnData.scenario:
            fdbckHandlesToPlot = list()
            fdbckHandlesToPlot = pgf.extract_intvl(setDict["fdbckPoi"], periodsOfInt, setDict["timePeriod"], scnData, fdbckHandlesToPlot)
        elif 'SCIRIS:Feedback' in scnData.scenario:
            scirisHandlesToPlot = list()
            scirisHandlesToPlot = pgf.extract_intvl(setDict["scirisPoi"], periodsOfInt, setDict["timePeriod"], scnData, scirisHandlesToPlot)
        elif 'SCIRIS:Control' in scnData.scenario:
            s245CntrlHandlesToPlot = list()
            s245CntrlHandlesToPlot = pgf.extract_intvl(setDict["s245CntrlPoi"], periodsOfInt, setDict["timePeriod"], scnData, s245CntrlHandlesToPlot)
        else:
            ic(scnData.scenario)
            #No sys.exit(), want to know what the error is if it fails here
    handlesToPlot = cntrlHandlesToPlot + fdbckHandlesToPlot + scirisHandlesToPlot + s245CntrlHandlesToPlot
    if not setDict["dimOfVrblty"]['timeBool']:
        handlesToPlot = [h.mean(dim='time') for h in handlesToPlot]
    if not setDict["dimOfVrblty"]['rlzBool']:
        handlesToPlot = [h.mean(dim='realization') for h in handlesToPlot]

    # Flatten data so dimensions don't confuse plotting code
    for ind, h in enumerate(handlesToPlot):
        handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    if baselineFlag:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"])-1,len(setDict["fdbckPoi"]),len(setDict["scirisPoi"]),len(setDict["s245CntrlPoi"]))
    else:
        colorsToPlot = plt_tls.select_colors(baselineFlag,len(setDict["cntrlPoi"]),len(setDict["fdbckPoi"]),len(setDict["scirisPoi"]),len(setDict["s245CntrlPoi"]))
    labelsToPlot = list()
    labelsToPlot = plt_tls.generate_labels(labelsToPlot, setDict, ensPrp, baselineFlag)

    unit = rlzToPlot[0].attrs['units']
    md = pgf.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot)
    titleStr = md['varStr'] + ' ' + md['levStr'] + ' ' + locTitleStr + ' ' + 'PDF:' + pgf.make_dov_title(setDict["dimOfVrblty"])
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
