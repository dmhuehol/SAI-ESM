''' fun_ens_plot
Contains functions to plot data with ensemble visualizations.

--TIMESERIES--
plot_ens_spaghetti_timeseries: plots a timeseries with each member visualized
as its own line, this is the classic "spaghetti plot"
plot_ens_spread_timeseries: plots timeseries with ensemble mean shown as line
and spread between max and min member shaded

--PDFs--
plot_ens_pdf: plot pdfs with different realizations as a degree of freedom

Written by Daniel Hueholt
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

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  # the location of the font file
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import fun_process_data as fpd
import fun_plot_tools as fpt
import fun_convert_unit as fcu
import region_library as rlib

## GLOBAL VARIABLES
ensPrp = {
    "dscntntyYrs": [2030],
    "drc": [21,4], #GLENS Control
    "drf": [21,21], #GLENS Feedback
    "drsci": [10,10], #ARISE
    "drs245": [4,5] #SSP2-4.5 Control
}

def plot_ens_spaghetti_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make a simple timeseries of output variable. Ensemble members are
    visualized in a familiar, basic spaghetti plot. '''
    # Set up data: Isolate time, level, and area of interest
    setYear = [2010, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzToi = rDarr.sel(time=timeSlice)
        rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzToPlot.append(rlzAoi)

    # Make timeseries
    plt.figure()
    for rsc,rsDarr in enumerate(rlzList):
        for rc in rsDarr['realization'].data:
            dataToPlot, locStr, locTitleStr = fpd.manage_area(rsDarr.sel(realization=rc), setDict["regOfInt"], areaAvgBool=True)
            md = fpd.meta_book(setDict, dataDict, dataToPlot, labelsToPlot=None)
            if 'GLENS:Control' in dataToPlot.scenario:
                activeColor = '#D93636'
                activeLabel = md['cntrlStr']
            elif 'GLENS:Feedback' in dataToPlot.scenario:
                activeColor = '#8346C1'
                activeLabel = md['fdbckStr']
            elif 'ARISE:Feedback' in dataToPlot.scenario:
                activeColor = '#12D0B2'
                activeLabel = md['ariseStr']
            elif 'ARISE:Control' in dataToPlot.scenario:
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
    plt.title(md['varStr'] + ' ' + md['levStr'] + ': ' + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'spaghetti')

    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['spg']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_ens_spread_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make a timeseries of output variable. Ensemble variability is visualized
    as the spread between max and min at each timestep. '''
    # Set up data: Isolate time, level, and area of interest
    setYear = [2010,2095]#[2010, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    rlzToPlot = list()
    for rc,rDarr in enumerate(rlzList):
        rlzToi = rDarr.sel(time=timeSlice)
        rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzToPlot.append(rlzAoi)

    # Make timeseries
    if setDict["insetFlag"] == 2:
        plt.rcParams.update({'font.size': 18})
        plt.rcParams.update({'font.family': 'Fira Sans'})
        plt.rcParams.update({'font.weight': 'bold'})
    fig, ax = plt.subplots()
    for rsc,rsDarr in enumerate(rlzList):
        dataAoi, locStr, locTitleStr = fpd.manage_area(rsDarr, setDict["regOfInt"], areaAvgBool=True)
        rlzMax = dataAoi.max(dim='realization')
        rlzMin = dataAoi.min(dim='realization')
        rlzMn = dataAoi[len(dataAoi['realization'])-1] #last member is ensemble mean
        md = fpd.meta_book(setDict, dataDict, rlzMn, labelsToPlot=None)
        yearsOfInt = rlzMn['time'].dt.year.data #bndDct['mtchYrs']
        if 'GLENS:Control' in rsDarr.scenario:
            activeColor = '#D93636'
            activeLabel = md['cntrlStr']
            # plt.plot(yearsOfInt[0:20],rlzMn.data[0:20],color=activeColor,label=activeLabel)
            # plt.plot(yearsOfInt[20:],rlzMn.data[20:],color=activeColor,linestyle='dotted')
        elif 'GLENS:Feedback' in rsDarr.scenario:
            activeColor = '#8346C1'
            activeLabel = md['fdbckStr']
            # plt.plot(yearsOfInt,rlzMn.data,color=activeColor,label=activeLabel)
        elif 'ARISE:Feedback' in rsDarr.scenario:
            activeColor = '#12D0B2'
            activeLabel = md['ariseStr']
            # plt.plot(yearsOfInt,rlzMn.data,color=activeColor,label=activeLabel)
        elif 'ARISE:Control' in rsDarr.scenario:
            activeColor = '#F8A53D'
            activeLabel = md['s245Cntrl']
            # plt.plot(yearsOfInt[0:175],rlzMn.data[0:175],color=activeColor,linestyle='dotted')
            # plt.plot(yearsOfInt[175:],rlzMn.data[175:],color=activeColor,label=activeLabel)
        else:
            sys.exit('Unknown scenario cannot be plotted!')

        # plt.plot(yearsOfInt,rlzMax.data,color=activeColor,linewidth=0.3)
        # plt.plot(yearsOfInt,rlzMin.data,color=activeColor,linewidth=0.3)
        plt.plot(yearsOfInt,rlzMn.data,color=activeColor,label=activeLabel,linewidth=2)
        ax.fill_between(yearsOfInt, rlzMax.data, rlzMin.data, color=activeColor, alpha=0.3, linewidth=0)

    b,t = plt.ylim()
    # b = 0 #Override automatic b
    # t = 0.45 #Override automatic t
    # plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed')
    # plt.plot(2015,b+(abs(b-t))*0.01,color='#F8A53D',marker='v')
    plt.plot(2015,3.07,color='#F8A53D',marker='v')
    # plt.plot(2030,b+(abs(b-t))*0.01,color='#D93636',marker='v')
    plt.plot(2030,3.07,color='#D93636',marker='v')
    plt.plot([2020,2020],[b,t], color='#8346C1', linewidth=0.7, linestyle='dashed')
    plt.plot([2035,2035],[b,t], color='#12D0B2', linewidth=0.7, linestyle='dashed')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.xlim(setYear[0],setYear[1])
    if setDict["insetFlag"] == 0:
        plt.ylabel(md['unit'])
        leg = plt.legend()
        plt.title(md['varStr'] + ' ' + md['levStr'] + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'spread')
        # plt.title('SST 2010-2095 ' + locTitleStr + ' spread') #Override automatic title generation
        savePrfx = ''
        # plt.ylim([0.6,0.9])
        plt.ylim([b,t])
    elif setDict["insetFlag"] == 1:
        savePrfx = 'INSETQUAL_'
        ax.axis('off')
    else:
        savePrfx = 'INSET_'
        plt.xticks([2010,2030,2050,2070,2090])
        # plt.yticks(np.arange(-30,100,2))
        plt.yticks(np.arange(0,100,10)) #Ice thickness
        plt.yticks(np.arange(0,100,2)) #Air temperature
        plt.yticks(np.arange(100,300,30)) #Monsoon precip
        plt.yticks(np.arange(25,33,2)) #Tropical SST
        plt.yticks(np.arange(10,50,1)) #Global SST
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # plt.ylim([b,t])
        # plt.ylim([10,30])
        plt.ylim([18.5,22.5])
        plt.ylabel('\u00B0C')
        # plt.ylabel('cm')
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)

    savePrfx = savePrfx + ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['sprd']
    # saveStr = 'SST' + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['sprd']
    savename = outDict["savePath"] + savePrfx + saveStr + '.pdf'
    # plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.savefig(savename,format='pdf')
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
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        rlzAoi, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], not setDict["dimOfVrblty"]['spcBool'])
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
    if not setDict["dimOfVrblty"]['timeBool']:
        handlesToPlot = [h.mean(dim='time') for h in handlesToPlot]
    if not setDict["dimOfVrblty"]['rlzBool']:
        handlesToPlot = [h.mean(dim='realization') for h in handlesToPlot]

    # Flatten data so dimensions don't confuse plotting code
    for ind, h in enumerate(handlesToPlot):
        handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    if baselineFlag:
        colorsToPlot = fpt.select_colors(baselineFlag,len(setDict["cntrlPoi"])-1,len(setDict["fdbckPoi"]),len(setDict["arisePoi"]),len(setDict["s245CntrlPoi"]))
    else:
        colorsToPlot = fpt.select_colors(baselineFlag,len(setDict["cntrlPoi"]),len(setDict["fdbckPoi"]),len(setDict["arisePoi"]),len(setDict["s245CntrlPoi"]))
    labelsToPlot = list()
    labelsToPlot = fpt.generate_labels(labelsToPlot, setDict, ensPrp, baselineFlag)

    unit = rlzToPlot[0].attrs['units']
    md = fpd.meta_book(setDict, dataDict, rlzToPlot[0], labelsToPlot)
    titleStr = md['varStr'] + ' ' + md['levStr'] + ' ' + locTitleStr + ' ' + 'PDF:' + fpd.make_dov_title(setDict["dimOfVrblty"])
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
