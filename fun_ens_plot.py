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
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
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
    "drari": [10,10], #ARISE
    "drs245": [4,5] #SSP2-4.5 Control
}

def plot_ens_spaghetti_timeseries(darrList, dataDict, setDict, outDict):
    ''' Make a simple timeseries of output variable. Ensemble members are
    visualized in a familiar, basic spaghetti plot. '''
    setYear = [2010, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0), cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))

    # Plot timeseries
    fig,ax = plt.subplots()
    scnToPlot = list() #Make list of all scenarios to be plotted
    for scnDarr in darrList:
        rlzInScn = scnDarr['realization'].data #Number of realizations in scenario
        scnToPlot.append(scnDarr.scenario) #Add scenario to list
        for rc in rlzInScn:
            ic(rc)
            rlzToi = scnDarr.sel(realization=rc, time=timeSlice) #Single rlz at time of interest
            rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"]) #Level of interest
            rlzToPlot, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=True) #Area of interest
            md = fpd.meta_book(setDict, dataDict, rlzToPlot, labelsToPlot=None) #Extract metadata
            activeColor, activeLabel = fpt.line_from_scenario(rlzToPlot.scenario, md)
            yrsToPlot = rlzToPlot['time'].dt.year.data #bndDct['mtchYrs']
            if rc==len(rlzInScn)-1:
                ic(rc)
                plt.plot(yrsToPlot, rlzToPlot, color=activeColor, label=activeLabel) #Ens mean
            else:
                plt.plot(yrsToPlot, rlzToPlot, color=activeColor, linewidth=0.3) #Individual rlz

    # Plot metadata and settings
    b,t = plt.ylim()
    fpt.plot_metaobjects(scnToPlot, fig)
    leg = plt.legend()
    plt.ylabel(md['unit'])
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.xlim(setYear[0], setYear[1])
    plt.title(md['varStr'] + ' ' + md['levStr'] + ': ' + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'spaghetti')

    # Save image
    savePrfx = ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['spg']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_ens_spread_timeseries(darrList, dataDict, setDict, outDict):
    ''' Make a timeseries of output variable. Ensemble variability is visualized
    as the spread between max and min at each timestep. '''
    setYear = [2010,2095]#[2010, 2095]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))

    # Make timeseries
    if setDict["insetFlag"] == 2:
        plt.rcParams.update({'font.size': 22})
        plt.rcParams.update({'font.family': 'Fira Sans'})
        plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
    fig, ax = plt.subplots()
    scnToPlot = list()
    for darr in darrList:
        darrToi = darr.sel(time=timeSlice)
        darrLoi = fpd.obtain_levels(darrToi, setDict["levOfInt"])
        dataToPlot, locStr, locTitleStr = fpd.manage_area(darrLoi, setDict["regOfInt"], areaAvgBool=True)
        rlzMax = dataToPlot.max(dim='realization')
        rlzMin = dataToPlot.min(dim='realization')
        rlzMn = dataToPlot[len(dataToPlot['realization'])-1] #Last member is ensemble mean
        md = fpd.meta_book(setDict, dataDict, rlzMn, labelsToPlot=None)
        yrsToPlot = rlzMn['time'].dt.year.data #bndDct['mtchYrs']
        scnToPlot.append(darr.scenario)
        activeColor, activeLabel = fpt.line_from_scenario(darr.scenario, md)
        # plt.plot(yearsOfInt,rlzMax.data,color=activeColor,linewidth=0.3) #Border on top of spread
        # plt.plot(yearsOfInt,rlzMin.data,color=activeColor,linewidth=0.3) #Border on bottom of spread
        plt.plot(yrsToPlot,rlzMn.data,color=activeColor,label=activeLabel,linewidth=2)
        ax.fill_between(yrsToPlot, rlzMax.data, rlzMin.data, color=activeColor, alpha=0.3, linewidth=0)

    # Plot metadata and settings
    b,t = plt.ylim()
    # b = 27 #Override automatic b
    # t = 32 #Override automatic t
    fpt.plot_metaobjects(scnToPlot, fig, b, t)
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.xlim(setYear[0],setYear[1])
    if setDict["insetFlag"] == 0: #Standard plot
        plt.ylabel(md['unit'])
        # plt.ylabel('cm')
        leg = plt.legend()
        plt.title(md['varStr'] + ' ' + md['levStr'] + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'spread')
        # plt.title('SST 2010-2095 ' + locTitleStr + ' spread') #Override automatic title generation
        savePrfx = ''
        plt.ylim([b,t])
        # plt.ylim([0.6,0.9]) #Override automatic y-axis limits
    elif setDict["insetFlag"] == 1: #Lines only
        savePrfx = 'INSETQUAL_'
        ax.axis('off')
    else: #Lines and nice axes usually used in posters
        savePrfx = 'INSET_'
        plt.xticks([2015,2040,2065,2090])
        plt.yticks(np.arange(-100, 100, 2))
        # plt.yticks(np.arange(100, 200, 5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylim([b,t])
        # plt.ylabel('cm', fontweight='light')
        # plt.ylabel('\u00B0C', fontweight='light')
        # plt.xlabel('years', fontweight='light')
        # ax.axes.xaxis.set_ticklabels([])
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)

    # Save image
    savePrfx = savePrfx + ''
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['sprd']
    # saveStr = 'SST' + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['sprd']
    savename = outDict["savePath"] + savePrfx + saveStr + '.pdf'
    # fig.set_size_inches(8, 6)
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    # plt.savefig(savename,format='pdf')
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
    rfrncFlag = True #True if plotting any data from before 2020 (during the reference period), False otherwise
    if (dataDict["idGlensCntrl"] is None) and (dataDict["idS245Hist"] is None):
        ic("No control run available--disabling 2011-2030 reference period")
        rfrncFlag = False

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
    glensCntrlHndls = list()
    glensFdbckHndls = list()
    ariseCntrlHndls = list()
    ariseFdbckHndls = list()
    for scnData in rlzToPlot:
        yrsInScn = scnData['time'].dt.year.data
        if 'GLENS:Control' in scnData.scenario:
            glensCntrlHndls = fpd.extract_intvl(setDict["cntrlPoi"], yrsInScn, setDict["timePeriod"], scnData, glensCntrlHndls)
        elif 'GLENS:Feedback' in scnData.scenario:
            glensFdbckHndls = fpd.extract_intvl(setDict["fdbckPoi"], yrsInScn, setDict["timePeriod"], scnData, glensFdbckHndls)
        elif 'ARISE:Control' in scnData.scenario:
            ariseCntrlHndls = fpd.extract_intvl(setDict["s245CntrlPoi"], yrsInScn, setDict["timePeriod"], scnData, ariseCntrlHndls)
        elif 'ARISE:Feedback' in scnData.scenario:
            ariseFdbckHndls = fpd.extract_intvl(setDict["arisePoi"], yrsInScn, setDict["timePeriod"], scnData, ariseFdbckHndls)
        else:
            ic(scnData.scenario)
            #No sys.exit(), want to know what the error is if it fails here
    nGlensCntrlPoi = len(glensCntrlHndls)
    nGlensFdbckPoi = len(glensFdbckHndls)
    nAriseCntrlPoi = len(ariseCntrlHndls)
    nAriseFdbckPoi = len(ariseFdbckHndls)
    handlesToPlot = glensCntrlHndls + glensFdbckHndls + ariseCntrlHndls + ariseFdbckHndls
    if not setDict["dimOfVrblty"]['timeBool']:
        handlesToPlot = [h.mean(dim='time') for h in handlesToPlot]
    if not setDict["dimOfVrblty"]['rlzBool']:
        handlesToPlot = [h.mean(dim='realization') for h in handlesToPlot]

    # Flatten data so dimensions don't confuse plotting code
    for ind, h in enumerate(handlesToPlot):
        handlesToPlot[ind] = h.data.flatten()

    # Generate colors and strings for plots and filenames
    labelsToPlot = list()
    colorsToPlot = list()
    labelsToPlot, colorsToPlot = fpt.generate_labels_colors(labelsToPlot, colorsToPlot, dataDict, setDict, ensPrp, rfrncFlag, nGlensCntrlPoi, nGlensFdbckPoi, nAriseFdbckPoi, nAriseCntrlPoi)
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
