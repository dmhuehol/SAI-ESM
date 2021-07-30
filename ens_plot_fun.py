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

def plot_ens_spaghetti_timeseries(rlzList, dataDict, setDict, outDict):
    ''' Make timeseries of GLENS output variable for RCP8.5 ("Control") and
    SAI/GEO8.5 ("Feedback"). Ensemble members are visualized in a familiar,
    basic spaghetti plot. '''

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
            if 'GLENS1:Control' in dataToPlot.scenario:
                activeColor = '#DF8C20'
                activeLabel = md['cntrlStr']
            elif 'GLENS1:Feedback' in dataToPlot.scenario:
                activeColor = '#20DFCC'
                activeLabel = md['fdbckStr']
            elif 'GLENS2:Feedback' in dataToPlot.scenario:
                activeColor = '#5C68FF'
                activeLabel = md['fdbckStrG2']
            else:
                sys.exit('Unknown scenario cannot be plotted!')
            yearsOfInt = dataToPlot['time'].dt.year.data #bndDct['mtchYrs']
            if rc==len(rsDarr['realization'].data)-1:
                plt.plot(yearsOfInt,dataToPlot.data,color=activeColor)
            else:
                plt.plot(yearsOfInt,dataToPlot.data,color=activeColor,linewidth=0.3)

    # plt.plot(yearsOfInt,cntrlToPlot.data,color='#DF8C20',label=md['cntrlStr'])
    # plt.plot(yearsOfInt,fdbckToPlot.data,color='#20DFCC',label=md['fdbckStr'])
    b,t = plt.ylim()
    plt.plot([ensPrp["dscntntyYrs"],ensPrp["dscntntyYrs"]],[b,t], color='#36454F', linewidth=0.5, linestyle='dashed', label='RCP8.5 ens: 21 to 2030, 4 to 2095')
    leg = plt.legend()
    # l1,l2,l3 = leg.get_texts()
    # ic(l1,l2,l3) #troubleshooting if the size is changed for the wrong entry
    # l3._fontproperties = l2._fontproperties.copy()
    # l3.set_fontsize(7)
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
