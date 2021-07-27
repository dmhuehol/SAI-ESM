''' the_arctic_ifrac_story
Makes Arctic ice timeseries with maximum aesthetics, used on Barnes Group
website in 2021
Written by Daniel Hueholt | July 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import xarray as xr
xr.set_options(keep_attrs=True)
import matplotlib
import matplotlib.font_manager as fm
matplotlib.font_manager._rebuild()
fontPath = '/Users/dhueholt/Library/Fonts/'  # the location of the font file
firaSansPath = {"Thin": 'FiraSans-Thin.ttf',
                "Medium": 'FiraSans-Medium.ttf'}
robotoPath = {"Regular": 'Roboto-Regular.ttf'}
Roboto = fm.FontProperties(fname=fontPath+robotoPath['Regular'])  # get the font based on the font_path
FiraSansMed = fm.FontProperties(fname=fontPath+firaSansPath['Medium'])
FiraSansThin = fm.FontProperties(fname=fontPath+firaSansPath['Thin'])
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy
import cartopy.crs as ccrs
import numpy as np
import cftime

import difference_over_time as dot
import process_glens_fun as pgf
import plotting_tools as plt_tls
import fun_convert_unit as fcu
import region_library as rlib

# Settings for wrap_basicplots_script
# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/sept_IFRAC/emrg/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*'
}
setDict = {
    "convert": None #tuple of converter(s), or None if using default units
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210726_newfigAndG2/',
    "dpiVal": 400
}
loopDict = {
    "realizations": ('mean',),
    "levels": (None,),
    "regions": (rlib.Arctic(),),
    "aaBools": (True,)
}

def plot_zchfd_timeseries_sept_arc_ifrac(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Make timeseries of GLENS output variable for RCP8.5 ("Control") and SAI/GEO8.5 ("Feedback") '''

    setYear = [2020, 2096]
    timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    glensCntrlPoi = glensCntrlRlz.sel(time=timeSlice)
    glensFdbckPoi = glensFdbckRlz.sel(time=timeSlice)

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlPoi, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckPoi, setDict["levOfInt"])

    # Deal with area
    cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlLoi, setDict["regOfInt"], areaAvgBool=True)
    fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckLoi, setDict["regOfInt"], areaAvgBool=True)

    # Make timeseries
    fig,ax = plt.subplots()
    yearsOfInt = glensCntrlPoi['time'].dt.year.data
    plt.plot(yearsOfInt,cntrlToPlot.data,color='#DF8C20')
    plt.plot(yearsOfInt,fdbckToPlot.data,color='#20DFCC')
    plt.ylabel('Ice fraction',fontproperties=FiraSansThin)
    for label in ax.get_xticklabels():
        label.set_fontproperties(FiraSansThin)
    for label in ax.get_yticklabels():
        label.set_fontproperties(FiraSansThin)
    ax.annotate('September Arctic sea ice fraction',(0.02,0.935),(0.02,0.935),'axes fraction',fontproperties=FiraSansMed,fontsize=16,color='#000000')
    # ax.annotate('RCP8.5',(0.08,0.38),(0.08,0.38),'axes fraction',fontproperties=FiraSansMed,fontsize=14,color='#DF8C20')
    # ax.annotate('SAI',(0.43,0.82),(0.43,0.82),'axes fraction',fontproperties=FiraSansMed,fontsize=14,color='#20DFCC')
    # ax.annotate('RCP8.5: Ice free by mid-century',(0.365,0.031),(0.365,0.031),'axes fraction',fontproperties=FiraSansMed,fontsize=11,color='#DF8C20')
    # ax.annotate('SAI: Ice is maintained or even recovers',(0.235,0.7),(0.235,0.7),'axes fraction',fontproperties=FiraSansMed,fontsize=11,color='#20DFCC')
    ax.annotate('RCP8.5',(0.08,0.38),(0.08,0.38),'axes fraction',fontproperties=FiraSansMed,fontsize=12,color='#DF8C20')
    ax.annotate('Stratospheric Aerosol Injection',(0.53,0.75),(0.53,0.75),'axes fraction',fontproperties=FiraSansMed,fontsize=12,color='#20DFCC')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.ylim(-0.005,0.55)
    plt.xlim(setYear[0],setYear[1])

    savePrfx = ''
    savename = outDict["savePath"] + 'SeptArctic_IFRAC_20202095' + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

# Make images
for rlz in loopDict["realizations"]:
    setDict["realization"] = 'mean'
    glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}
    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        for reg in loopDict["regions"]:
            setDict["regOfInt"] = reg
            plot_zchfd_timeseries_sept_arc_ifrac(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
