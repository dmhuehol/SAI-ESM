''' wrap_ensplots_script
Runs ensemble plotting functions in fun_ens_plot.

dataDict is for data inputs
setDict sets settings related to plotting
outDict is for outputs
loopDict determines which images are made

Written by Daniel Hueholt | February 2022
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean
import numpy as np

import fun_ens_plot as fep
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
seaIcyRegions = rlib.atlas_seaicy_regions()
planetary = ('global', rlib.Tropics(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere())
insets = rlib.atlas_insets()

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/extreme_MHW/roll5_centerFalseShift4/',
    "idGlensCntrl": 'control_*', #'control_*' or None
    "idGlensFdbck": 'feedback_*', #'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*', #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', #'*BWSSP245*' or None
    "idS245Hist": None, #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": None, #None or 'land'
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "arisePoi": [2041], #pdf
    "s245CntrlPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "areaAvgBool": True, #True=mean, 'sum'=sum (e.g. ice extent), False=functionality that hasn't been tested in two years, and yes I know that's not a boolean
    "dimOfVrblty": {'rlzBool':True,'timeBool':True,'spcBool':False}, #pdf
    "convert": (fcu.roll5_to_percentdays,), #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "insetFlag": 2, #ts 0 default, 1 lines only, 2 AMS style
    "mute": False, #ts True/False to use image muting on parts of time period
    "ylim": [0,100], #ts None for automatic, [min,max] for manual
    "ylabel": '', #ts None for automatic
    "yticks": np.arange(0,101,30), #ts None for automatic, np.arange(mn,mx,step) for manual
    "xticks": True
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20220401_mhwLastSteps/9_aesthetic/',
    "dpiVal": 400
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": (rlib.WesternAustraliaMHW_point(), )
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

# Make images
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
for lev in loopDict["levels"]:
    setDict["levOfInt"] = lev
    for reg in loopDict["regions"]:
        setDict["regOfInt"] = reg
        # fep.plot_ens_spaghetti_timeseries(scnList, dataDict, setDict, outDict)
        fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'step'
        # fep.plot_ens_pdf(scnList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'kde'
        # fep.plot_ens_pdf(scnList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'hist'
        # fep.plot_ens_pdf(scnList, dataDict, setDict, outDict)
