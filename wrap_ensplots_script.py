''' wrap_ensplots_script
Runs ensemble plotting functions in fun_ens_plot.

dataDict: defines input data
setDict: settings for analysis/visualization
    areaAvgBool: sets area averaging behavior, more accurately named areaAvgFlag
    styleFlag: choose preset visual style for ens_spread_timeseries figures
outDict: output image settings
loopDict: determines which images are made

Supplementary Fig. 1:
  Run with default settings

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

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Desktop/OSF/CSHF/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*DEFAULT*',  # '*DEFAULT*' or None
    "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
    "idArise1p0": None, # '*ARISE1P0*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": None, # None or 'land'
    "areaAvgBool": True, # see docstring for valid inputs
    "convert": (fcu.kel_to_cel,), #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "styleFlag": 5, # see docstring for valid inputs
    "mute": False, #True/False to use image muting on parts of time period
    "ylim": [14.8,16.7], #None for automatic, [min,max] for manual
    "ylabel": '\u00b0C', #None for automatic
    "yticks": None, #None for automatic, np.arange(mn,mx,step) for manual
    "xticks": True, #True/False to enable/disable x-axis tick labels
}
outDict = {
    "savePath": '/Users/dhueholt/Desktop/OSF/images/',
    "addToSaveStr": '',
    "dpiVal": 'pdf'
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global', ) #rlib.region() for single (rlib.region(), rlib.region()) to concatenate)
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
