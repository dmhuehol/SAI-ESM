''' wrap_basicplots_script
Runs plotting functions in fun_basic_plot.

dataDict is for inputs
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
import cmasher
import seaborn
import numpy as np

import fun_basic_plot as fbp
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
testAllTypes = rlib.atlas_all_types() #testAllTypes["allRegions"]
gnsht = ('global', rlib.Arctic(), rlib.HudsonBay(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere(),)

# Specials
tropicalPal = seaborn.diverging_palette(133, 324, as_cmap=True)

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/eimnsn_PRECT/anmn/',
    "idGlensCntrl": 'control_*', #'control_*' or None
    "idGlensFdbck": 'feedback_*', #'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*', #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', #'*BWSSP245*' or None
    "idS245Hist": '*BWHIST*', #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc'
}
setDict = {
    "landmaskFlag": None,
    "startIntvl": [2015,2020,2030,2035], #dg [2015,2020,2030,2035]
    "endIntvl": [2025,2030,2040,2045], #dg [2025,2030,2040,2045]
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "arisePoi": [2041], #pdf
    "s245CntrlPoi": [2011, 2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "convert": (fcu.m_to_cm, fcu.persec_peryr), #TUPLE of converter(s), or None if using default units
    "cmap": cmocean.cm.curl_r, #None for default cmocean "balance" or choose colormap here
    "cbVals": [-50, 50] #None for automatic or [min,max] to override #dg
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20220209_polishMergeFinal/1_polish/',
    "dpiVal": 400
}
loopDict = {
    "realizations": ('mean',), #number for individual member, 'mean' for ens mean of all available members
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global',),#('global',rlib.Arctic(),rlib.EastNorthAmerica()),
    "aaBools": (True,)
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

# Make images
for rlz in loopDict["realizations"]:
    setDict["realization"] = rlz
    scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        # fbp.plot_basic_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_six_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_single_basic_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_basic_difference_polar(scnList, dataDict, setDict, outDict)
        # fbp.plot_glens_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_arise_difference_globe(scnList, dataDict, setDict, outDict)
