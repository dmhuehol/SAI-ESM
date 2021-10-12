''' wrap_basicplots_script
Runs plotting functions in fun_basic_plot. This is used by run_plots_script.sh
to submit jobs through the NCAR queue.

dataDict is for inputs
setDict sets settings related to plotting
outDict is for outputs
loopDict determines which images are made

Written by Daniel Hueholt | July 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean
import numpy as np

import fun_basic_plot as fbp
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
testAllTypes = rlib.atlas_all_types() #testAllTypes["allRegions"]
gnsht = ('global', rlib.Arctic(), rlib.HudsonBay(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere(),)

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/feb_hi/',
    "idGlensCntrl": 'control_*', #'control_*'
    "idGlensFdbck": 'feedback_*', #'feedback_*'
    "idSciris": '*SSP245-TSMLT-GAUSS*', #'*SSP245*'
    "idS245Cntrl": '*BWSSP245*', #'*ssp245*'
    "idS245Hist": '*BWHIST*', #'*historical*'
    "idMask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc'
}
setDict = {
    "landmaskFlag": None,
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "scirisPoi": [2041], #pdf
    "s245CntrlPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "convert": None #TUPLE of converter(s), or None if using default units
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20211012_syncWithCasperPlusExtremes/1_testBeforeSync/',
    "dpiVal": 400
}
loopDict = {
    "realizations": ('mean',), #number for individual member, 'mean' for ens mean of all available members
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": (rlib.Antarctica(),),#('global',rlib.Arctic(),rlib.EastNorthAmerica()),
    "aaBools": (True,)
}

# Verify inputs (troubleshooting)
ic(dataDict)
ic(setDict)
ic(outDict)
# ic(setDict["convert"])

# Make images
for rlz in loopDict["realizations"]:
    setDict["realization"] = rlz
    rlzList, cmnDict = fpd.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}
    # fbp.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict) #comment if running variable with no levels
    # fbp.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict) #comment out if running variable with no levels

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        # fbp.plot_basic_difference_globe(rlzList, dataDict, setDict, outDict)
        fbp.plot_basic_difference_polar(rlzList, dataDict, setDict, outDict)
        # fbp.plot_single_basic_difference_globe(rlzList, dataDict, setDict, outDict)
    #
    #     for reg in loopDict["regions"]:
    #         setDict["regOfInt"] = reg
    #         fbp.plot_timeseries(rlzList, dataDict, setDict, outDict)
    #
    #         for aab in loopDict["aaBools"]:
    #             setDict["areaAvgBool"] = aab
    #             setDict["plotStyle"] = 'step'
    #             fbp.plot_pdf(rlzList, dataDict, setDict, outDict)
    #             setDict["plotStyle"] = 'kde'
    #             fbp.plot_pdf(rlzList, dataDict, setDict, outDict)
    #             setDict["plotStyle"] = 'hist'
    #             fbp.plot_pdf(rlzList, dataDict, setDict, outDict)
