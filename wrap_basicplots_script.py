''' wrap_basicplots_script
Runs plotting functions in basic_plot_fun. This is used by run_plots_script.sh
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

import basic_plot_fun as bpf
import ens_plot_fun as epf
import fun_convert_unit as fcu
import process_glens_fun as pgf
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
gnsht = ('global', rlib.Arctic(), rlib.SouthernHemisphere(),)

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/sept_IFRAC/AllVersions/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*',
    "fnameGlens2": '*SSP245*'#'*SSP245*',
}
setDict = {
    "toPlot": ('GLENS1:Control', 'GLENS1:Feedback', 'GLENS2:Feedback'),
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "glens2Poi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "convert": (fcu.fraction_to_percent,) #TUPLE of converter(s), or None if using default units
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210730_addGlens2/',
    "dpiVal": 400
}
loopDict = {
    "realizations": ('mean',), #number for individual member, 'mean' for ens mean of all available members, or 'ensplot' for ensemble spaghetti
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global',rlib.Arctic(),),
    "aaBools": (True,)
}

# Verify inputs (troubleshooting)
ic(setDict["convert"])

# Make images
for rlz in loopDict["realizations"]:
    setDict["realization"] = rlz
    rlzList, cmnDict = pgf.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}
    # bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict) #comment if running variable with no levels
    # bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict) #comment out if running variable with no levels

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        # bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
        # bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

        for reg in loopDict["regions"]:
            setDict["regOfInt"] = reg
            # if rlz == 'ensplot':
            #     epf.plot_ens_spaghetti_timeseries(rlzList, dataDict, setDict, outDict)
            # else:
            #     bpf.plot_timeseries(rlzList, dataDict, setDict, outDict)

            for aab in loopDict["aaBools"]:
                setDict["areaAvgBool"] = aab
                setDict["plotStyle"] = 'step'
                bpf.plot_pdf(rlzList, dataDict, setDict, outDict)
                setDict["plotStyle"] = 'kde'
                bpf.plot_pdf(rlzList, dataDict, setDict, outDict)
                setDict["plotStyle"] = 'hist'
                bpf.plot_pdf(rlzList, dataDict, setDict, outDict)
