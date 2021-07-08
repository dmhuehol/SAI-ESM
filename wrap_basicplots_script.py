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
import fun_convert_unit as fcu
import process_glens_fun as pgf
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
gnsh = ('global', rlib.NorthernHemisphere(), rlib.SouthernHemisphere())

# Dictionaries
dataDict = {
    "dataPath": '/glade/scratch/dhueholt/annual_T/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*'
}
setDict = {
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "convert": fcu.kel_to_cel, #relevant unit converter or None if using default units
}
outDict = {
    "savePath": '/glade/work/dhueholt/20210708_plots/',
    "dpiVal": 400
}
loopDict = {
    "realizations": (1,2,3,21,'mean',), #number for individual member or 'mean' for ens mean of all available members
    "levels": (1000,), #'stratosphere', 'troposphere', 'total', or numeric level(s)
    "regions": gnsh,
    "aaBools": (False,)
}

# Verify inputs (troubleshooting)
ic(setDict["convert"])

# Make images
for rlz in loopDict["realizations"]:
    setDict["realization"] = rlz
    glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}
    bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
        bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

        # for reg in loopDict["regions"]:
            # setDict["regOfInt"] = reg
            # bpf.plot_timeseries(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

            # for aab in loopDict["aaBools"]:
                # setDict["areaAvgBool"] = aab
                # setDict["plotStyle"] = 'step'
                # bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
                # setDict["plotStyle"] = 'kde'
                # bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
                # setDict["plotStyle"] = 'hist'
                # bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
