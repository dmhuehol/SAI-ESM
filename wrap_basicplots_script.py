''' wrap_basicplots_script
Runs plotting functions in basic_plot_fun. This is used by run_plots_script.sh
to submit jobs through the NCAR queue.

Written by Daniel Hueholt | June 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean
import numpy as np

import plotting_tools as plt_tls
import basic_plot_fun as bpf
import process_glens_fun as pgf
import region_library as rlib

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_Q/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*'
}
setDict = {
    "realization": 'mean', #number for individual member or 'mean' for ensemble mean | TODO: array entry to choose particular members
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "levOfInt": 1000, #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
    "regOfInt": rlib.AlaskaNorthwestCanada(), #ts, pdf
    "areaAvgBool": False, #pdf
    "plotStyle": 'step', #pdf
    "quantileOfInt": 0.67 #dg
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210629_refEnsAndNewPlots/2_fullCheck/',
    "dpiVal": 400
}

# Batch using loops
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5()
loopDict = {
    "realizations": (3,4,'mean',),#(1,2,3,21,'mean'),
    "levels": (1000,),
    "regions": (rlib.Amazon(),),
    "aaBools": (True,False)
}

for rzc in loopDict["realizations"]:
    setDict["realization"] = rzc
    glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}
    bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        for reg in loopDict["regions"]:
            setDict["regOfInt"] = reg
            glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
            dataDict = {**dataDict, **cmnDict}
            bpf.plot_timeseries(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
            for aab in loopDict["aaBools"]:
                setDict["areaAvgBool"] = aab
                setDict["plotStyle"] = 'step'
                glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
                dataDict = {**dataDict, **cmnDict}

                bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
                setDict["plotStyle"] = 'kde'
                bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
                setDict["plotStyle"] = 'hist'
                bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

# One at a time
# glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}

# bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
# bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
# bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
# bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

# bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
# setDict["plotStyle"] = 'kde'
# bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
# setDict["plotStyle"] = 'hist'
# bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

# ic('Plots completed!')
