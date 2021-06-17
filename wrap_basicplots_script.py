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

import plotting_tools as plt_tls
import basic_plot_fun as bpf
import process_glens_fun as pgf
import region_library as rlib

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_U/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*'
}
setDict = {
    "realization": 3, #number for individual member or 'mean' for ensemble mean | TODO: array entry to choose particular members
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2076], #pdf
    "fdbckPoi": [2076], #pdf
    "timePeriod": 20, #pdf
    "levOfInt": 50, #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
    "regOfInt": rlib.NoLandLatitude(), #ts, pdf
    "areaAvgBool": False, #pdf
    "plotStyle": 'step', #pdf
    "quantileOfInt": 0.67 #dg
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210616_finalPrep/2_fnamesCompleteCheck/',
    "dpiVal": 400
}

glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}

bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

bpf.plot_timeseries(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)

bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
setDict["plotStyle"] = 'kde'
bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
setDict["plotStyle"] = 'hist'
bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
