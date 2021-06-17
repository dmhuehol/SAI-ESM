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
    "dataPath": '/glade/scratch/dhueholt/annual_T/',
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
    "regOfInt": 'global', #ts, pdf
    "areaAvgBool": False, #pdf
    "plotStyle": 'step', #pdf
    "quantileOfInt": 0.67 #dg
}
outDict = {
    "savePath": '/glade/work/dhueholt/20210617_begin/',
    "dpiVal": 400
}

# Batch using loops
for rzc in (1,2,3,21,'mean'):
    setDict["realization"] = rzc
    # glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
    # dataDict = {**dataDict, **cmnDict}
    # bpf.plot_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    # bpf.plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    # bpf.plot_vertical_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    # bpf.plot_vertical_baseline_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
    for region in ('global'):
        setDict["regOfInt"] = region
        setDict["plotStyle"] = 'step'
        glensCntrlRlz, glensFdbckRlz, cmnDict = pgf.call_to_open(dataDict, setDict)
        dataDict = {**dataDict, **cmnDict}
        try:
            bpf.plot_timeseries(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
            bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
            setDict["plotStyle"] = 'kde'
            bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
            setDict["plotStyle"] = 'hist'
            bpf.plot_pdf(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict)
        except:
            ic('Failed on: ' + str(rzc))
            continue

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
