''' wrap_basicplots_script
Runs plotting functions in basic_plot_fun. This is used by run_plots_script.sh
to submit jobs through the NCAR queue.

Written by Daniel Hueholt | June 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean

import plotting_tools as plt_tls
import basic_plot_fun as bpf
import region_library as rlib

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_U/',
    "fnameCntrl": 'control_*',
    "fnameFdbck": 'feedback_*'
}
# control_*
# feedback_*
# control_003_U_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc
# feedback_003_U_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc

setDict = {
    "realization": 'mean', #number for individual member or 'mean' for ensemble mean
    "startIntvl": [2010,2019],
    "endIntvl": [2090,2099],
    "cntrlPoi": [2011,2076],
    "fdbckPoi": [2076],
    "timePeriod": 20,
    "levOfInt": 50, #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
    "regOfInt": 'global',
    "areaAvgBool": False,
    "plotStyle": 'step',
    "quantileOfInt": 0.67
}

outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210615_polishingEns/',
    "dpiVal": 400
}

# bpf.plot_basic_difference_globe(dataDict, setDict, outDict)
# bpf.plot_single_basic_difference_globe(dataDict, setDict, outDict)
# bpf.plot_vertical_difference_globe(dataDict, setDict, outDict)
# bpf.plot_vertical_baseline_difference_globe(dataDict, setDict, outDict)
#
# bpf.plot_timeseries(dataDict, setDict, outDict)

bpf.plot_pdf(dataDict, setDict, outDict)
# setDict["plotStyle"] = 'kde'
# bpf.plot_pdf(dataDict, setDict, outDict)
# setDict["plotStyle"] = 'hist'
# bpf.plot_pdf(dataDict, setDict, outDict)
