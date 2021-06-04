from icecream import ic

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean

import plotting_tools as plt_tls
import plot_difference_globes

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_o3/',
    "fnameCntrl": 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc',
    "fnameFdbck": 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
}

setDict = {
    "startIntvl": [2010,2019],
    "endIntvl": [2090,2099],
    "levOfInt": 200, #'stratosphere', 'troposphere', 'total', numeric level, or list of numeric levels
    "quantileOfInt": 0.67
}

outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20210604_refactoringZonalAndMore/',
    "dpiVal": 400
}


plot_difference_globes.plot_basic_difference_globe(dataDict, setDict, outDict)
plot_difference_globes.plot_single_basic_difference_globe(dataDict, setDict, outDict)
plot_difference_globes.plot_vertical_difference_globe(dataDict, setDict, outDict)
plot_difference_globes.plot_vertical_baseline_difference_globe(dataDict, setDict, outDict)

ic('Completed! :D')
