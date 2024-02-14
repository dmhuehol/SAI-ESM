''' wrap_wrae_script
Plots warming rate vs. area exposed to given value of climate speed

Lots of the inputs in this script are NOT necessary to make a wrae plot, 
but are rather inherited from its construction from wrap_plot_slice_globe_script.
If I rewrote this from the ground up, I would make different decisions (and this
might happen someday)--but it works well enough for the current purposes!

dataDict: defines input data
setDict: settings for analysis/visualization
    plotScenarios: determines panel to plot, valid entries given below
        'RCP8.5': RCP8.5
        'SSP2-4.5': CESM2 SSP2-4.5 (from CESM2-ARISE)
        'UKESM-SSP2-4.5': UKESM SSP2-4.5 (from UKESM-ARISE)
        'GLENS-SAI': GLENS-SAI
        'ARISE-SAI-1.5': CESM2-ARISE-SAI-1.5
        'ARISE-SAI-DelayedStart': CESM2-ARISE-SAI-1.37-DelayedStart
        'ARISE-SAI-1.0': CESM2-ARISE-SAI-1.0
        'UKESM-ARISE-SAI-1.5': UKESM-ARISE-SAI-1.5
        'CESM2-WACCM:PreindustrialControl': CESM2 preindustrial control
outDict: output image settings
loopDict: determines which images are made

Fig. 3:
  Run script with default settings.
  
Supplementary Fig. 10:
  This plot requires some recoding within fun_calc_var, fun_derive_data, and fun_special_plot
  In fun_calc_var: See in-line comments. Uncomment desired line and comment others. (search on "Supplementary Fig. 8")
  In fun_derive_data: See in-line comments. Uncomment desired line and comment others. (search on "Supplementary Fig. 8")
  In fun_special_plot: Uncomment plt.xlim and plt.ylim settings (search on "Supplementary Fig. 8")

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from matplotlib import cm
import numpy as np
import cmasher, cmocean, seaborn  # Colormap packages

import fun_convert_unit as fcu
import fun_calc_var as fcv
import fun_process_data as fpd
import fun_special_plot as fsp
import get_periods as gp

per = gp.periods()
# Dictionaries to define inputs
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_2mTemp/',
    "idGlensCntrl": 'control_*',  # 'control_*' or None
    "idGlensFdbck": 'feedback_*',  # 'feedback_*' or None
    "idArise": '*DEFAULT*',  # '*DEFAULT*' or None
    "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    "idS245Hist": '*BWHIST*',  # '*BWHIST*' or None
    "idUkesmNoSai": '*ssp245*', #'*ssp245*' or None
    "idUkesmArise": '*arise-sai-1p5*', #'*arise-sai-1p5*' or None
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
    "idArise1p0": '*ARISE1P0*', # '*ARISE1P0*' or None
    "idPiControl": None, #'*piControl*'
    "idLastma": '*past1000*', #'*past1000*' or None
    "idS126": '*ssp126*', # '*ssp126' or None
    "idEra5": '*era5*', # '*era5*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_2mTemp/', # Need this info in both dictionaries for this script
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": ([2020, 2039],),
        "RCP8.5": ([2045, 2064],),
        "CESM2-ARISE": ([2035, 2054], ),
        "CESM2-SSP245": ([2045, 2064],),
        "CESM2-ARISE-DelayedStart": ([2045, 2064],),
        "UKESM-ARISE": ([2035, 2044],),
        "piControl": (
            [0, 19], [20, 39], [40, 59], [60, 79], [80, 99],
            [100, 119], [120, 139], [140, 159], [160, 179], [180, 199],
            [200, 219], [220, 239], [240, 259], [260, 279], [280, 299],
            [300, 319], [320, 339], [340, 359], [360, 379], [380, 399],
            [400, 419], [420, 439], [440, 459], [460, 479], [480, 499]),
        "LastMillennium": (per['per_all']),
        "Historical": ([1996, 2015],),
        "CESM2-SSP126": ([2045, 2064],),
        "ERA5": ([1996, 2015],),
        },
    "convert": 'See loop',  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-51,51],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'SSP2-4.5', 'Historical',
        'ARISE-SAI-1.5',
        'ARISE-SAI-DelayedStart', 'ARISE-SAI-1.0',
        'GLENS-SAI',
        'RCP8.5',
        'UKESM-ARISE-SAI-1.5',
        'UKESM-SSP2-4.5',
        # 'CESM2-WACCM:PreindustrialControl',
        'LastMillennium',
        'SSP1-2.6',
        'ERA5',
        ), # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20240213_periodsAndPaperFigs/',
    "dpiVal": 'pdf'
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    # TODO: Refactor so allplot is only behavior
    "levels": (None,),  # None for single-level variable
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at terminal

# Make images
wrCsList = list()
wrCsDictList = list()
for rlz in loopDict["rlzs"]:
    setDict["realization"] = rlz
    for dax in (
        (fcu.kel_to_cel, fcv.calc_warming_rate,), 
        (fcv.calc_climate_speed,)):
        setDict["convert"] = dax
        scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
        dataDict["landmaskFlag"] = setDict["landmaskFlag"] # Needed by wrae
        # You can put custom opening functions for datasets that can't
        # plug into call_to_open here. Remember this is "bespoke"--it 
        # needs to be replicable, not flexible.
        #
        # Note that plot_warmrate_areaexposed RELIES ON THE ORDER OF 
        # DATASETS IN THE INPUT LISTS. DO NOT MODIFY THIS.
        dataDict = {**dataDict, **cmnDict}
        wrCsList.append(scnList)
        wrCsDictList.append(dataDict)

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        fsp.plot_warmrate_areaexposed(
            wrCsList, wrCsDictList, setDict, outDict)

ic('Completed! :D')