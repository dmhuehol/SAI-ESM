''' wrap_wrae_script
Plots warming rate vs. area exposed to given value (as for climate velocity)

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
    rlzs: 
        number or index for member(s), integer for CESM raw or 'r8i1p1f2' for CMIP
        'mean' ens mean all members
        'ensplot' for both member information and mean (i.e. for robustness)

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
import region_library as rlib
#
# Special color palettes
import fun_plot_tools as fpt
duskPink = np.array([241/256, 191/256, 202/256, 1])
tri_map = fpt.get_trifurcate_colormap('cmo.gray', 'cmr.amber', duskPink)
toWinter = [0, 0, 1, 1]
toSummer = [1, 0, 0, 1]
ssnPal = seaborn.diverging_palette(
    247, 321, s=100, l=50, as_cmap=True)
tri_div = fpt.get_trifurcate_div_colormap(toWinter, ssnPal, toSummer)
zmzmPal = seaborn.diverging_palette(
    285, 16, s=100, l=50, as_cmap=True)
    
zmzmDisc = fpt.get_cspd_colormap('zmzm')

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
    "idPiControl": '*piControl*', #'*piControl*'
    "idS126": '*ssp126*', # '*ssp126' or None
    "idEra5": '*era5*', # '*era5*' or None
    "idCruTs4": None, # '*cruts4*' or None
    # "idGlensCntrl": None,  # 'control_*' or None
    # "idGlensFdbck": None,  # 'feedback_*' or None
    # "idArise": None,  # '*DEFAULT*' or None
    # "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    # "idS245Hist": '*BWHIST*',  # '*BWHIST*' or None
    # "idUkesmNoSai": None, #'*ssp245*' or None
    # "idUkesmArise": None, #'*arise-sai-1p5*' or None
    # "idDelayedStart": None, # '*DELAYED*' or None
    # "idArise1p0": None, # '*ARISE1P0*' or None
    # "idPiControl": None, #'*piControl*'
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": ([2020, 2039],),
        "RCP8.5": ([2045, 2064],),
        "CESM2-ARISE": ([2035, 2054], ),
        "CESM2-SSP245": ([2045, 2064],),
        "CESM2-ARISE-DelayedStart": ([2045, 2064],),
        "UKESM-ARISE": ([2035, 2044],),
        # "piControl": (
        #     [10, 12], [48,50])
        # "piControl": (
        #     [10, 19], [48, 57], [100, 109],
        #     [129, 138], [169, 178], [264, 273],
        #     [285, 294], [341, 350], [384, 393], [471, 480])
        # "piControl": (
        #     [10, 29], [48, 67], [100, 119],
        #     [129, 148], [169, 188], [264, 283],
        #     [285, 304], [341, 360], [384, 403], [471, 490]),
        "piControl": (
            [0, 19], [20, 39], [40, 59], [60, 79], [80, 99],
            [100, 119], [120, 139], [140, 159], [160, 179], [180, 199],
            [200, 219], [220, 239], [240, 259], [260, 279], [280, 299],
            [300, 319], [320, 339], [340, 359], [360, 379], [380, 399],
            [400, 419], [420, 439], [440, 459], [460, 479], [480, 499]),
        "Historical": ([1996, 2015],),
        "CESM2-SSP126": ([2045, 2064],),
        "ERA5": ([1996, 2015],),
        "CRU_TS4": ([1996, 2015],),
        },
    "convert": 'See loop',  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmDisc,  # None for default (cmocean balance) or choose colormap
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
        'CESM2-WACCM:PreindustrialControl',
        'SSP1-2.6',
        # 'ERA5',
        ), # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230824_forScri/',
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
        # needs to be replicable, not flexible. (Kind of like your old
        # ice diagram code!)
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