''' wrap_plot_slice_globe_script
Run map plotting functions for a single time slice.

dataDict: defines input data
setDict: settings for analysis/visualization
    plotPanel: determines panel to plot, valid entries given below
        'RCP85': RCP8.5
        'CESMS245': CESM2 SSP2-4.5 (from CESM2-ARISE)
        'UKESMS245': UKESM SSP2-4.5 (from UKESM-ARISE)
        'GLENS': GLENS-SAI
        'ARISE15': CESM2-ARISE-SAI-1.5
        'UKESMARISE15': UKESM-ARISE-SAI-1.5
        any other value: make a blank map
outDict: output image settings
loopDict: determines which images are made

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from matplotlib import cm
import cmasher, cmocean, seaborn  # Colormap packages

import fun_plot_slice_globe as fpsg
import fun_convert_unit as fcu
import fun_calc_var as fcv
import fun_process_data as fpd
import region_library as rlib
#
# Special color palettes

# Dictionaries to define inputs
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_tas/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": '*arise-sai-1p5*', #'*arise-sai-1p5*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": [2015,2020],
        "CESM2-ARISE": [2035,2044],
        "UKESM-ARISE": [2035,2044]
        },
    "convert": (fcu.kel_to_cel, fcv.calc_decadal_climate_distance),  # TUPLE of converter(s) or calculators
    "cmap": cmasher.cm.dusk,  # None for default (cmocean balance) or choose colormap
    "cbVals": [0,10],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotPanel": 'UKESMARISE15' # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230210_climvel/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('ensplot',),  # number(s) for member(s), 'mean' ens mean all members, 'ensplot' for both member information and mean (i.e. for robustness)
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line

# Make images
for rlz in loopDict["rlzs"]:
    setDict["realization"] = rlz
    scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        fpsg.plot_single_slice_globe(
            scnList, dataDict, setDict, outDict)

ic('Completed! :D')
