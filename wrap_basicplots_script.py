''' wrap_basicplots_script
Runs map plotting functions in fun_basic_plot.

dataDict: defines input data
setDict: settings for analysis/visualization
outDict: output image settings
loopDict: determines which images are made

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from matplotlib import cm
import cmasher, cmocean, seaborn  # Colormap packages

import fun_basic_plot as fbp
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Special color palettes
tropicalPal = seaborn.diverging_palette(133, 324, as_cmap=True)
precipPal = seaborn.diverging_palette(58, 162, s=100, l=45, as_cmap=True)

# Dictionaries to define inputs
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*',  # '*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    "idS245Hist": '*BWHIST*',  # '*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' # Landmask file location
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "startIntvl": [2015,2020,2030,2035],  # Window years [glens,glens,arise,arise]
    "endIntvl": [2025,2030,2040,2045],  # Window years [glens,glens,arise,arise]
    "convert": None,  # TUPLE of converter(s), None for default units
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-2,2],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": True  # True/False to run robustness
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20230201_revisitRefactor/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('ensplot',),  # number(s) for member(s), 'mean' ens mean all members, 'ensplot' for both member information and mean (i.e. for robustness)
    "levels": (None,),  # 'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global',),  # 'global' only implemented for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line

# Make images
for rlz in loopDict["rlzs"]:
    setDict["realization"] = rlz
    scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
    dataDict = {**dataDict, **cmnDict}

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        # fbp.plot_basic_difference_globe(scnList, dataDict, setDict, outDict)
        fbp.plot_single_basic_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_single_robust_globe(scnList, dataDict, setDict, outDict)

        # Will be removed for paper version of record
        # fbp.plot_six_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_basic_difference_polar(scnList, dataDict, setDict, outDict)
        # fbp.plot_glens_difference_globe(scnList, dataDict, setDict, outDict)
        # fbp.plot_arise_difference_globe(scnList, dataDict, setDict, outDict)
