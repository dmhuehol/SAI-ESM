# Empty script reserved for testing code
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
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_tas/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": '*ssp245*', #'*ssp245*' or None
    "idUkesmArise": '*arise-sai-1p5*', #'*arise-sai-1p5*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": 'land',  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "startIntvl": [2015,2020,2030,2035],  # Window years [glens,glens,arise,arise]
    "endIntvl": [2025,2030,2040,2045],  # Window years [glens,glens,arise,arise]
    "convert": None,  # TUPLE of converter(s), None for default units
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-2,2],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": True,  # True/False to run robustness
    "plotPanel": "intiGLENS" # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20230209_incorpUkesmClimvel/',
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
    ic(scnList, cmnDict)
    dataDict = {**dataDict, **cmnDict}

    # for lev in loopDict["levels"]:
    #     setDict["levOfInt"] = lev
    #     fbp.plot_single_basic_difference_globe(
    #         scnList, dataDict, setDict, outDict)
        # fbp.plot_single_robust_globe(
        #    scnList, dataDict, setDict, outDict)
