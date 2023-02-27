# Try calculating seasonal shifts
from icecream import ic
import sys

from matplotlib import cm
import numpy as np
import cmasher, cmocean, seaborn  # Colormap packages

import fun_plot_slice_globe as fpsg
import fun_convert_unit as fcu
import fun_calc_var as fcv
import fun_process_data as fpd
import region_library as rlib

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/monthly_tas/mergetime/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": '*ssp245*', #'*ssp245*' or None
    "idUkesmArise": '*arise-sai-1p5*', #'*arise-sai-1p5*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": [2020, 2029],
        "CESM2-ARISE": [2045, 2054],
        "UKESM-ARISE": [2035, 2044]
        },
    "convert": (fcu.kel_to_cel, fcv.calc_seasonal_shift),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-51,51],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotPanel": 'UKESMS245', # See docstring for valid inputs
    "plotEnsType": 4 #'mean', 'max'/'min' pointwise max/min, number for single member
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230224_seasonShift/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    # TODO: Refactor so allplot is only behavior
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line

# Make images
for rlz in loopDict["rlzs"]:
    setDict["realization"] = rlz
    scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
    ic(scnList, cmnDict)
    sys.exit('STOP')
    dataDict = {**dataDict, **cmnDict}


ic('Hooray!')