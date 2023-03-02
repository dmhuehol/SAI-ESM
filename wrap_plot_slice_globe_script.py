''' wrap_plot_slice_globe_script
Run map plotting functions for a single time slice.

dataDict: defines input data
setDict: settings for analysis/visualization
    plotPanel: determines panel to plot, valid entries given below
        'RCP85': RCP8.5
        'CESMS245': CESM2 SSP2-4.5 (from CESM2-ARISE)
        'UKESMS245': UKESM SSP2-4.5 (from UKESM-ARISE)
        'GLENS': GLENS-SAI
        'CESMARISE15': CESM2-ARISE-SAI-1.5
        'UKESMARISE15': UKESM-ARISE-SAI-1.5
        'CESM2PIC': CESM2 preindustrial control
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

import fun_plot_slice_globe as fpsg
import fun_convert_unit as fcu
import fun_calc_var as fcv
import fun_process_data as fpd
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
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_tas/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*',  # '*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": '*arise-sai-1p5*', #'*arise-sai-1p5*' or None
    "idPiControl": '*piControl*', #'*piControl*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landtropicmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": 'land',  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": [2020, 2029],
        "CESM2-ARISE": [2035, 2044],
        "UKESM-ARISE": [2035, 2044],
        "piControl": [10, 19],
        },
    "convert": (fcu.kel_to_cel, fcv.calc_decadal_climate_distance,),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmDisc,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-51,51],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotPanel": 'UKESMARISE15', # See docstring for valid inputs
    "plotEnsType": 2 #'mean', 'max'/'min' pointwise max/min, number for single member
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230302_pi/',
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
    dataDict = {**dataDict, **cmnDict}

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        fpsg.plot_single_slice_globe(
            scnList, dataDict, setDict, outDict)

ic('Completed! :D')
