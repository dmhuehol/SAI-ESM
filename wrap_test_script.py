''' wrap_plot_cvel_stats
Run plotting function for statistics of climate velocities,
e.g., range plots.

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
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*DEFAULT*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idDelayedStart": None, # '*DELAYED*' or None
    "idArise1p0": None, # '*ARISE1P0*' or None
    "idPiControl": '*piControl*', #'*piControl*'
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": 'See loop',  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": ([2020, 2039],),
        "CESM2-ARISE": ([2035, 2054], ),
        "CESM2-SSP245": ([2045, 2064],),
        "CESM2-ARISE-DelayedStart": ([2045, 2064],),
        "UKESM-ARISE": ([2035, 2054],),
        "piControl": (
            [10, 12], [48,50])
        # "piControl": (
        #     [10, 12], [48, 50], [100, 102],
        #     [129, 131], [169, 171], [264, 266],
        #     [285, 287], [341, 343], [384, 386], [471, 473])
        # "piControl": (
        #     [10, 29], [48, 67], [100, 119],
        #     [129, 148], [169, 188], [264, 283],
        #     [285, 304], [341, 360], [384, 403], [471, 490]),
        },
    "convert": (fcu.kel_to_cel, fcv.calc_climate_speed,),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmDisc,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-51,51],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'SSP2-4.5', 'ARISE-SAI-1.5', 'ARISE-SAI-DelayedStart', 
        'CESM2-WACCM:PreindustrialControl',
        ), # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230511_rangeArea/',
    "dpiVal": 'pdf'
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
    landoceanScnList = list()
    landoceanDdList = list()
    for landocean in ('land', 'ocean'):
        setDict["landmaskFlag"] = landocean
        setDict["realization"] = rlz
        scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
        dataDict = {**dataDict, **cmnDict}
        dataDict["landmaskFlag"] = setDict["landmaskFlag"] # Required info for rangeplot
        landoceanScnList.append(scnList)
        landoceanDdList.append(dataDict)
    
    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        fsp.plot_rangeplot(
            landoceanScnList, landoceanDdList, setDict, outDict)

ic('Completed! :D')