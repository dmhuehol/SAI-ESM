''' wrap_rangeplot_script
Run plotting function for the range of variability of climate speeds.
Lots of the inputs in this script are NOT necessary to make a range plot, 
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
    rlzs: 
        number or index for member(s), integer for CESM raw or 'r8i1p1f2' for CMIP
        'mean' ens mean all members
        'ensplot' for both member information and mean (i.e. for robustness)

Fig. 2a:
    setDict["magBool"]: True

Supplementary Fig. 7
    setDict["magBool"]: False

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

# Dictionaries to define inputs
dataDict = {
    "dataPath": '/Users/dhueholt/Desktop/OSF/CSHF/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*DEFAULT*',  # '*DEFAULT*' or None
    "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
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
            [10, 29], [48, 67], [100, 119],
            [129, 148], [169, 188], [264, 283],
            [285, 304], [341, 360], [384, 403], [471, 490]),
        },
    "convert": (fcu.kel_to_cel, fcv.calc_climate_speed,),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-51,51],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'SSP2-4.5', 'ARISE-SAI-1.5', 'ARISE-SAI-DelayedStart', 
        'CESM2-WACCM:PreindustrialControl',
        ), # See docstring for valid inputs
    "magBool": True #Plot magnitude true/false
}
outDict = {
    "savePath": '/Users/dhueholt/Desktop/OSF/images/',
    "dpiVal": 'pdf'
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    "levels": (None,),  # None for single-level variable
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