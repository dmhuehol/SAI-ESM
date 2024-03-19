''' wrap_plot_slice_globe_script
Run vector map plotting functions for a single time slice. This could
be handled as part of wrap_plot_slice_globe_script, but requires enough
retuning that it's easier to write separately.

dataDict: defines input data
setDict: settings for analysis/visualization
    plotPanel: determines panel to plot, valid entries given below
        'RCP8.5': RCP8.5
        'SSP2-4.5': CESM2 SSP2-4.5 (from CESM2-ARISE)
        'UKESM-SSP2-4.5': UKESM SSP2-4.5 (from UKESM-ARISE)
        'GLENS-SAI': GLENS-SAI
        'ARISE-SAI-1.5': CESM2-ARISE-SAI-1.5
        'ARISE-SAI-DelayedStart': CESM2-ARISE-SAI-1.37-DelayedStart
        'ARISE-SAI-1.0': CESM2-ARISE-SAI-1.0
        'UKESM-ARISE-SAI-1.5': UKESM-ARISE-SAI-1.5
        'LastMillennium': CESM2(WACCM6ma) Last Millennium
        'CESM2-WACCM:PreindustrialControl': CESM2 preindustrial control
outDict: output image settings
loopDict: determines which images are made
    rlzs: 
        number or index for member(s), integer for CESM raw or 'r8i1p1f2' for CMIP
        'mean' ens mean all members
        'ensplot' for both member information and mean (i.e. for robustness)

Dictionaries to define inputs
  
To plot vector maps (Supplementary Fig. 11abc):
Supplementary Fig. 11a:
  dataDict["idS245Cntrl"]: '*BWSSP245*'
  All other "id" keys: None
  setDict["landmaskFlag"]: None
  setDict["plotScenarios"]: ('SSP2-4.5',)
  
Supplementary Fig. 11b:
  dataDict["idArise"]: '*DEFAULT*'
  All other "id" keys: None
  setDict["landmaskFlag"]: None
  setDict["plotScenarios"]: ('ARISE-SAI-1.5',)
  
Supplementary Fig. 11c:
  dataDict["idArise"]: '*DELAYED*'
  All other "id" keys: None
  setDict["landmaskFlag"]: None
  setDict["plotScenarios"]: ('ARISE-SAI-DelayedStart',)
  

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

import cmasher, cmocean, seaborn  # Colormap packages
from icecream import ic
from matplotlib import cm
import numpy as np

import fun_convert_unit as fcu
import fun_calc_var as fcv
import fun_plot_slice_globe as fpsg
import fun_plot_tools as fpt
import fun_process_data as fpd
import get_periods as gp

zmzmVec = fpt.get_cspd_vector_colormap('zmzm')
per = gp.periods()

dataDict = {
    "dataPath": '/Users/dhueholt/Desktop/OSF/ClimateSpeeds/',
    "idArise": None,  # '*DEFAULT*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
    "mask": '/Users/dhueholt/Desktop/OSF/ClimateSpeeds/Masks/cesm_atm_mask.nc', # Landmask file location (CESM)
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "CESM2-ARISE": ([2035, 2054], ),
        "CESM2-SSP245": ([2045, 2064],),
        "CESM2-ARISE-DelayedStart": ([2045, 2064],),
        },
    "convert": (fcu.kel_to_cel, fcv.calc_climate_velocity),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmVec,  # None for default (cmocean balance) or choose colormap
    "cbVals": None,  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'ARISE-SAI-DelayedStart',
        # 'ARISE-SAI-1.5',
        # 'SSP2-4.5',
        ), # See docstring for valid inputs
    "plotEnsType": '' # 'mean', 'max'/'min' pointwise max/min, number for single member
}
outDict = {
    "savePath": 'ENTER PATH',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line

# Make images
for pet in ('mean',): # For ensemble mean figures
    setDict["plotEnsType"] = pet
    for rlz in loopDict["rlzs"]:
        setDict["realization"] = rlz
        scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
        dataDict = {**dataDict, **cmnDict}
    
        for lev in loopDict["levels"]:
            setDict["levOfInt"] = lev
            fpsg.plot_single_slice_vector_globe(
                scnList, dataDict, setDict, outDict)

ic('Completed! :D')
