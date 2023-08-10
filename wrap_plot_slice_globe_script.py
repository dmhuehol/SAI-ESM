''' wrap_plot_slice_globe_script
Run map plotting functions for a single time slice.

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
import region_library as rlib

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
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
    "idArise1p0": None, # '*ARISE1P0*' or None
    "idPiControl": None, #'*piControl*'
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": 'land',  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "calcIntvl": { # Years to calculate
        "GLENS": ([2020, 2039],),
        "RCP8.5": ([2045, 2064],),
        "CESM2-ARISE": ([2035, 2054], ),
        "CESM2-SSP245": ([2045, 2064],),
        "CESM2-ARISE-DelayedStart": ([2045, 2064],),
        "UKESM-ARISE": ([2035, 2044],),
        "piControl": ( #main
            [5, 24], [43, 62], [95, 114],
            [124, 143], [164, 183], [259, 278],
            [280, 299], [336, 355], [379, 398], [465, 484]),
        # "piControl": ( #alternate
        #     [17, 36], [93, 112], [132, 151],
        #     [234, 253], [300, 319], [321, 340],
        #     [367, 386], [392, 411], [412, 431], [438, 457]),
        },
    "convert": (fcu.kel_to_cel, fcv.calc_climate_speed_forced_response,),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmDisc,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-5051, 5051],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'ARISE-SAI-DelayedStart',
        # 'ARISE-SAI-1.5',
        # 'ARISE-SAI-1.0',
        # 'UKESM-ARISE-SAI-1.5',
        # 'CESM2-WACCM:PreindustrialControl',
        # 'SSP2-4.5',
        # 'GLENS-SAI',
        # 'RCP8.5',
        ), # See docstring for valid inputs
    "plotEnsType": '' # 'mean', 'max'/'min' pointwise max/min, number for single member
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230803_orderAndCode/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line
# ic(np.diff(setDict["calcIntvl"]["piControl"])) # Check to ensure correct time periods

# Make images
for pet in ('mean',):
    setDict["plotEnsType"] = pet
    for rlz in loopDict["rlzs"]:
        setDict["realization"] = rlz
        scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
        dataDict = {**dataDict, **cmnDict}
    
        for lev in loopDict["levels"]:
            setDict["levOfInt"] = lev
            fpsg.plot_single_slice_globe(
                scnList, dataDict, setDict, outDict)
            # fpsg.plot_single_slice_vector_globe(
            #     scnList, dataDict, setDict, outDict)

ic('Completed! :D')
