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
        'LastMillennium': CESM2(WACCM6ma) Last Millennium
        'CESM2-WACCM:PreindustrialControl': CESM2 preindustrial control
outDict: output image settings
loopDict: determines which images are made
    rlzs: 
        number or index for member(s), integer for CESM raw or 'r8i1p1f2' for CMIP
        'mean' ens mean all members
        'ensplot' for both member information and mean (i.e. for robustness)

Dictionaries to define inputs
Fig. 1a: 
  dataDict["idS245Cntrl"]: '*BWSSP245*'
  All other "id" keys: None
  setDict["landmaskFlag"]: 'land'
  setDict["plotScenarios"]: ('SSP2-4.5',)
Fig. 1b:
  Same as Fig. 1a, but:
  setDict["landmaskFlag"]: 'ocean'
Fig. 1c:
  dataDict["idLastma"]: '*past1000*'
  All other "id" keys: None
  setDict["landmaskFlag"]: 'land'
  setDict["plotScenarios"]: ('LastMillennium',)
Fig. 1d:
  Same as Fig. 1c, but:
  setDict["landmaskFlag"]: 'ocean'
Fig. 1e:
  dataDict["idArise"]: '*DEFAULT*'
  All other "id" keys: None
  setDict["landmaskFlag"]: 'land'
  setDict["plotScenarios"]: ('ARISE-SAI-1.5',)
Fig. 1f:
  Same as Fig. 1e, but:
  setDict["landmaskFlag"]: 'ocean'
Fig. 1g:
  dataDict["idDelayedStart"]: '*DELAYED*'
  All other "id" keys: None
  setDict["landmaskFlag"]: 'land'
  setDict["plotScenarios"]: ('ARISE-SAI-DelayedStart',)
Fig. 1h:
  Same as Fig. 1g, but:
  setDict["landmaskFlag"]: 'ocean'

To plot ARISE-DelayedStart - ARISE-1.5 difference (Supplementary Fig. 2ab)
  Supplementary Fig. 2a:
    dataDict["idArise"]: '*DEFAULT*'
    dataDict["idDelayedStart"]: '*DELAYED*'
    All other "id" keys: None
    setDict["landmaskFlag"]: 'land'
    setDict["plotScenarios"]: ('ARISE-SAI-DelayedStart', 'ARISE-SAI-1.5', 'DS-15')
  Supplementary Fig. 2b:
    Same as 2a, except:
    setDict["landmaskFlag"]: 'ocean'

To plot Unforced (Supplementary Fig. 3ab):
  Supplementary Fig. 3a:
    dataDict["idPiControl"]: '*piControl*'
    All other ["id"] keys: None
    setDict["landmaskFlag"]: 'land'
    setDict["plotScenarios"]: 'CESM2-WACCM:PreindustrialControl'
  Supplementary Fig. 3b:
    Same as 3a, except:
    setDict["landmaskFlag"]: 'ocean'
  
To plot all ensemble members (Supplementary Fig. 4-7): 
  Uncomment "for pet in (0,1,2,3,4,5,6,7,8,9):"
  Comment "for pet in ('mean',)"
Supplementary Fig. 4:
  Top half: Same as Fig. 1a
  Bottom half: Same as Fig. 1b
Supplementary Fig. 5: 
  Top half: Same as Fig. 1c
  Bottom half: Same as Fig. 1d
Supplementary Fig. 6:
  Top half: Same as Fig. 1e
  Bottom half: Same as Fig. 1f
Supplementary Fig. 7: 
  Top half: Same as Fig. 1g
  Bottom half: Same as Fig. 1h
Supplementary Fig. 8
  Top half: Same as Supplementary Fig. 3a
  Top half: Same as Supplementary Fig. 3b
Fig. 2b, 2c, 2d correspond to members 4, 2, and 6 of ARISE-DelayedStart (Supplementary Fig. 6)
  Note that because of Python's zero indexing, this corresponds to indices
  3, 1, and 5.

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

zmzmDisc = fpt.get_cspd_colormap('zmzm')
per = gp.periods()

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
    "idPiControl": None, #'*piControl*' or None
    "idLastma": '*past1000*', #'*past1000*' or None
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
        "piControl": (
            [5, 24], [43, 62], [95, 114],
            [124, 143], [164, 183], [259, 278],
            [280, 299], [336, 355], [379, 398], [465, 484]),
        "LastMillennium": ([1646, 1665],)#(per['per_ens']),
        },
    "convert": (fcu.kel_to_cel, fcv.calc_climate_speed),  # TUPLE of converter(s) or calculators from fun_convert_unit or fun_calc_var
    "cmap": zmzmDisc,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-5051, 5051],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotScenarios": (
        'LastMillennium',
        # 'ARISE-SAI-DelayedStart',
        # 'ARISE-SAI-1.5',
        # 'DS-15',
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
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20240213_periodsAndPaperFigs/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('allplot',),  # See docstring for valid inputs
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line
# ic(np.diff(setDict["calcIntvl"]["LastMillennium"])) # Check to ensure correct time periods

# Make images
# for pet in (0,1,2,3,4,5,6,7,8,9): # For ensemble member figures (Supplementary Fig. 3, 4, 5, 6)
for pet in ('mean',): # For ensemble mean figures
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
                # scnList, dataDict, setDict, outDict)

ic('Completed! :D')
